#!/usr/bin/env python3
"""Run both REGENIE steps on deterministic synthetic data."""

import argparse
import gzip
import logging
import math
import random
import shutil
import subprocess
import sys
from pathlib import Path

LOG = logging.getLogger("smoke")
RUNNER = Path(__file__).resolve().with_name("run_regenie.py")


def run(command, work):
    LOG.info("Run %s in %s", Path(str(command[0])).name, work.name)
    subprocess.run([str(value) for value in command], cwd=work, check=True)


def generate_data(work):
    rng = random.Random(20261013)
    samples = [f"sample_{i:04d}" for i in range(240)]
    genotypes = []
    vcf = work / "synthetic.vcf"
    with vcf.open("w") as handle:
        handle.write('##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=2>\n')
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for i in range(200):
            frequency = 0.1 + (i % 8) * 0.05
            dosage = [int(rng.random() < frequency) + int(rng.random() < frequency) for _ in samples]
            genotypes.append(dosage)
            calls = [("0/0", "0/1", "1/1")[value] for value in dosage]
            handle.write(f"{1 + i // 100}\t{1000 + (i % 100) * 10}\tv{i + 1}\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")
    with (work / "phenotypes.txt").open("w") as pheno, (work / "covariates.txt").open("w") as covar:
        pheno.write("FID IID qt1 qt2 bt01 bt12\n")
        covar.write("FID IID age batch\n")
        for i, sample in enumerate(samples):
            age = rng.uniform(25, 75)
            batch = i % 3
            qt1 = 0.2 * genotypes[0][i] + 0.02 * age + rng.gauss(0, 1)
            qt2 = 0.2 * genotypes[105][i] + rng.gauss(0, 1)
            probability = 1 / (1 + math.exp(-(-0.5 + 0.2 * genotypes[5][i] + 0.005 * age)))
            bt = int(rng.random() < probability)
            pheno.write(f"{sample} {sample} {qt1:.8f} {qt2:.8f} {bt} {bt + 1}\n")
            covar.write(f"{sample} {sample} {age:.8f} {batch}\n")
    (work / "variants,common.txt").write_text("".join(f"v{i + 1}\n" for i in range(200)))
    (work / "keep,samples.txt").write_text("".join(f"{sample} {sample}\n" for sample in samples))
    (work / "remove,samples.txt").write_text("".join(f"{sample} {sample}\n" for sample in samples[-5:]))
    run(["plink2", "--vcf", vcf, "--double-id", "--make-pgen", "--out", work / "whole", "--threads", "2"], work)
    for chrom in (1, 2):
        run(["plink2", "--pfile", work / "whole", "--chr", chrom, "--make-pgen", "--out", work / f"chr{chrom}", "--threads", "2"], work)


def common_arguments(data, prefix, trait_type, phenotype_columns, encoding):
    arguments = []
    for suffix in ("pgen", "pvar", "psam"):
        arguments.extend(["--" + suffix, str(data / f"{prefix}.{suffix}")])
    arguments.extend([
        "--phenotypes", str(data / "phenotypes.txt"), "--phenotype-columns", phenotype_columns,
        "--trait-type", trait_type, "--binary-encoding", encoding,
        "--covariates", str(data / "covariates.txt"), "--covariate-columns", "age,batch",
        "--categorical-covariates", "batch", "--keep", str(data / "keep,samples.txt"),
        "--remove", str(data / "remove,samples.txt"), "--threads", "2", "--block-size", "50",
    ])
    return arguments


def check_associations(directory, traits, chrom):
    for trait in traits.split(","):
        path = directory / "results" / f"assoc_{trait}.regenie.gz"
        with gzip.open(path, "rt") as handle:
            header = handle.readline().split()
            rows = [line.split() for line in handle if line.strip()]
        if not rows or not {"CHROM", "ID", "LOG10P"}.issubset(header):
            raise ValueError(f"Missing association rows or required columns: {path}")
        chromosome_index = header.index("CHROM")
        p_index = header.index("LOG10P")
        if any(row[chromosome_index] != str(chrom) for row in rows):
            raise ValueError(f"Association output contains the wrong chromosome: {path}")
        if not any(row[p_index] != "NA" and math.isfinite(float(row[p_index])) for row in rows):
            raise ValueError(f"No finite association statistics: {path}")


def test_model(data, work, trait_type, traits, encoding):
    fit = work / "fit"
    fit.mkdir(parents=True)
    extra = ["--rint"] if trait_type == "quantitative" else []
    run([sys.executable, RUNNER, "step1", *common_arguments(data, "whole", trait_type, traits, encoding),
         "--extract", str(data / "variants,common.txt"), "--lowmem", "--cv-folds", "3", *extra], fit)
    # Copy files into a separate task directory. The old paths become cloud
    # metadata. Step 2 must use the explicitly localized prediction files.
    localized = work / "localized"
    localized.mkdir()
    predictions = sorted((fit / "results").glob("fit_*.loco"))
    if len(predictions) != len(traits.split(",")):
        raise ValueError("Unexpected number of LOCO prediction files")
    for path in predictions:
        shutil.copyfile(path, localized / path.name)
    file_list = localized / "files.txt"
    file_list.write_text("".join(str(localized / path.name) + "\n" for path in predictions))
    rows = (fit / "results/fit_pred.list").read_text().splitlines()
    pred_list = localized / "pred.list"
    pred_list.write_text("".join(row.split()[0] + " gs://synthetic-bucket/" + Path(row.split()[1]).name + "\n" for row in rows))
    for chrom in (1, 2):
        association = work / f"test_chr{chrom}"
        association.mkdir()
        run([sys.executable, RUNNER, "step2", *common_arguments(data, f"chr{chrom}", trait_type, traits, encoding),
             "--prediction-list", str(pred_list), "--prediction-files", str(file_list), "--min-mac", "5",
             "--firth-threshold", "1", *extra], association)
        check_associations(association, traits, chrom)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    data = work / "data"
    data.mkdir()
    generate_data(data)
    for trait_type, traits, encoding in (("quantitative", "qt1,qt2", "01"), ("binary", "bt01", "01"), ("binary", "bt12", "12")):
        test_model(data, work / f"{trait_type}_{encoding}", trait_type, traits, encoding)
    LOG.info("Smoke test passed: two quantitative traits, both binary encodings, and two chromosome shards")


if __name__ == "__main__":
    main()

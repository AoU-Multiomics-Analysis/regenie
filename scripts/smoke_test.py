#!/usr/bin/env python3
"""Run both REGENIE steps on deterministic synthetic data."""

import argparse
import gzip
import json
import logging
import math
import random
import shutil
import subprocess
import sys
from pathlib import Path

LOG = logging.getLogger("smoke")
RUNNER = Path(__file__).resolve().with_name("run_regenie.py")
SPLITTER = Path(__file__).resolve().with_name("split_genotypes.py")
BATCHER = Path(__file__).resolve().with_name("prepare_batches.py")
COLLECTOR = Path(__file__).resolve().with_name("collect_predictions.py")


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
        pheno.write("FID IID qt1 qt2 bt01 bt12 bt02\n")
        covar.write("FID IID age batch\n")
        for i, sample in enumerate(samples):
            age = rng.uniform(25, 75)
            batch = i % 3
            qt1 = 0.2 * genotypes[0][i] + 0.02 * age + rng.gauss(0, 1)
            qt2 = 0.2 * genotypes[105][i] + rng.gauss(0, 1)
            probability = 1 / (1 + math.exp(-(-0.5 + 0.2 * genotypes[5][i] + 0.005 * age)))
            bt = int(rng.random() < probability)
            bt2 = int(rng.random() < 0.4)
            pheno.write(f"0 {sample} {qt1:.8f} {qt2:.8f} {bt} {bt + 1} {bt2}\n")
            covar.write(f"0 {sample} {age:.8f} {batch}\n")
    (work / "variants,common.txt").write_text("".join(f"v{i + 1}\n" for i in range(200)))
    (work / "keep,samples.txt").write_text("".join(f"0 {sample}\n" for sample in samples))
    (work / "remove,samples.txt").write_text("".join(f"0 {sample}\n" for sample in samples[-5:]))
    run(["plink2", "--vcf", vcf, "--double-id", "--make-pgen", "--out", work / "whole", "--threads", "2"], work)
    # Use a standard IID-only PSAM with no SEX. The splitter and REGENIE model
    # must add local FID=0 and SEX=0 without changing this source file.
    psam = work / "whole.psam"
    rows = [line.split() for line in psam.read_text().splitlines()]
    if rows[0][:2] != ["#FID", "IID"]:
        raise ValueError("Unexpected PLINK synthetic PSAM header")
    psam.write_text("#IID\n" + "".join(row[1] + "\n" for row in rows[1:]))


def split_data(data):
    psam_before = (data / "whole.psam").read_bytes()
    arguments = [value for suffix in ("pgen", "pvar", "psam")
                 for value in ("--" + suffix, data / f"whole.{suffix}")]
    all_work = data / "split"
    all_work.mkdir()
    run([sys.executable, SPLITTER, *arguments, "--threads", "2", "--memory-mb", "1024"], all_work)
    manifest = json.loads((all_work / "chromosome_genotypes.json").read_text())
    if [entry["chromosome"] for entry in manifest] != ["1", "2"]:
        raise ValueError("Split manifest has unexpected chromosomes")
    source_variants = [line.split()[:5] for line in (data / "whole.pvar").read_text().splitlines() if not line.startswith("#")]
    split_variants = [line.split()[:5] for entry in manifest for line in Path(entry["pvar"]).read_text().splitlines() if not line.startswith("#")]
    if source_variants != split_variants or (data / "whole.psam").read_bytes() != psam_before:
        raise ValueError("Splitting changed source samples, variant IDs, order, or alleles")
    # Export hard calls from the source and both outputs. Compare every value
    # and sample order, so the test also checks the split PGEN content.
    matrices = []
    for prefix in [data / "whole", *(Path(entry["pgen"]).with_suffix("") for entry in manifest)]:
        out = prefix.parent / (prefix.name + "_dosages")
        run(["plink2", "--pfile", prefix, "--export", "A", "--out", out, "--threads", "2", "--memory", "1024"], data)
        matrices.append([line.split() for line in Path(str(out) + ".raw").read_text().splitlines()])
    if matrices[0][0][6:] != matrices[1][0][6:] + matrices[2][0][6:]:
        raise ValueError("Split genotype columns differ from the source")
    for whole, first, second in zip(*[matrix[1:] for matrix in matrices], strict=True):
        if whole[1] != first[1] or whole[1] != second[1] or whole[6:] != first[6:] + second[6:]:
            raise ValueError("Splitting changed genotype values or sample order")
    for name, variants, expected_chromosomes in (
        ("selected", "v2\nv102\n", ["1", "2"]), ("one_chromosome", "v102\n", ["2"]),
    ):
        work = data / name
        work.mkdir()
        selected = work / "selected.txt"
        selected.write_text(variants)
        run([sys.executable, SPLITTER, *arguments, "--extract", selected, "--threads", "2", "--memory-mb", "1024"], work)
        entries = json.loads((work / "chromosome_genotypes.json").read_text())
        ids = [line.split()[2] for entry in entries for line in Path(entry["pvar"]).read_text().splitlines() if not line.startswith("#")]
        if ids != variants.split() or [entry["chromosome"] for entry in entries] != expected_chromosomes:
            raise ValueError("Split variant selection is incorrect")


def common_arguments(data, prefix, trait_type, phenotype_columns, encoding, phenotype_file=None):
    arguments = []
    for suffix in ("pgen", "pvar", "psam"):
        arguments.extend(["--" + suffix, str(data / f"{prefix}.{suffix}")])
    arguments.extend([
        "--phenotypes", str(phenotype_file or data / "phenotypes.txt"),
        "--trait-type", trait_type, "--binary-encoding", encoding,
        "--covariates", str(data / "covariates.txt"), "--covariate-columns", "age,batch",
        "--categorical-covariates", "batch", "--keep", str(data / "keep,samples.txt"),
        "--remove", str(data / "remove,samples.txt"), "--threads", "2", "--block-size", "50",
    ])
    if phenotype_columns is not None:
        arguments.extend(["--phenotype-columns", phenotype_columns])
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


def test_model(data, work, trait_type, traits, encoding, all_phenotypes=False):
    fit = work / "fit"
    fit.mkdir(parents=True)
    extra = ["--rint"] if trait_type == "quantitative" else []
    phenotype_file = None
    selection = traits
    if all_phenotypes:
        # All columns in an automatic run must use the selected model/encoding.
        rows = [line.split() for line in (data / "phenotypes.txt").read_text().splitlines()]
        indices = [0, 1, *(rows[0].index(trait) for trait in traits.split(","))]
        phenotype_file = work / "all_phenotypes.txt"
        phenotype_file.write_text("".join(" ".join(row[i] for i in indices) + "\n" for row in rows))
        selection = None if trait_type == "quantitative" else ""
    run([sys.executable, RUNNER, "step1", *common_arguments(data, "whole", trait_type, selection, encoding, phenotype_file),
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
    reference = {}
    for chrom in (1, 2):
        association = work / f"test_chr{chrom}"
        association.mkdir()
        selection = ("" if chrom == 1 else "   ") if all_phenotypes else traits
        run([sys.executable, RUNNER, "step2", *common_arguments(data, f"split/chromosomes/chr{chrom}", trait_type, selection, encoding, phenotype_file),
             "--prediction-list", str(pred_list), "--prediction-files", str(file_list), "--min-mac", "5",
             "--firth-threshold", "0.99", *extra], association)
        check_associations(association, traits, chrom)
        reference[chrom] = association
    test_batches(data, work, trait_type, traits, encoding, phenotype_file, extra, fit, reference)


def compare_associations(reference, batched, traits):
    for trait in traits.split(","):
        tables = []
        for directory in (reference, batched):
            with gzip.open(directory / "results" / f"assoc_{trait}.regenie.gz", "rt") as handle:
                header = handle.readline().split()
                tables.append([dict(zip(header, line.split())) for line in handle if line.strip()])
        if len(tables[0]) != len(tables[1]):
            raise ValueError("Phenotype batching changed the association variant count")
        for first, second in zip(*tables):
            for key in ("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "N", "TEST"):
                if first[key] != second[key]:
                    raise ValueError(f"Phenotype batching changed association field {key}")
            for key in ("A1FREQ", "BETA", "SE", "CHISQ", "LOG10P"):
                if first[key] == second[key] == "NA":
                    continue
                if not math.isclose(float(first[key]), float(second[key]), rel_tol=1e-4, abs_tol=1e-6):
                    raise ValueError(f"Phenotype batching changed {trait} {key}: {first[key]} vs {second[key]}")


def test_batches(data, work, trait_type, traits, encoding, phenotype_file, extra, shared_fit, reference):
    for separate_models in ((False, True) if len(traits.split(",")) > 1 else (False,)):
        batch_work = work / f"batches_step1_{separate_models}"
        plan = batch_work / "plan"
        plan.mkdir(parents=True)
        selection = ["--phenotype-columns", traits] if phenotype_file is None else []
        run([sys.executable, BATCHER, "--phenotypes", phenotype_file or data / "phenotypes.txt",
             *selection, "--batch-size", "1", *(["--batch-step1"] if separate_models else [])], plan)
        models = json.loads((plan / "model_batches.json").read_text())
        selected_traits = (plan / "selected_phenotypes.txt").read_text().strip()
        expected_count = len(traits.split(",")) if separate_models else 1
        if len(models) != expected_count:
            raise ValueError("Batch planner produced the wrong model count")
        fits = []
        for model in models:
            if separate_models:
                fit = batch_work / f"fit_{model['batch_index']}"
                fit.mkdir()
                run([sys.executable, RUNNER, "step1", *common_arguments(data, "whole", trait_type, model["phenotype_columns"], encoding, phenotype_file),
                     "--extract", data / "variants,common.txt", "--lowmem", "--cv-folds", "3", *extra], fit)
            else:
                fit = shared_fit
            fits.append(fit)
        collector = batch_work / "collect"
        collector.mkdir()
        lists = collector / "prediction_lists.txt"
        files = collector / "prediction_files.txt"
        lists.write_text("".join(str(fit / "results/fit_pred.list") + "\n" for fit in fits))
        model_files = [sorted((fit / "results").glob("fit_*.loco")) for fit in fits]
        files.write_text("".join(str(path) + "\n" for group in model_files for path in group))
        counts = ",".join(str(len(group)) for group in model_files)
        run([sys.executable, COLLECTOR, "--prediction-lists", lists, "--prediction-files", files,
             "--prediction-file-counts", counts, "--phenotype-columns", selected_traits, "--batch-size", "1"], collector)
        batches = json.loads((collector / "results/prediction_batches.json").read_text())
        if [row["phenotype_columns"] for row in batches] != traits.split(","):
            raise ValueError("Prediction collector lost or reordered a phenotype batch")
        for batch in batches:
            for chrom in (1, 2):
                association = batch_work / f"test_chr{chrom}_batch{batch['batch_index']}"
                association.mkdir()
                localized_files = association / "files.txt"
                localized_files.write_text("".join(path + "\n" for path in batch["prediction_files"]))
                run([sys.executable, RUNNER, "step2", *common_arguments(data, f"split/chromosomes/chr{chrom}", trait_type, batch["phenotype_columns"], encoding, phenotype_file),
                     "--prediction-list", batch["prediction_list"], "--prediction-files", localized_files,
                     "--min-mac", "5", "--firth-threshold", "0.99", *extra], association)
                check_associations(association, batch["phenotype_columns"], chrom)
                compare_associations(reference[chrom], association, batch["phenotype_columns"])


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
    split_data(data)
    for trait_type, traits, encoding, automatic in (("quantitative", "qt1,qt2", "01", False),
                                                   ("quantitative", "qt1,qt2", "01", True),
                                                   ("binary", "bt01,bt02", "01", True), ("binary", "bt12", "12", False)):
        test_model(data, work / f"{trait_type}_{encoding}_{automatic}", trait_type, traits, encoding, automatic)
    LOG.info("Smoke test passed: chromosome splitting, shared and separate model batches, both trait types, and joint/batched association agreement")


if __name__ == "__main__":
    main()

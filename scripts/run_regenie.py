#!/usr/bin/env python3
"""Check localized inputs, stage a PGEN prefix, and run one REGENIE step."""

from __future__ import annotations

import argparse
import logging
import math
import os
import re
import subprocess
from pathlib import Path

LOG = logging.getLogger("regenie")


def local_file(value: str | Path, label: str) -> Path:
    text = str(value)
    if re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*://", text):
        raise ValueError(f"File localization error: {label} is an unresolved cloud URI: {text}")
    if not text or any(char in text for char in "\r\n\t"):
        raise ValueError(f"File localization error: {label} has an invalid local path")
    path = Path(text).resolve()
    try:
        with path.open("rb") as handle:
            if not handle.read(1):
                raise ValueError(f"{label} is empty: {path}")
    except OSError as error:
        raise ValueError(f"File localization error: {label} is not readable: {path}") from error
    return path


def psam_adjustments(source: Path) -> tuple[bool, bool]:
    header = None
    seen = set()
    with source.open(encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##") and header is None:
                continue
            fields = line.split()
            if header is None:
                if not fields or (fields[0] != "#IID" and fields[:2] != ["#FID", "IID"]):
                    raise ValueError("PSAM header must start with #IID or #FID IID")
                header = fields
                missing_fid = fields[0] == "#IID"
                continue
            if len(fields) != len(header):
                raise ValueError(f"PSAM row {number} has a different number of columns from the header")
            pair = ("0", fields[0]) if missing_fid else (fields[0], fields[1])
            if pair in seen:
                raise ValueError(f"Duplicate sample ID pair in PSAM row {number}")
            seen.add(pair)
    if header is None or not seen:
        raise ValueError("PSAM must contain a header and at least one sample")
    return missing_fid, "SEX" not in header


def stage_genotypes(pgen: str | Path, pvar: str | Path, psam: str | Path, directory: Path) -> Path:
    sources = [local_file(value, suffix.upper()) for value, suffix in zip((pgen, pvar, psam), ("pgen", "pvar", "psam"))]
    with sources[1].open("rb") as handle:
        magic = handle.read(4)
        if magic[:2] == b"\x1f\x8b" or magic == b"\x28\xb5\x2f\xfd":
            raise ValueError("PVAR must be uncompressed. Decompress .pvar.gz or .pvar.zst before the run.")
    directory.mkdir(parents=True, exist_ok=True)
    missing_fid, missing_sex = psam_adjustments(sources[2])
    prefix = directory / "genotypes"
    for source, suffix in zip(sources, ("pgen", "pvar", "psam")):
        target = Path(str(prefix) + "." + suffix)
        if target.is_symlink():
            target.unlink()
        elif target.exists():
            raise ValueError(f"Refusing to replace an existing input file: {target}")
        if suffix == "psam" and (missing_fid or missing_sex):
            with source.open(encoding="utf-8-sig") as reader, target.open("w") as writer:
                for line in reader:
                    if line.startswith("##"):
                        writer.write(line)
                    else:
                        fields = line.split()
                        is_header = fields[0] in ("#FID", "#IID")
                        if missing_fid:
                            fields = ["#FID", "IID", *fields[1:]] if is_header else ["0", *fields]
                        if missing_sex:
                            fields.append("SEX" if is_header else "0")
                        writer.write("\t".join(fields) + "\n")
            LOG.info("Prepared task-local PSAM: added zero FID=%s, added unknown SEX=%s", missing_fid, missing_sex)
        else:
            target.symlink_to(source)
    return prefix


def columns(value: str) -> list[str]:
    names = value.split(",")
    if any(not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name) for name in names):
        raise ValueError("Column names must contain letters, digits, underscores, dots, or hyphens. Separate names with commas.")
    if len(set(names)) != len(names):
        raise ValueError("Duplicate column names are not allowed")
    return names


def sample_pairs(path: str) -> list[tuple[str, str]]:
    rows = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"Sample filter row {number} must contain FID IID with no header")
        rows.append((fields[0], fields[1]))
    if not rows or rows[0] == ("FID", "IID"):
        raise ValueError("Sample filter must contain FID IID pairs with no header")
    return rows


def stage_auxiliary(args: argparse.Namespace, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for name in ("phenotypes", "covariates", "keep", "remove", "extract"):
        value = getattr(args, name)
        if value is None:
            continue
        source = local_file(value, name)
        target = directory / f"{name}.txt"
        if target.is_symlink():
            target.unlink()
        elif target.exists():
            raise ValueError(f"Refusing to replace input file: {target}")
        target.symlink_to(source)
        setattr(args, name, str(target))
    if args.keep:
        keep = sample_pairs(args.keep)
    if args.remove:
        remove = set(sample_pairs(args.remove))
    if args.keep and args.remove:
        selected = [pair for pair in keep if pair not in remove]
        if not selected:
            raise ValueError("The keep/remove filters leave no selected samples")
        target = directory / "selected_samples.txt"
        target.write_text("".join(f"{fid} {iid}\n" for fid, iid in selected))
        args.keep = str(target)
        args.remove = None
        LOG.info("Combined sample filters: %d sample IDs retained", len(selected))


def verify_association_outputs(phenotype_columns: str, directory: Path) -> None:
    for trait in columns(phenotype_columns):
        path = directory / f"assoc_{trait}.regenie.gz"
        if not path.is_file():
            raise ValueError(f"Step 2 did not produce an association file for {trait}. Check sample overlap, trait values, and case counts.")
        local_file(path, f"association output for {trait}")


def prediction_rows(path: Path) -> list[tuple[str, str]]:
    rows = []
    traits = set()
    for line in local_file(path, "prediction list").read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or fields[0] in traits:
            raise ValueError("Prediction list must contain one trait and one path per row, with no duplicate traits")
        columns(fields[0])
        traits.add(fields[0])
        rows.append((fields[0], fields[1]))
    if not rows:
        raise ValueError("Prediction list has no traits")
    return rows


def rebind_predictions(prediction_list: Path, files: list[Path], directory: Path) -> Path:
    sources = {}
    for value in files:
        path = local_file(value, "LOCO prediction")
        if path.name in sources:
            raise ValueError(f"Duplicate prediction filename: {path.name}")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", path.name):
            raise ValueError(f"Invalid prediction filename: {path.name}")
        sources[path.name] = path
    rows = prediction_rows(prediction_list)
    directory.mkdir(parents=True, exist_ok=True)
    rebound = []
    for trait, old_path in rows:
        name = old_path.rsplit("/", 1)[-1]
        if name not in sources:
            raise ValueError(f"No explicitly localized prediction file for {trait}: {name}")
        target = directory / name
        if target.is_symlink():
            target.unlink()
        elif target.exists():
            raise ValueError(f"Refusing to replace prediction file: {target}")
        target.symlink_to(sources[name])
        rebound.append(f"{trait} {target.as_posix()}\n")
    output = directory / "predictions.list"
    output.write_text("".join(rebound))
    return output


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("step", choices=("step1", "step2"))
    for name in ("pgen", "pvar", "psam", "phenotypes", "trait-type"):
        result.add_argument("--" + name, required=True)
    result.add_argument("--phenotype-columns", help="Comma-separated trait names. Omit or leave blank to use all phenotype columns.")
    for name in ("covariates", "covariate-columns", "categorical-covariates", "keep", "remove", "extract", "prediction-list", "prediction-files"):
        result.add_argument("--" + name)
    result.add_argument("--binary-encoding", choices=("01", "12"), default="01")
    result.add_argument("--threads", type=int, default=4)
    result.add_argument("--block-size", type=int, default=1000)
    result.add_argument("--cv-folds", type=int, default=5)
    result.add_argument("--min-mac", type=float, default=5)
    result.add_argument("--min-case-count", type=int, default=10)
    result.add_argument("--firth-threshold", type=float, default=0.01)
    result.add_argument("--lowmem", action="store_true")
    result.add_argument("--rint", action="store_true")
    result.add_argument("--no-firth", action="store_true")
    return result


def validate_options(args: argparse.Namespace) -> None:
    if args.trait_type not in ("quantitative", "binary"):
        raise ValueError("Trait type must be quantitative or binary")
    if args.phenotype_columns and args.phenotype_columns.strip():
        columns(args.phenotype_columns)
    for value in (args.covariate_columns, args.categorical_covariates):
        if value:
            columns(value)
    if args.rint and args.trait_type != "quantitative":
        raise ValueError("Inverse normal transformation requires quantitative traits")
    if args.trait_type == "quantitative" and args.binary_encoding != "01":
        raise ValueError("Binary encoding 12 requires binary traits")
    if any(value < 1 for value in (args.threads, args.block_size, args.min_case_count)) or args.cv_folds < 2:
        raise ValueError("Threads, block size, and minimum case count must be positive. CV folds must be at least 2.")
    if not math.isfinite(args.min_mac) or args.min_mac < 0:
        raise ValueError("Minimum MAC must be a finite, nonnegative number")
    if not 0 < args.firth_threshold < 1:
        raise ValueError("Firth threshold must be greater than zero and less than one")
    if (args.covariate_columns or args.categorical_covariates) and not args.covariates:
        raise ValueError("Selected covariates require a covariate file")
    if args.step == "step1" and not args.extract:
        raise ValueError("Step 1 requires a selected common-variant list through --extract")
    if args.step == "step2" and not (args.prediction_list and args.prediction_files):
        raise ValueError("Step 2 requires a prediction list and explicitly localized prediction files")


def table_header(path: str, label: str) -> list[str]:
    with local_file(path, label).open() as handle:
        header = handle.readline().strip().split()
    if header[:2] != ["FID", "IID"] or len(header) < 3 or len(set(header)) != len(header):
        raise ValueError(f"{label} needs a unique header that starts with FID IID")
    return header[2:]


def validate_tables(args: argparse.Namespace) -> None:
    header = table_header(args.phenotypes, "phenotype file")
    if not args.phenotype_columns or not args.phenotype_columns.strip():
        args.phenotype_columns = ",".join(header)
        LOG.info("Selected all %d phenotype columns from the localized file header", len(header))
    missing = set(columns(args.phenotype_columns)) - set(header)
    if missing:
        raise ValueError(f"Unknown phenotype columns: {', '.join(sorted(missing))}")
    if args.covariates:
        header = table_header(args.covariates, "covariate file")
        selected = columns(args.covariate_columns) if args.covariate_columns else header
        categorical = columns(args.categorical_covariates) if args.categorical_covariates else []
        missing = (set(selected) | set(categorical)) - set(header)
        if missing or not set(categorical).issubset(selected):
            raise ValueError("Selected and categorical covariates must exist in the covariate file. Categorical covariates must also be selected.")


def regenie_command(args: argparse.Namespace, prefix: str | Path, prediction_list: str | Path | None = None) -> list[str]:
    command = ["regenie", "--step", args.step[-1], "--pgen", str(prefix), "--phenoFile", args.phenotypes,
               "--phenoColList", args.phenotype_columns, "--bsize", str(args.block_size), "--threads", str(args.threads),
               "--out", "results/fit" if args.step == "step1" else "results/assoc"]
    command.append("--qt" if args.trait_type == "quantitative" else "--bt")
    if args.trait_type == "binary":
        command.extend(["--minCaseCount", str(args.min_case_count)])
        if args.binary_encoding == "12":
            command.append("--cc12")
    for attribute, flag in (("covariates", "--covarFile"), ("covariate_columns", "--covarColList"),
                            ("categorical_covariates", "--catCovarList"), ("keep", "--keep"),
                            ("remove", "--remove"), ("extract", "--extract")):
        value = getattr(args, attribute)
        if value is not None:
            command.extend([flag, value])
    if args.rint:
        command.append("--apply-rint")
    if args.step == "step1":
        command.extend(["--cv", str(args.cv_folds), "--use-relative-path"])
        if args.lowmem:
            command.extend(["--lowmem", "--lowmem-prefix", "results/lowmem"])
    else:
        command.extend(["--pred", str(prediction_list), "--minMAC", str(args.min_mac), "--gz"])
        if args.trait_type == "binary" and not args.no_firth:
            command.extend(["--firth", "--approx", "--pThresh", str(args.firth_threshold)])
    return command


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parser().parse_args()
    try:
        LOG.info("Check options and localized input files for %s", args.step)
        validate_options(args)
        validate_tables(args)
        stage_auxiliary(args, Path("inputs"))
        prefix = stage_genotypes(args.pgen, args.pvar, args.psam, Path("inputs"))
        prediction_list = None
        if args.step == "step2":
            file_list = local_file(args.prediction_files, "localized prediction file list")
            paths = [Path(line) for line in file_list.read_text().splitlines() if line]
            prediction_list = rebind_predictions(Path(args.prediction_list), paths, Path("inputs/predictions"))
            if not set(columns(args.phenotype_columns)).issubset(trait for trait, _ in prediction_rows(prediction_list)):
                raise ValueError("Step 1 predictions are missing for a selected phenotype")
            LOG.info("Rebuilt prediction paths from %d localized files", len(paths))
        Path("results").mkdir(exist_ok=True)
        os.environ["OMP_NUM_THREADS"] = str(args.threads)
        os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
        LOG.info("Start REGENIE %s with %s traits", args.step, args.trait_type)
        subprocess.run(regenie_command(args, prefix, prediction_list), check=True)
        if args.step == "step1":
            rows = prediction_rows(Path("results/fit_pred.list"))
            if set(columns(args.phenotype_columns)) != {trait for trait, _ in rows}:
                raise ValueError("Step 1 did not produce predictions for every selected phenotype. Check trait values and case counts.")
            for _, path in rows:
                local_file(path, "step 1 LOCO prediction")
        else:
            verify_association_outputs(args.phenotype_columns, Path("results"))
        LOG.info("REGENIE %s completed", args.step)
        return 0
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        LOG.error("%s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

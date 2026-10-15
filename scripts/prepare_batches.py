#!/usr/bin/env python3
"""Validate phenotype selection and plan REGENIE phenotype batches."""

import argparse
import json
import logging
from pathlib import Path

from run_regenie import columns, table_header

LOG = logging.getLogger("regenie")


def partition(selection: str, batch_size: int) -> list[dict]:
    if batch_size < 0:
        raise ValueError("Phenotype batch size must be zero or positive")
    names = columns(selection)
    width = batch_size or len(names)
    return [{"batch_index": i // width, "phenotype_columns": ",".join(names[i:i + width])}
            for i in range(0, len(names), width)]


def phenotype_batches(path: str | Path, selection: str | None, batch_size: int) -> tuple[str, list[dict]]:
    header = table_header(str(path), "phenotype file")
    selected = selection if selection and selection.strip() else ",".join(header)
    batches = partition(selected, batch_size)
    missing = set(columns(selected)) - set(header)
    if missing:
        raise ValueError(f"Unknown phenotype columns: {', '.join(sorted(missing))}")
    return selected, batches


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--phenotype-columns")
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--batch-step1", action="store_true")
    args = parser.parse_args()
    try:
        LOG.info("Start phenotype batch preparation")
        selected, batches = phenotype_batches(args.phenotypes, args.phenotype_columns, args.batch_size)
        models = batches if args.batch_step1 else [{"batch_index": 0, "phenotype_columns": selected}]
        Path("selected_phenotypes.txt").write_text(selected + "\n")
        Path("phenotype_batches.json").write_text(json.dumps(batches, indent=2) + "\n")
        Path("model_batches.json").write_text(json.dumps(models, indent=2) + "\n")
        LOG.info("Phenotype batch preparation completed: %d traits, %d association batches, %d model fits",
                 len(columns(selected)), len(batches), len(models))
        return 0
    except (ValueError, OSError) as error:
        LOG.error("%s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

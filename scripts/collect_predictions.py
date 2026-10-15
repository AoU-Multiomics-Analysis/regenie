#!/usr/bin/env python3
"""Collect localized model outputs and group LOCO predictions by phenotype batch."""

import argparse
import json
import logging
from pathlib import Path
import shutil

from prepare_batches import partition
from run_regenie import columns, local_file, prediction_rows

LOG = logging.getLogger("regenie")


def collect_predictions(prediction_lists, prediction_files, counts, selection, batch_size, directory):
    batches = partition(selection, batch_size)
    if not prediction_lists or len(counts) != len(prediction_lists) or any(count < 1 for count in counts) or sum(counts) != len(prediction_files):
        raise ValueError("Prediction file counts must match each model and the localized File array")
    by_trait = {}
    start = 0
    for listing, count in zip(prediction_lists, counts):
        rows = prediction_rows(local_file(listing, "model prediction list"))
        files = [local_file(value, "LOCO prediction") for value in prediction_files[start:start + count]]
        start += count
        sources = {path.name: path for path in files}
        names = [old.rsplit("/", 1)[-1] for _, old in rows]
        if len(sources) != count or len(names) != count or set(names) != set(sources):
            raise ValueError("Model prediction list does not match its localized prediction files")
        for trait, old in rows:
            if trait in by_trait:
                raise ValueError(f"Duplicate model phenotype: {trait}")
            by_trait[trait] = sources[old.rsplit("/", 1)[-1]]
    selected = columns(selection)
    if set(by_trait) != set(selected):
        raise ValueError("Model predictions do not match all selected phenotypes")
    directory = Path(directory).resolve()
    prediction_dir = directory / "predictions"
    prediction_dir.mkdir(parents=True)
    destinations = {}
    for index, trait in enumerate(selected, 1):
        target = prediction_dir / f"trait_{index}.loco"
        shutil.copyfile(by_trait[trait], target)
        destinations[trait] = str(target)
    (directory / "all_pred.list").write_text("".join(f"{trait} {destinations[trait]}\n" for trait in selected))
    result = []
    for batch in batches:
        traits = columns(batch["phenotype_columns"])
        listing = directory / f"batch_{batch['batch_index']}_pred.list"
        listing.write_text("".join(f"{trait} {destinations[trait]}\n" for trait in traits))
        result.append({**batch, "prediction_list": str(listing), "prediction_files": [destinations[trait] for trait in traits]})
    # This output manifest contains typed File outputs for the next scatter.
    (directory / "prediction_batches.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-lists", required=True, help="Task-local newline list of localized prediction-list Files")
    parser.add_argument("--prediction-files", required=True, help="Task-local newline list of localized LOCO Files")
    parser.add_argument("--prediction-file-counts", required=True, help="Comma-separated LOCO counts in model order")
    parser.add_argument("--phenotype-columns", required=True)
    parser.add_argument("--batch-size", type=int, default=50)
    args = parser.parse_args()
    try:
        LOG.info("Start prediction collection")
        lists = local_file(args.prediction_lists, "localized prediction-list file list").read_text().splitlines()
        files = local_file(args.prediction_files, "localized LOCO file list").read_text().splitlines()
        counts = [int(value) for value in args.prediction_file_counts.split(",")]
        result = collect_predictions(lists, files, counts, args.phenotype_columns, args.batch_size, Path("results"))
        LOG.info("Prediction collection completed: %d model fits, %d traits, %d batches", len(lists), len(files), len(result))
        return 0
    except (ValueError, OSError) as error:
        LOG.error("%s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

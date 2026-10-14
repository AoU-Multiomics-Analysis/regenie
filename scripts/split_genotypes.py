#!/usr/bin/env python3
"""Split one localized PLINK 2 file set into chromosome file sets."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
from pathlib import Path

from run_regenie import local_file, stage_genotypes

LOG = logging.getLogger("regenie")
CHROMOSOMES = {str(value): value for value in range(1, 23)} | {"X": 23, "Y": 24, "XY": 25, "MT": 26}


def chromosome_name(value: str) -> str:
    label = value.upper().removeprefix("CHR")
    label = {"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"}.get(label, label)
    if label not in CHROMOSOMES:
        raise ValueError(f"Unsupported human chromosome in PVAR: {value}")
    return label


def chromosome_counts(path: Path, selected: set[str] | None = None) -> list[tuple[str, int]]:
    counts = {}
    header = None
    with path.open(encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##") and header is None:
                continue
            fields = line.split()
            if header is None:
                if fields[:5] != ["#CHROM", "POS", "ID", "REF", "ALT"]:
                    raise ValueError("PVAR header must start with #CHROM POS ID REF ALT")
                header = fields
                continue
            if len(fields) < 5:
                raise ValueError(f"PVAR row {number} has fewer than five columns")
            if selected is None or fields[2] in selected:
                chromosome = chromosome_name(fields[0])
                counts[chromosome] = counts.get(chromosome, 0) + 1
    if header is None:
        raise ValueError("PVAR has no header")
    return sorted(counts.items(), key=lambda item: CHROMOSOMES[item[0]])


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    for suffix in ("pgen", "pvar", "psam"):
        result.add_argument("--" + suffix, required=True)
    result.add_argument("--extract")
    result.add_argument("--threads", type=int, default=4)
    result.add_argument("--memory-mb", type=int, default=13107)
    return result


def plink_command(args: argparse.Namespace, prefix: Path, chromosome: str) -> list[str]:
    command = ["plink2", "--pfile", str(prefix), "--chr", chromosome, "--make-pgen",
               "--out", f"chromosomes/chr{chromosome}", "--threads", str(args.threads),
               "--memory", str(args.memory_mb)]
    if args.extract:
        command.extend(["--extract", args.extract])
    return command


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parser().parse_args()
    try:
        if args.threads < 1 or args.memory_mb < 640:
            raise ValueError("Split task needs at least one thread and 640 MB of PLINK memory")
        prefix = stage_genotypes(args.pgen, args.pvar, args.psam, Path("inputs"))
        selected = None
        if args.extract:
            source = local_file(args.extract, "association variant IDs")
            selected = set()
            for number, line in enumerate(source.read_text().splitlines(), 1):
                fields = line.split()
                if len(fields) != 1:
                    raise ValueError(f"Variant filter row {number} must contain one ID with no header")
                selected.add(fields[0])
            target = Path("inputs/extract.txt")
            target.symlink_to(source)
            args.extract = str(target)
        counts = chromosome_counts(Path(str(prefix) + ".pvar"), selected)
        if not counts:
            raise ValueError("No association variants remain in this genotype file set. Check association_variant_ids.")
        LOG.info("Start chromosome split: %d variants across %d chromosomes", sum(n for _, n in counts), len(counts))
        Path("chromosomes").mkdir()
        outputs = []
        for chromosome, count in counts:
            LOG.info("Split chromosome %s: %d selected variants", chromosome, count)
            subprocess.run(plink_command(args, prefix, chromosome), check=True)
            files = {suffix: str(local_file(f"chromosomes/chr{chromosome}.{suffix}", f"chromosome {chromosome} {suffix}"))
                     for suffix in ("pgen", "pvar", "psam")}
            local_file(f"chromosomes/chr{chromosome}.log", "PLINK log")
            if chromosome_counts(Path(files["pvar"])) != [(chromosome, count)]:
                raise ValueError(f"Split output for chromosome {chromosome} has unexpected variants")
            outputs.append({"chromosome": chromosome, **files})
        # This is an output manifest. WDL reads its entries as typed File outputs.
        Path("chromosome_genotypes.json").write_text(json.dumps(outputs, indent=2) + "\n")
        LOG.info("Chromosome split completed: %d file sets", len(outputs))
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        LOG.error("%s", error)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()

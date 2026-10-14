#!/usr/bin/env python3
"""Write a selected PLINK 2 file set before a REGENIE task."""

from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path

from run_regenie import local_file, stage_auxiliary, stage_genotypes

LOG = logging.getLogger("subset")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    for name in ("pgen", "pvar", "psam"):
        result.add_argument("--" + name, required=True)
    for name in ("extract", "keep", "remove"):
        result.add_argument("--" + name)
    result.add_argument("--threads", type=int, default=4)
    result.add_argument("--memory-mb", type=int, default=13107)
    return result


def plink_command(args: argparse.Namespace, prefix: Path) -> list[str]:
    command = ["plink2", "--pfile", str(prefix), "--make-pgen", "--out", "results/subset",
               "--threads", str(args.threads), "--memory", str(args.memory_mb)]
    for name in ("extract", "keep", "remove"):
        value = getattr(args, name)
        if value is not None:
            command.extend(["--" + name, value])
    return command


def check_outputs(prefix: Path) -> tuple[int, int]:
    paths = {suffix: local_file(str(prefix) + "." + suffix, "PLINK " + suffix.upper())
             for suffix in ("pgen", "pvar", "psam")}
    counts = []
    for suffix in ("pvar", "psam"):
        with paths[suffix].open() as handle:
            counts.append(sum(1 for line in handle if line.strip() and not line.startswith("#")))
    if not all(counts):
        raise ValueError("PLINK subset has no variants or samples")
    return counts[0], counts[1]


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parser().parse_args()
    try:
        if args.threads < 1 or args.memory_mb < 640:
            raise ValueError("Threads must be positive. PLINK memory must be at least 640 MB.")
        LOG.info("Check localized genotype and filter files")
        # Reuse the REGENIE staging checks and keep-minus-remove semantics.
        args.phenotypes = None
        args.covariates = None
        stage_auxiliary(args, Path("inputs"))
        prefix = stage_genotypes(args.pgen, args.pvar, args.psam, Path("inputs"))
        Path("results").mkdir(exist_ok=True)
        LOG.info("Start PLINK subset: variant filter=%s, keep filter=%s, remove filter=%s",
                 bool(args.extract), bool(args.keep), bool(args.remove))
        subprocess.run(plink_command(args, prefix), check=True)
        variants, samples = check_outputs(Path("results/subset"))
        LOG.info("PLINK subset completed: %d variants, %d samples", variants, samples)
        return 0
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        LOG.error("%s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Supply PLINK-shaped outputs for File mapping tests; not a genotype test."""

import shutil
import sys
from pathlib import Path

args = sys.argv[1:]
prefix = args[args.index("--pfile") + 1]
out = args[args.index("--out") + 1]
chromosome = args[args.index("--chr") + 1]
selected = None
if "--extract" in args:
    selected = set(Path(args[args.index("--extract") + 1]).read_text().split())
rows = []
for line in Path(prefix + ".pvar").read_text().splitlines():
    if line.startswith("#"):
        rows.append(line)
    else:
        fields = line.split()
        if fields[0].upper().removeprefix("CHR") == chromosome and (selected is None or fields[2] in selected):
            rows.append(line)
Path(out + ".pvar").write_text("\n".join(rows) + "\n")
for suffix in ("pgen", "psam"):
    shutil.copyfile(prefix + "." + suffix, out + "." + suffix)
Path(out + ".log").write_text("FILE_MAPPING_TEST_ONLY\n")

"""Select one ordered cohort and write analysis tables."""
import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

from .contracts import LOG, add_scientific_arguments, options_from_args, require_local_file, setup_logging

MISSING = {"", "NA", "NaN", "nan", "."}

@dataclass(frozen=True)
class PreparedFiles:
    sample_map: Path
    keep_samples: Path
    phenotypes: Path
    covariates: Path
    batches: Path
    report: Path

def sample_token(fid, iid):
    raw = json.dumps([fid, iid], ensure_ascii=False, separators=(",", ":")).encode()
    return "s" + hashlib.sha256(raw).hexdigest()

def read_ids(path):
    table = pd.read_csv(require_local_file(path), sep=r"\s+", dtype=str, keep_default_na=False)
    table = table.rename(columns={c: c.lstrip("#") for c in table.columns})
    if "IID" not in table:
        raise ValueError("Input table must contain IID")
    if "FID" not in table:
        table.insert(0, "FID", "0")
    if table[["FID", "IID"]].eq("").any().any() or table.duplicated(["FID", "IID"]).any():
        raise ValueError("Empty or duplicate sample pair")
    return table.set_index(["FID", "IID"], drop=False)

def read_filter(path):
    if path is None:
        return None
    rows = [line.split() for line in require_local_file(path).read_text().splitlines() if line.strip()]
    if any(len(row) != 2 for row in rows):
        raise ValueError("Sample filter must have FID IID pairs")
    return set(map(tuple, rows))

def select_columns(table, requested, label):
    selected = list(requested) if requested else [c for c in table if c not in ("FID", "IID")]
    if len(selected) != len(set(selected)) or any(c not in table or not re.fullmatch(r"[A-Za-z0-9_.-]+", c) for c in selected):
        raise ValueError(f"Invalid {label} columns")
    return selected

def prepare_cohort(psam, phenotypes, covariates, keep, remove, options, out_dir):
    LOG.info("Validate phenotype and sample inputs")
    samples, raw = read_ids(psam), read_ids(phenotypes)
    traits = select_columns(raw, options.phenotype_columns, "phenotype")
    if not traits:
        raise ValueError("No phenotype columns")
    keep_ids, remove_ids = read_filter(keep), read_filter(remove) or set()
    selected = [pair for pair in samples.index if pair in raw.index and (keep_ids is None or pair in keep_ids) and pair not in remove_ids]
    if not selected:
        raise ValueError("Empty analysis cohort")
    pheno = raw.loc[selected, traits].replace(list(MISSING), np.nan).apply(pd.to_numeric, errors="raise")
    cov_raw = read_ids(covariates) if covariates else None
    columns = select_columns(cov_raw, options.covariate_columns, "covariate") if cov_raw is not None else []
    if options.categorical_covariates and not set(options.categorical_covariates) <= set(columns):
        raise ValueError("Categorical covariate is not selected")
    cov = cov_raw.reindex(pheno.index)[columns].replace(list(MISSING), np.nan) if cov_raw is not None else pd.DataFrame(index=pheno.index)
    complete = pheno.notna().all(axis=1) & cov.notna().all(axis=1)
    removed = int((~complete).sum())
    if removed and options.missing_data_policy == "error":
        raise ValueError(f"Selected cohort has missing values in {removed} samples")
    pheno, cov = pheno.loc[complete], cov.loc[complete]
    if len(pheno) < options.sample_block_count:
        raise ValueError("Too few samples for the sample blocks")
    pheno = pheno.iloc[np.random.default_rng(options.seed).permutation(len(pheno))]
    cov = cov.reindex(pheno.index)
    report = {"sample_count": len(pheno), "removed_incomplete": removed,
              "model_samples_without_phenotype_row": len(samples.index.difference(raw.index)), "traits": {}, "categories": {}}
    for name in traits:
        values = pheno[name].to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite phenotype value")
        if len(set(values)) < 2:
            raise ValueError("Selected phenotype is constant")
        if options.trait_type == "binary":
            allowed = (0, 1) if options.binary_encoding == "01" else (1, 2)
            if not set(values) <= set(allowed):
                raise ValueError("Invalid binary phenotype value")
            values = values - allowed[0]
            cases = int(values.sum())
            controls = len(values) - cases
            if cases < options.min_cases or controls < options.min_controls:
                raise ValueError("Insufficient binary cases or controls")
            report["traits"][name] = {"cases": cases, "controls": controls}
            pheno[name] = values
        else:
            if options.inverse_normal_transform:
                pheno[name] = norm.ppf((rankdata(values, method="average") - .5) / len(values))
            report["traits"][name] = {"samples": len(values)}
    encoded = pd.DataFrame(index=pheno.index)
    for name in columns:
        if name in options.categorical_covariates:
            levels = sorted(cov[name].unique())
            report["categories"][name] = {"reference": levels[0], "levels": levels}
            for level in levels[1:]:
                encoded[f"cov_{len(encoded.columns)}"] = (cov[name] == level).astype(float)
        else:
            encoded[f"cov_{len(encoded.columns)}"] = pd.to_numeric(cov[name], errors="raise")
    design = np.column_stack([np.ones(len(pheno)), encoded.to_numpy(float)])
    if not np.isfinite(design).all() or np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Covariate design is not finite or has deficient rank")
    pairs = list(pheno.index)
    identifiers = [sample_token(fid, iid) for fid, iid in pairs]
    mapping = pd.DataFrame({"sample_id": identifiers, "FID": [p[0] for p in pairs], "IID": [p[1] for p in pairs]})
    pheno.index = encoded.index = pd.Index(identifiers, name="sample_id")
    batch_size = options.batch_size or len(traits)
    batches = pd.DataFrame([{"batch_index": i // batch_size, "trait_index": i, "trait_name": trait} for i, trait in enumerate(traits)])
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = PreparedFiles(*(out_dir / n for n in ("sample_map.tsv", "keep_samples.tsv", "phenotypes.tsv", "covariates.tsv", "phenotype_batches.tsv", "preparation.json")))
    mapping.to_csv(files.sample_map, sep="\t", index=False)
    mapping[["FID", "IID"]].to_csv(files.keep_samples, sep="\t", index=False, header=False)
    pheno.to_csv(files.phenotypes, sep="\t")
    encoded.to_csv(files.covariates, sep="\t")
    batches.to_csv(files.batches, sep="\t", index=False)
    files.report.write_text(json.dumps(report, indent=2) + "\n")
    LOG.info("Prepared %d samples and %d traits", len(pheno), len(traits))
    return files

def main(argv=None):
    setup_logging()
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("psam", "phenotypes", "out-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("covariates", "keep", "remove"):
        parser.add_argument("--" + name, type=Path)
    add_scientific_arguments(parser)
    args = parser.parse_args(argv)
    prepare_cohort(args.psam, args.phenotypes, args.covariates, args.keep, args.remove, options_from_args(args), args.out_dir)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

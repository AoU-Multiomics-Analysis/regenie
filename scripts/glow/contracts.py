"""Shared typed options and command-line interfaces."""
import argparse
import logging
from dataclasses import asdict, dataclass, fields
from pathlib import Path

LOG = logging.getLogger("regenie.spark")
GLOW_COMMIT = "72637ad17ff267ca6bbc63cfc77c78db619f38b5"

@dataclass(frozen=True)
class ScientificOptions:
    trait_type: str
    phenotype_columns: tuple[str, ...] | None = None
    covariate_columns: tuple[str, ...] | None = None
    categorical_covariates: tuple[str, ...] = ()
    binary_encoding: str = "01"
    missing_data_policy: str = "error"
    batch_size: int = 50
    variants_per_block: int = 1000
    sample_block_count: int = 10
    min_mac: float = 5
    min_cases: int = 10
    min_controls: int = 10
    firth_threshold: float = 0.01
    seed: int = 1
    inverse_normal_transform: bool = False

    def __post_init__(self):
        if self.trait_type not in ("quantitative", "binary"):
            raise ValueError("trait_type must be quantitative or binary")
        if self.binary_encoding not in ("01", "12") or self.missing_data_policy not in ("error", "complete_case"):
            raise ValueError("Invalid encoding or missing-data policy")
        if self.batch_size < 0 or min(self.variants_per_block, self.min_cases, self.min_controls) < 1 or self.sample_block_count < 2:
            raise ValueError("Invalid batch size, block count, or sample counts")
        import math
        if not math.isfinite(self.min_mac) or self.min_mac < 0 or not 0 <= self.firth_threshold <= 1:
            raise ValueError("Invalid MAC or Firth threshold")
        if self.trait_type == "binary" and self.inverse_normal_transform:
            raise ValueError("Binary traits cannot use inverse normal transformation")

def require_local_file(path: Path) -> Path:
    path = Path(path)
    if str(path).startswith(("gs:", "s3:", "http:", "https:")):
        raise ValueError("File localization error: unresolved cloud URI")
    if "\n" in str(path) or "\r" in str(path):
        raise ValueError("Newline in input path")
    if not path.is_file():
        raise ValueError(f"Input file is not readable: {path}")
    with path.open("rb"):
        pass
    return path

def add_scientific_arguments(parser):
    parser.add_argument("--trait-type", choices=("quantitative", "binary"), required=True)
    for name in ("phenotype-columns", "covariate-columns", "categorical-covariates"):
        parser.add_argument("--" + name)
    for name in ("binary-encoding", "missing-data-policy"):
        parser.add_argument("--" + name, default=getattr(ScientificOptions("quantitative"), name.replace("-", "_")))
    defaults = ScientificOptions("quantitative")
    for name in ("batch_size", "variants_per_block", "sample_block_count", "min_cases", "min_controls", "seed"):
        parser.add_argument("--" + name.replace("_", "-"), type=int, default=getattr(defaults, name))
    for name in ("min_mac", "firth_threshold"):
        parser.add_argument("--" + name.replace("_", "-"), type=float, default=getattr(defaults, name))
    parser.add_argument("--inverse-normal-transform", action="store_true")

def options_from_args(args):
    values = {f.name: getattr(args, f.name) for f in fields(ScientificOptions)}
    for key in ("phenotype_columns", "covariate_columns", "categorical_covariates"):
        raw = values[key]
        values[key] = tuple(raw.split(",")) if raw and raw.strip() else (() if key == "categorical_covariates" else None)
    return ScientificOptions(**values)

def scientific_arguments(options):
    result = []
    for key, value in asdict(options).items():
        if value is None or value == ():
            continue
        flag = "--" + key.replace("_", "-")
        if isinstance(value, bool):
            if value:
                result.append(flag)
        else:
            result.extend([flag, ",".join(value) if isinstance(value, tuple) else str(value)])
    return result

def setup_logging():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

@dataclass(frozen=True)
class InputArtifact:
    name: str
    path: Path

@dataclass(frozen=True)
class StagedArtifact:
    name: str
    uri: str
    sha256: str
    size_bytes: int
    generation: int

@dataclass(frozen=True)
class OutputArtifact:
    kind: str
    uri: str
    relative_path: str
    sha256: str
    size_bytes: int

# Glow reshapes NumPy arrays through pandas when it creates Spark tables.
SPARK_ARROW_PROPERTIES = {
    "spark.sql.execution.arrow.pyspark.enabled": "true",
    "spark.sql.execution.arrow.pyspark.fallback.enabled": "false",
}

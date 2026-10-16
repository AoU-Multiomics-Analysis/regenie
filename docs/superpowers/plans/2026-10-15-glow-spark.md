# Glow Spark Workflow Implementation Plan

> **October 16, 2026 update:** The user selected Spark on one Terra VM, like
> MTtoVCF, as the default. `GlowRegenie` now uses `local[analysis_cpu]` and normal
> Cromwell File localization. It has no Dataproc setup inputs. The original
> Dataproc implementation described below remains a separate optional workflow,
> `GlowRegenieDataproc`. The original implementation plan is complete; do not
> restart it. See [the current Terra guide](../../glow-terra.md) and
> [the Dataproc guide](../../glow-dataproc.md) for current use and validation status.


> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a separate Terra WDL that runs Glow model fitting and association tests on Dataproc Serverless.

**Architecture:** WDL tasks prepare a common sample set and convert PGEN inputs to BGEN. A controller stages explicit File inputs, submits one Spark batch, and collects verified outputs. The Spark application processes phenotype batches in sequence and distributes each batch's genotype calculations across workers.

**Tech Stack:** WDL 1.0, Python 3.11.13, PLINK 2 2.0.0a.6.9, Glow 2.0.3, Spark 3.5.1, Scala 2.12.19, Dataproc runtime 1.2, micromamba 2.3.3, GitHub Actions, unittest, miniwdl, Cromwell 87.

**Spec:** [Approved design](../specs/2026-10-15-glow-spark-design.md).

## Global Constraints

- Target a standard Terra workspace on Google Cloud and managed Cromwell.
- Use a separate WDL 1.0 workflow and image. Preserve the existing C++ workflow.
- Support biallelic, diploid variants on human autosomes 1–22.
- Preserve typed File inputs through command rendering. Create file lists within tasks after localization.
- Use safely quoted named CLI arguments. Do not use JSON argument wrappers.
- Reject workflow-scope file-writing functions and manifest-derived WDL File outputs.
- Use missing-data policy `error` by default; `complete_case` requires an explicit selection.
- Use 50 traits per batch, 1,000 variants per model block, 10 sample blocks, and minimum MAC 5.
- Require at least 10 cases and 10 controls per binary trait. Default Firth threshold: 0.01.
- Export BGEN 1.2 with 16-bit probabilities and reference allele first. Dosage tolerance: `1e-4`.
- Submit one remote batch per run. Default remote lifetime: 21,600 seconds; minimum polling interval: 30 seconds.
- Pin Glow source commit `72637ad17ff267ca6bbc63cfc77c78db619f38b5`, sbt 1.9.9, and Scala 2.12.19.
- Build and smoke-test images in GitHub Actions. Automatic image builds respond only to included script changes.
- Do not install Spark, PySpark, or Java in the production image. Supply them in a separate CI test target.
- Log task commands and remote state changes. Do not log credentials or phenotype values.
- Do not submit cloud jobs or claim Terra validation during implementation.
- Create an isolated worktree when execution starts, after checking existing task attachments.

## Review Focus

1. Missing FID, leading-zero IDs, and ambiguous ID pairs must retain exact sample identity. Task 1 owns these tests.
2. Source sets with identical basenames must retain distinct outputs and input identities. Tasks 2 and 5 own these tests.
3. Lost submission responses and controller retries must not launch duplicate remote jobs. Task 6 owns these tests.
4. Remote success without complete, matching outputs must fail collection. Task 5 owns these tests.
5. A stale or incompatible Glow artifact must fail CI before image publication. Task 8 owns these tests.

## File Structure

| Path | Responsibility |
| --- | --- |
| `scripts/glow/contracts.py` | Typed options, artifact records, and local path validation |
| `scripts/glow/prepare.py` | Sample selection, covariate design, and phenotype preparation |
| `scripts/glow/convert.py` | PGEN validation, local PSAM normalization, and BGEN export |
| `scripts/glow/genotypes.py` | Spark BGEN loading, dosage calculation, and variant checks |
| `scripts/glow/analysis.py` | Ridge models, chromosome offsets, and association tests |
| `scripts/glow/storage.py` | Staging identity, remote files, output manifests, and checksums |
| `scripts/glow/controller.py` | Dataproc submission, retries, monitoring, and result collection |
| `scripts/glow/spark_main.py` | Remote CLI, data persistence, phenotype batches, and output completion |
| `scripts/glow/smoke_test.py` | Deterministic synthetic data and Spark integration checks |
| `scripts/glow/__init__.py` | Package marker; no import of optional Spark or cloud dependencies |
| `workflows/glow_regenie.wdl` | Terra task graph and explicit output collection |
| `envs/glow/Dockerfile` | Source builder, production image, and CI test target |
| `envs/glow/environment.yml` | Exact production dependency versions |
| `envs/glow/build-environment.yml` | Exact source-build dependency versions |
| `.github/workflows/glow-image.yml` | Separate image build, smoke test, and publication |
| `tests/test_glow_*.py` | Focused unit, localization, cloud-state, and integration checks |
| `docs/glow-terra.md` | Setup, inputs, resource limits, recovery, and validation status |
| `examples/glow.*.inputs.json` | Separate quantitative and binary Terra input examples |

Keep the project package at `scripts.glow` and the upstream package at `glow`.
Install project code under `/opt/regenie-spark` and set its parent directory on `PYTHONPATH`.
Do not create a top-level `glow.py` module that hides the upstream library.

## Task 1: Prepare a stable analysis cohort

**Files:** Create `scripts/glow/__init__.py`, `scripts/glow/contracts.py`, `scripts/glow/prepare.py`, `tests/test_glow_inputs.py`.

**Interfaces:** Define frozen `ScientificOptions` with these typed fields:

```python
trait_type: str  # required: quantitative or binary
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
```

Define `PreparedFiles(sample_map: Path, keep_samples: Path, phenotypes: Path, covariates: Path, batches: Path, report: Path)`.
Provide `sample_token(fid: str, iid: str) -> str` and `require_local_file(path: Path) -> Path`.
Provide `prepare_cohort(psam: Path, phenotypes: Path, covariates: Path | None, keep: Path | None, remove: Path | None, options: ScientificOptions, out_dir: Path) -> PreparedFiles`.

- [ ] **Step 1: Write failing identity and cohort tests.** Include these assertions:

  ```python
  assert sample_token("a_b", "c") != sample_token("a", "b_c")
  assert sample_token("0", "001") != sample_token("0", "1")
  assert options.batch_size == 50
  assert options.missing_data_policy == "error"
  # test_complete_case_uses_all_selected_traits_before_batching
  assert batch_0_sample_ids == batch_1_sample_ids == retained_sample_ids
  ```

  Add named tests for PSAM without FID, duplicate pairs, keep/remove precedence, empty trait selection,
  invalid binary values in both encodings, case/control limits, categorical reference levels,
  constant traits, rank-deficient covariates, empty cohorts, unresolved `gs://` paths, and deterministic transformation.

- [ ] **Step 2: Run the failing tests.** Run `python3 -m unittest discover -s tests -p test_glow_inputs.py -v`. Confirm a missing implementation causes failure.

- [ ] **Step 3: Implement the interfaces and preparation CLI.** Hash a canonical UTF-8 JSON array `[FID, IID]` with SHA-256 for internal IDs; use prefix `s` and the full lowercase digest. Serialize with `ensure_ascii=False` and compact separators. Preserve the source strings in `sample_map.tsv`. Shuffle the retained model sample order with NumPy's seeded generator, then record that order explicitly. Use average ranks and `norm.ppf((rank - 0.5) / n)` for the optional quantitative transformation. Write `keep_samples.tsv`, `phenotypes.tsv`, `covariates.tsv`, `phenotype_batches.tsv`, and `preparation.json`.
  The prepared phenotype and covariate tables begin with `sample_id`; the batch table has `batch_index`, `trait_index`, and `trait_name`.
  An absent covariate input produces an index-only prepared covariate table.
  The CLI accepts `--psam`, `--phenotypes`, optional table/filter paths, individual scientific options, and `--out-dir`.

- [ ] **Step 4: Run the same tests.** Require all cases to pass. Check that input tables remain byte-identical and that report counts match retained rows.

- [ ] **Step 5: Commit only Task 1 files.** Commit message: `Add Glow cohort preparation and sample identity checks`.

## Task 2: Convert explicit PGEN inputs without losing dosages

**Files:** Create `scripts/glow/convert.py`, `tests/test_glow_conversion.py`, `tests/fake_glow_plink2.py`.

**Interfaces:** Define `ConvertedFiles(bgen: Path, sample: Path, variants: Path, report: Path, log: Path)`.
Provide `convert_pgen(pgen: Path, pvar: Path, psam: Path, sample_map: Path, variant_ids: Path | None, source_index: int, role: str, out_dir: Path, threads: int, memory_mb: int) -> ConvertedFiles`.
`role` is `model` or `association`. `variants.tsv` records source index, chromosome, position, variant ID, REF, and ALT.

- [ ] **Step 1: Write failing conversion tests.** Assert the real argument vector contains `--export bgen-1.2 bits=16 ref-first id-paste=iid`. Assert conversion preserves source indices and selected sample order. Add failures for duplicate variant IDs, multiallelic records, unsupported chromosomes, missing retained samples, mismatched triplets, and unreadable local files. Include two source sets with identical basenames.

- [ ] **Step 2: Run the failing tests.** Run `python3 -m unittest discover -s tests -p test_glow_conversion.py -v`.

- [ ] **Step 3: Implement the converter and CLI.** Check PVAR/PSAM metadata before calling PLINK. Make a local PSAM copy in the source sample order, with prepared internal IDs for retained samples and distinct temporary IDs for excluded samples. Use local keep and sample-order files to filter and reorder the PGEN. Invoke PLINK through an argument list. Verify the exported sample order and selected variant metadata, then write explicit outputs and logs. Do not modify source files or copy their paths into an argument JSON.

- [ ] **Step 4: Run conversion tests.** The fake executable proves argument handling and file safety only. Reserve real dosage-preservation evidence for Task 8's image test.

- [ ] **Step 5: Commit only Task 2 files.** Commit message: `Add dosage-preserving BGEN conversion for Glow`.

## Task 3: Load and validate Spark genotype arrays

**Files:** Create `scripts/glow/genotypes.py`, `tests/test_glow_genotypes.py`.

**Interfaces:** Provide `read_dosages(spark: SparkSession, bgen_uri: str, sample_uri: str, variants_uri: str, expected_sample_ids: list[str], source_index: int) -> DataFrame`.
Return `contigName`, `start`, `variant_id`, `referenceAllele`, `alternateAllele`, `source_index`, `values`, `observed_count`, `alternate_allele_count`, and `minor_allele_count`.
Provide `retain_usable_variants(frame: DataFrame, min_mac: float | None) -> tuple[DataFrame, dict[str, int]]`.

- [ ] **Step 1: Write failing dosage tests.** Include pure arithmetic tests and Spark cases:

  ```python
  assert dosage_from_probabilities([0.1, 0.6, 0.3]) == 1.2
  assert observed_mac([0.2, 1.8, None]) == 2.0
  # test_missing_dosage_is_imputed_after_observed_mac
  assert imputed_values == [0.2, 1.8, 1.0]
  ```

  Define these scalar helper functions in `genotypes.py` as reference implementations for the Spark expressions.
  Add Spark failures for phased or nondiploid records, wrong sample order on any variant, absent chromosomes,
  missing IDs, nonfinite values, invalid probability sums, and disagreement with the exported variant table.

- [ ] **Step 2: Run the failing non-Spark cases.** Run `python3 -m unittest discover -s tests -p test_glow_genotypes.py -v`. Spark cases must require `GLOW_SPARK_TESTS=1`; CI must set it explicitly.

- [ ] **Step 3: Implement BGEN probability loading.** Set `emitHardCalls=false`, preserve sample IDs, and disable index use because conversion does not create BGI files. Validate unphased diploid biallelic probabilities before calculating `P1 + 2*P2`. Use Spark array operations for production calculations. Check positions against PVAR using `start + 1`. Calculate MAC from observed values, then impute from the retained samples. Report and remove all-missing and zero-variance rows.

- [ ] **Step 4: Run scalar tests now and Spark cases in Task 8 CI.** Apply MAC 5 only to association data. Keep all usable variants from the user-selected model panel.

- [ ] **Step 5: Commit only Task 3 files.** Commit message: `Validate Glow dosage arrays and variant metadata`.

## Task 4: Fit models and route chromosome offsets correctly

**Files:** Create `scripts/glow/analysis.py`, `tests/test_glow_analysis.py`.

**Interfaces:** Define `FittedBatch(regression: object, reduction: object, sample_ids: list[str], traits: list[str], trait_type: str)`.
Provide `fit_batch(model_frame: DataFrame, phenotypes: pandas.DataFrame, covariates: pandas.DataFrame, options: ScientificOptions) -> FittedBatch`.
Provide `predict_chromosome(fit: FittedBatch, chromosome: str) -> pandas.DataFrame`.
Provide `test_chromosome(frame: DataFrame, phenotypes: pandas.DataFrame, covariates: pandas.DataFrame, offsets: pandas.DataFrame, chromosome: str, options: ScientificOptions) -> DataFrame`.
Provide `release_fit(fit: FittedBatch) -> None` to release persisted model data between batches.

- [ ] **Step 1: Write failing model and offset tests.** Test quantitative and binary branches, mismatched sample/trait order, missing model chromosomes, an incorrect chromosome offset, and binary folds without both classes. Define independent quantitative least-squares and binary score-test reference cases with fixed covariates and offsets. Require relative tolerance `1e-6` and absolute tolerance `1e-8` for those small reference cases.

- [ ] **Step 2: Run failing interface tests without Spark.** Use spies to prove that chromosome and sample mismatches fail before a Glow call. Run `python3 -m unittest discover -s tests -p test_glow_analysis.py -v`.

- [ ] **Step 3: Implement the Glow adapter.** Use `block_variants_and_samples`, `RidgeReduction(..., cov_df=..., label_type=...)`, and its `fit_transform()` result. Use `RidgeRegression.from_ridge_reduction` or `LogisticRidgeRegression.from_ridge_reduction`, followed by `fit()`. Quantitative prediction uses `transform_loco(chromosomes=[chromosome])`; binary prediction also sets `response='linear'`. Reindex the returned predictions explicitly to the prepared sample order. Call Glow GWAS with matching `contigs`, `offset_df`, `values_column='values'`, and `intersect_samples=False`. Binary tests use `correction='approx-firth'` and the selected threshold. Preserve the statistics Glow actually returns; an uncorrected binary score test can lack an effect estimate.

- [ ] **Step 4: Run routing tests now and numerical cases in Task 8 CI.** Also force a corrected binary test in the smoke fixture and require a finite result. Do not compare a binary score statistic with a full-model Wald statistic.

- [ ] **Step 5: Commit only Task 4 files.** Commit message: `Add Glow model fitting and chromosome-specific tests`.

## Task 5: Stage inputs and verify complete remote outputs

**Files:** Create `scripts/glow/storage.py`, `tests/test_glow_storage.py`.

**Interfaces:** In `contracts.py`, add `InputArtifact(name: str, path: Path)`, `StagedArtifact(name: str, uri: str, sha256: str, size_bytes: int, generation: int)`, and `OutputArtifact(kind: str, uri: str, relative_path: str, sha256: str, size_bytes: int)`.
Provide `input_identity(inputs: list[InputArtifact], settings: dict[str, object]) -> dict[str, object]`.
Provide `stage_inputs(client: object, run_uri: str, identity: dict[str, object], inputs: list[InputArtifact]) -> list[StagedArtifact]`.
Provide `publish_outputs(client: object | None, run_uri: str, identity_sha256: str, outputs: list[OutputArtifact], software: dict[str, str]) -> None`.
Provide `collect_outputs(client: object | None, run_uri: str, identity_sha256: str, out_dir: Path) -> list[Path]`.
Inject the Cloud Storage client so tests can use an in-memory fake.
For CI completion tests, the last two functions also support `file://` prefixes through ordinary filesystem operations with `client=None`.
The production controller rejects `file://` staging prefixes.

- [ ] **Step 1: Write failing storage tests.** Cover concurrent identity creation, repeat staging, relocated local files, duplicate basenames, altered scientific settings, corrupt objects, missing output files, duplicate output names, and an absent completion marker. Require rejection of `../escape`, absolute paths, another bucket, and a neighboring run prefix.

- [ ] **Step 2: Run failing tests.** Run `python3 -m unittest discover -s tests -p test_glow_storage.py -v`.

- [ ] **Step 3: Implement immutable staging and completion.** Identify files by generated roles such as `model/bgen` and `association/0/bgen`, not local basenames. Hash canonical settings and named file content; exclude localized path strings and timestamps from identity. Use `if_generation_match=0` for identity creation and new immutable objects. Reuse existing objects only after checking content. Write `outputs.json` with `schema_version=1`, identity hash, software, and output records. Write `COMPLETE` last, containing the manifest SHA-256. During collection, verify identity, prefix, destination, size, and checksum before an atomic local rename. Collect association files into one local directory with globally unique names such as `src0_chr1_batch0_part00000.tsv.gz`; collect predictions into a separate directory with unique chromosome/batch names. Preserve remote partition metadata in the manifest. Do not delete staged data.

- [ ] **Step 4: Run storage tests.** A successful collection must return only verified local paths. Test a partial download followed by a retry and ensure no partial file enters a WDL output directory.

- [ ] **Step 5: Commit only Task 5 files.** Commit message: `Add immutable Spark staging and verified result collection`.

## Task 6: Submit and monitor one Dataproc batch

**Files:** Create `scripts/glow/controller.py`, `tests/test_glow_controller.py`.

**Interfaces:** Define `CloudOptions` with string fields `project`, `region`, `subnetwork_uri`, `service_account`, `staging_bucket`, `staging_prefix`, `run_id`, `terra_image`, and `spark_image`. Its integer fields are `driver_cores`, `driver_memory_mb`, `driver_overhead_mb`, `executor_cores`, `executor_memory_mb`, `executor_overhead_mb`, `min_executors`, `max_executors`, `driver_disk_gb`, `executor_disk_gb`, `ttl_seconds`, and `poll_seconds`. Defaults are specified in Step 3 and the Global Constraints.
Provide `make_batch(options: CloudOptions, staged: list[StagedArtifact], science: ScientificOptions, identity_sha256: str) -> google.cloud.dataproc_v1.Batch`.
Provide `run_controller(options: CloudOptions, science: ScientificOptions, inputs: list[InputArtifact], out_dir: Path, batch_client: object | None = None, storage_client: object | None = None) -> list[Path]`.

- [ ] **Step 1: Write failing controller tests.** Assert one create call after a timeout response followed by a successful lookup. Assert a retry attaches to the same job. Cover `AlreadyExists`, transient lookups, mismatched existing jobs, terminal failure, prior success, invalid output completion, deadline expiry, cancellation failure, and a catchable termination signal. Assert normal polling never sleeps less than 30 seconds; inject a clock and sleeper in tests.

- [ ] **Step 2: Run failing tests.** Run `python3 -m unittest discover -s tests -p test_glow_controller.py -v`.

- [ ] **Step 3: Implement CLI, request construction, and state handling.** Accept local file arrays through task-local newline lists and individual scalar CLI flags. Stage all input files and the image's `spark_main.py`. Use the staged `.py` URI as `main_python_file_uri`; project imports resolve from the custom image. Build typed API requests, not a CLI configuration file. Set runtime 1.2 and `execution_config.ttl=21600` seconds by default. Use `run_id + '-glow'` as the stable batch ID. Derive a UUID5 request ID from project, region, run ID, and input identity. Verify an existing batch's full request and identity before attachment. Require both image references to use the same content digest. Retry uncertain submissions with the same IDs and a bounded API backoff. Handle `PENDING`, `RUNNING`, `SUCCEEDED`, `FAILED`, `CANCELLING`, and `CANCELLED` explicitly. On timeout or catchable termination, cancel the long-running operation referenced by the batch; always report the batch resource and cancellation outcome. Do not substitute deletion of batch metadata for compute cancellation.

  Starting Spark settings: driver/executor 4 cores each, JVM memory 8 GiB each, overhead 4 GiB each,
  executor minimum 2 and maximum 4, and driver/executor disk 250 GiB each.
  Validate cores in `{4, 8, 16}`, total memory per core in `[1024, 7424]` MiB for the standard tier,
  executor bounds in `[2, 2000]`, disks at least 250 GiB, and lifetime in `[600, 1209600]` seconds.
  Set dynamic allocation bounds at submission time. Keep these as adjustable starting values.
  [Dataproc resource settings](https://docs.cloud.google.com/managed-spark/docs/concepts/spark-properties-serverless).

- [ ] **Step 4: Run controller tests.** Confirm no test uses live credentials or submits an API request. Verify that a failed job never produces a successful execution report.

- [ ] **Step 5: Commit only Task 6 files.** Commit message: `Add retry-safe Dataproc batch controller`.

## Task 7: Add the Spark entry point and Terra WDL

**Files:** Create `scripts/glow/spark_main.py`, `workflows/glow_regenie.wdl`, `tests/test_glow_wdl.py`, `tests/test_glow_cromwell.py`; extend `tests/test_terra_check.py`.

**Interfaces:** `spark_main.main(argv: list[str] | None = None) -> int` accepts `--model-bgen`, `--model-sample`, `--model-variants`, `--sample-map`, `--phenotypes`, `--covariates`, `--phenotype-batches`, `--output-prefix`, and `--identity-sha256`. Association inputs use aligned repeated `--association-bgen`, `--association-sample`, and `--association-variants` options. Scientific scalar flags use the same names as the preparation/controller CLIs. Define workflow `GlowRegenie`, struct `GlowPgenFiles` with `File pgen`, `File pvar`, `File psam`, and struct `GlowBgenFiles` with `File bgen`, `File sample`, `File variants`, and `Int source_index`.
Use tasks `PrepareGlowInputs`, `ConvertGlowModelGenotypes`, `ConvertGlowAssociationGenotypes`, and `RunGlowAnalysis`.

- [ ] **Step 1: Write failing WDL and entry-point tests.** For every task, start with cloud File values and render them to a local directory containing spaces, a quote, and literal shell substitution text. Do not pre-localize File members inside structs in the test helper. Test absent and present optional arguments, aligned file lists, empty association arrays, repeated basenames, and newline rejection. Check task logs and explicit `glob()` output collection. Extend static regression cases to include the new workflow.

- [ ] **Step 2: Run failing tests.** Run `python3 -m unittest discover -s tests -p test_glow_wdl.py -v` and `miniwdl check --strict workflows/glow_regenie.wdl`.

- [ ] **Step 3: Implement the Spark entry point.** Register Glow and validate runtime versions. Download small staged sample/phenotype/covariate/batch tables to the driver's working directory with the storage client before pandas reads them. Spark reads BGEN and large variant metadata directly from their URIs. Persist usable genotype frames to run-specific Parquet paths. Process phenotype batches in sequence. For each chromosome, write LOCO predictions and association parts with unique source/chromosome/batch identifiers. Release predictions between chromosomes and model data between batches. Inventory final objects, add counts and software metadata, then call Task 5's completion functions. Support a local file output prefix only for CI; production controller input validation requires Cloud Storage.

- [ ] **Step 4: Implement WDL task commands and outputs.** Pass top-level File task inputs, including genotype members split out at the call boundary. Use safely quoted named CLI arguments, task-local file lists, readable-input checks, and timestamped start/end messages. Commands invoke `python3 -m scripts.glow.prepare`, `.convert`, or `.controller` from `/opt/regenie-spark`. Calculate conversion disk from source bytes plus explicit export/temp allowance. Calculate controller disk from staged input sizes plus an explicit result allowance. Expose association parts and predictions with flat directory globs and unique basenames from Task 5. Expose metadata, reports, and logs as Files. Keep runtime image references required and validate digest syntax in the controller.

- [ ] **Step 5: Verify integration.** Run strict miniwdl checks, `python3 tools/check_terra_wdl.py workflows/regenie.wdl workflows/glow_regenie.wdl`, and the WDL localization tests. Run the production graph with Cromwell 87 using local computation/API probes; preserve its declarations, task input types, and output collection. Keep this evidence separate from cloud validation.

- [ ] **Step 6: Commit only Task 7 files and the required checker test changes.** Commit message: `Add Terra WDL and Spark analysis entry point`.

## Task 8: Build the source artifacts and run Spark tests in GitHub Actions

**Files:** Create `envs/glow/Dockerfile`, `envs/glow/environment.yml`, `envs/glow/build-environment.yml`, `scripts/glow/smoke_test.py`, `.github/workflows/glow-image.yml`, `tests/test_glow_image.py`; update `.github/workflows/wdl-validation.yml` and add a separate exact dependency file for Glow unit tests if needed.

**Interfaces:** Docker stages are `glow-builder`, `production`, and `ci-test`. The production image contains `/opt/regenie-spark`, the built Glow wheel, and `/opt/glow/jars/glow-assembly.jar`. Provide `smoke_test.main(argv: list[str] | None = None) -> int` with `--work-dir` and `--spark-master` arguments. Publish to `ghcr.io/aou-multiomics-analysis/regenie-glow:sha-<commit>` after successful smoke tests.

- [ ] **Step 1: Write failing artifact and CI checks.** Assert all direct dependencies have exact pins; production has no Spark/PySpark/Java installation; image scripts match automatic build path filters; no image publishes after a failed smoke test; and Python/JAR versions and source commits match. Include a vendor-import test that proves `import glow` resolves to the upstream package.

- [ ] **Step 2: Run failing checks.** Run `python3 -m unittest discover -s tests -p test_glow_image.py -v`. Do not build Docker locally.

- [ ] **Step 3: Implement the source build and image targets.** Fetch exactly Glow commit `72637ad17ff267ca6bbc63cfc77c78db619f38b5` with Git and verify its SHA. The published PyPI/Maven 2.0.3 URLs returned 404 during planning, so do not depend on those URLs. Use upstream sbt 1.9.9, Scala 2.12.19, and Spark 3.5.1 build settings. Generate the Python functions with `python/generateFunctions`, build the assembly with `core/assembly`, and build the wheel from `python/`. Keep JDK, sbt, source checkout, and build caches in the builder stage. Write `/opt/glow/build.json` with the source commit, wheel and JAR checksums, and build versions. Verify the wheel's 2.0.3 metadata and the artifact checksums before copying them into the final image.

  Pin NumPy 1.26.4, pandas 1.5.3, SciPy 1.11.4, statsmodels 0.14.4, PyArrow 15.0.2,
  nptyping 2.5.0, opt_einsum 3.4.0, and typeguard 2.13.3 for the initial compatibility test.
  Pin google-cloud-dataproc 5.9.0, google-cloud-storage 2.16.0, and google-api-core 2.19.1.
  Add exact pins for build helpers and required system utilities when resolving the image environment.
  Run `pip check` and record the full resolved environment. Do not silently change the approved Glow version.
  Configure UID/GID 1099, writable HOME, procps/tini, Python executable paths, and driver/executor JAR classpaths.
  The CI target adds Spark/PySpark 3.5.1 and Java 17; these additions must not enter the published production image.

- [ ] **Step 4: Implement synthetic integration tests.** Generate 240 samples, 200 variants on two chromosomes, two quantitative traits, and two binary traits with both encodings. Include non-integer dosages, missing genotype values, sample reorderings, and a nonconstant numerical covariate. Test real PLINK export, dosage differences at most `1e-4`, both model types, LOCO routing, MAC filtering, and complete outputs. Execute fixed-offset reference tests from Task 4. Compare canonical sorted results between one-worker and two-worker Spark clusters with `rtol=1e-6`, `atol=1e-8`; exclude timestamps and filenames. Require actual separate worker processes for the multiple-worker test. Allocate at least 1 GiB per test executor, cap BLAS threads, and run sessions sequentially.

- [ ] **Step 5: Implement the GitHub Actions workflows.** Upload versions, logs, counts, and failed outputs as CI artifacts.
  Automatic image build paths are `scripts/glow/**/*.py` only; manual dispatch handles environment or Dockerfile edits.
  Run CI on pull requests and the repository's existing main/`codex/**` branch pattern.
  Build the production and CI targets on the GitHub runner, test the CI target, and publish the same production content only after success.
  Extend WDL validation triggers to include nested Glow scripts and its exact unit-test dependencies.
  Extend the Cromwell job to execute both the existing and Glow localization probes.

- [ ] **Step 6: Verify CI definitions.** Run unit/static checks locally. Inspect actual GitHub Actions results after the implementation branch is pushed. If a required CI job cannot run, report that gap explicitly.

- [ ] **Step 7: Commit only Task 8 files and required CI changes.** Commit message: `Build Glow artifacts and test Spark workflow in GitHub Actions`.

## Task 9: Document setup and verify the delivered workflow

**Files:** Create `docs/glow-terra.md`, `examples/glow.quantitative.inputs.json`, `examples/glow.binary.inputs.json`, `tests/test_glow_examples.py`; update `README.md`, `.dockstore.yml`, and existing example-discovery tests only where required for multiple workflow names.

**Interfaces:** Dockstore adds workflow name `glow-regenie` with `/workflows/glow_regenie.wdl` and its two separate test parameter files. Existing REGENIE examples keep their existing workflow prefix.

- [ ] **Step 1: Write failing example/schema tests.** Check both new input examples against `GlowRegenie` inputs and types, required cloud settings, image digests, one trait type per run, and binary encoding. Assert the existing examples still validate against `Regenie`.

- [ ] **Step 2: Run failing tests.** Run `python3 -m unittest discover -s tests -p test_glow_examples.py -v`.

- [ ] **Step 3: Add examples and operator instructions.** Use clearly marked sample values for bucket, project, region, subnet, account, unique run ID, and image digest. Change `tests/test_examples.py` to read each workflow's `testParameterFiles` from its Dockstore entry instead of applying the `Regenie` schema to every JSON file. Check descriptor paths for all entries. Document the separate Terra task and Spark workload identities, permissions, APIs, network setup, registry copy and digest verification, call-caching disablement, resource inputs, complete-case behavior, autosome limits, controller/staging costs, and storage retention. Include remote job lookup and cancellation commands. Explain that complete-case selection and Glow calculations can differ from native REGENIE. State `Not tested on Terra` until a complete approved Terra run succeeds.

- [ ] **Step 4: Run final checks.** Run strict miniwdl checks for both workflows, the Terra static checker for both workflows, all Python tests, both Cromwell localization probes, and required GitHub Actions image/Spark tests. Run `git diff --check`. Require fresh results before claiming success; classify skips and cloud limitations explicitly.

- [ ] **Step 5: Commit and provide the review result.** Commit message: `Document Glow workflow setup and Terra inputs`. Use the required finishing/review skills at execution time. Provide the implementation branch, pull request if created, tested image digest, checks that passed, and the explicit Terra validation status. Attach every created pull request with the native artifact tool. Do not merge or submit cloud jobs.

## Plan Review and Execution

Tasks 1–4 define the scientific inputs and calculations. Tasks 5–6 define remote file and job control.
Task 7 connects those interfaces, Task 8 verifies the built environment, and Task 9 completes delivery.
The controller and scientific code share interfaces, so implement the plan in order.

Recommended execution method: native implementation in this chat, followed by an independent review when execution reaches that stage.
The written plan and the user's chosen execution method require review before implementation starts.

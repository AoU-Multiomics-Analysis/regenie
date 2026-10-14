# REGENIE for Terra

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/regenie/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/regenie/actions/workflows/docker-image.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/regenie/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/regenie/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/regenie/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/regenie/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/regenie/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/regenie/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->


Run REGENIE 4.1.2 with PLINK 2 PGEN, PVAR, and PSAM inputs in Terra.
The workflow uses WDL 1.0. It supports quantitative traits and binary traits.
Select one trait type for each run.

## Workflow

1. `FitModel` fits the whole-genome model with selected common variants.
2. `TestVariants` tests each association file set with the step 1 predictions.
3. Terra runs the association tasks in parallel. Each file set can contain one chromosome or the whole genome.

LOCO means leave one chromosome out. REGENIE omits the tested chromosome from each model prediction.
The pipeline performs single-variant association tests. Prepare genotype QC and the step 1 variant panel before use.

## Required inputs

| Input | Format | Purpose |
| --- | --- | --- |
| `step1_genotypes` | Struct with `pgen`, `pvar`, and `psam` File fields | Whole-genome model inputs |
| `step1_variant_ids` | Text file, one variant ID per line, no header | Selected common variants for step 1 |
| `association_genotypes` | Nonempty array of genotype structs | Variant file sets for step 2 |
| `phenotypes` | Text table with a header that starts with `FID IID` | Trait values |
| `phenotype_columns` | Comma-separated column names | Traits to analyze |
| `docker_image` | Image reference, preferably with a digest | Tested software environment |

Use uncompressed PVAR files. Each PGEN, PVAR, and PSAM triplet must describe the same genotype set.
The three input filenames do not need to share a prefix. The task creates a local prefix after localization.
Use valid, sorted, biallelic variants. Use the same reference build and allele definitions across the two steps.
The step 1 file set must cover the chromosomes tested in step 2.
Step 2 samples must have matching step 1 predictions and phenotype IDs.

Choose a suitable common-variant panel for step 1. Do not use the full rare-variant set for model fitting.
The supplied variant list selects IDs from `step1_genotypes`; it does not perform QC.
The REGENIE limit for step 1 is one million variants. This workflow does not override that limit.

Example phenotype file:

```text
FID IID bmi height
1001 1001 24.8 172.0
1002 1002 NA 165.2
```

Example binary phenotype file:

```text
FID IID disease
1001 1001 0
1002 1002 1
1003 1003 NA
```

Use IDs that match the PSAM identifiers as read by REGENIE. Do not change FID values between files.
If the PSAM file starts with `#IID` and omits FID, each task adds `FID=0` to its local PSAM copy.
Use `FID=0` with the same IID in phenotype, covariate, and sample-filter tables for these inputs.
Existing PSAM FID values are preserved. The source PSAM file is not changed.
If SEX is absent, the local copy adds `SEX=0` (unknown), which also permits minimal ID-only PSAM files with REGENIE 4.1.2.
Use `NA` for missing traits. The default binary encoding is `0=control, 1=case`.
Set `binary_encoding` to `12` for `1=control, 2=case`.
All selected phenotype columns must use the same trait type and binary encoding.
Column names can contain letters, digits, underscores, dots, and hyphens. Do not put spaces around commas.

## Optional inputs

- `covariates`: table with a header that starts with `FID IID`.
- `covariate_columns`: selected column names. If absent, REGENIE uses all covariate columns.
- `categorical_covariates`: names of selected categorical columns.
- `keep_samples` and `remove_samples`: text files with `FID IID` pairs and no header. Both steps use these filters. If both are supplied, the task removes the listed pairs from the keep list and passes one keep file to REGENIE.
- `association_variant_ids`: variant IDs to retain in step 2.
- `inverse_normal_transform`: apply rank inverse normal transformation to quantitative traits in both steps. Default: `false`.
- `approximate_firth`: use approximate Firth correction for binary step 2. Default: `true`.
- `firth_threshold`: apply Firth correction below this p-value. Default: `0.01`.
- `min_mac`: minimum minor allele count in step 2. Default: `5`.
- `min_case_count`: minimum binary case count. Default: `10`.

If a selected phenotype does not produce step 1 predictions or a step 2 association file, the task fails with a clear error.
Do not combine quantitative traits and binary traits in one run.

## Run in Terra

1. Use the tested image digest in the examples, or copy the digest from a later successful **Docker Image CI** run.
2. Make the image readable by Terra. A new GHCR package can be private. Set its visibility to public, or copy it to an accessible Google Artifact Registry.
3. Import `workflows/regenie.wdl` into Terra through Dockstore, or add it to the Terra method repository.
4. Start with [quantitative inputs](examples/quantitative.inputs.json) or [binary inputs](examples/binary.inputs.json).
5. Replace the example `gs://YOUR_BUCKET` paths, phenotype names, covariate names, and image digest.
6. Use the input JSON in the Terra method configuration. Review the compute settings before submitting a run.

The example JSON files are workflow inputs. Tasks use typed WDL inputs and named CLI arguments.
The examples contain placeholder paths. They do not refer to a public test dataset.

| Resource input | Step 1 default | Step 2 default |
| --- | --- | --- |
| CPU | 8 | 4 |
| Memory | 32 GB | 16 GB |
| Extra disk | 100 GB | 50 GB |
| Block size | 1000 | 400 |

Each task adds the genotype and main input file sizes to its extra disk setting.
Step 2 also includes the prediction file sizes. Increase extra disk for large temporary files and outputs.
`lowmem=true` reduces step 1 memory use by writing temporary predictions to disk.
Adjust memory and disk for the sample count, variant count, and number of traits.
The defaults are starting values, not validated All of Us capacity estimates.
The default settings use no preemptible attempts and allow one task retry.

## Outputs

| Output | Contents |
| --- | --- |
| `associations_by_shard` | Association File arrays in the same order as `association_genotypes` |
| `associations` | Flat array of gzip-compressed `.regenie.gz` tables |
| `prediction_list` | Step 1 trait-to-prediction list |
| `prediction_files` | Step 1 LOCO files |
| `step1_log`, `step2_logs` | REGENIE logs |
| `step1_task_log`, `step2_task_logs` | Task command logs |

Results remain separate by shard and phenotype. Use `associations_by_shard` to retain the input-to-output mapping.
Files from different shards can have the same basename. Keep the Terra task directory paths when downloading results.
The prediction list contains task-local path metadata. To reuse predictions, supply the LOCO files as explicit File inputs and rebuild the list.
Files referenced inside a text or JSON file are not automatically localized by Terra.

## Container and tests

The [image build and smoke test](https://github.com/AoU-Multiomics-Analysis/regenie/actions/runs/38519424381) passed, including PSAM inputs with no FID or SEX columns.
This image is readable without registry credentials:

```text
ghcr.io/aou-multiomics-analysis/regenie@sha256:4232a95edf04477deab443032cf8cff01a8fc3cd200e6a9505721230a2795c43
```

The container uses `mambaorg/micromamba:2.3.3` and these direct dependency versions:

- REGENIE `4.1.2` from bioconda.
- PLINK 2 `2.0.0a.6.9` from bioconda, used to create synthetic test data.
- Python `3.11.13` from conda-forge.

Direct versions are pinned. Transitive package builds are resolved during the image build.
Use the published image digest to retain the complete tested environment.
GitHub Actions builds the Linux AMD64 image and runs both REGENIE steps before it publishes the image.
The smoke test uses 240 synthetic samples, 200 variants, two chromosomes, two quantitative traits, and both binary encodings.
It also tests covariates, sample filters, inverse normal transformation, Firth correction, and prediction path replacement.
PR builds test the image without publishing it. Push and manual builds publish a `sha-<commit>` tag after the tests pass.

Automatic image builds run only when `scripts/*.py` changes. These scripts are copied into the image.
Use **Run workflow** after a Dockerfile or dependency-only change.
WDL, test, and documentation changes do not cause automatic image builds.
No local Docker build is required.

Run the checks locally with Python 3.11:

```bash
python -m pip install -r requirements-dev.txt
miniwdl check --strict workflows/regenie.wdl
python tools/check_terra_wdl.py workflows/regenie.wdl
python -m unittest discover -s tests -v
```

The scope checker rejects workflow-scope file-writing functions, including functions in calls, scatters, conditionals, and outputs.
Command tests simulate cloud-to-local File mapping in both tasks. They test safe quoting with special characters in localized paths.
The image smoke test runs real REGENIE. It uses different task directories and cloud metadata in the prediction list.

**Terra validation status:** The complete workflow has not been run on Terra. Syntax checks and synthetic tests do not prove Terra execution.
No Terra job has been submitted as part of this implementation.

## References

- [REGENIE input formats and options](https://rgcgithub.github.io/regenie/options/)
- [REGENIE two-step method](https://rgcgithub.github.io/regenie/overview/)
- [PLINK 2 file formats](https://www.cog-genomics.org/plink/2.0/formats)
- [Bioconda REGENIE package](https://bioconda.github.io/recipes/regenie/README.html)

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

1. `PrepareBatches` validates phenotype names and divides them into batches.
2. `FitModel` fits all selected traits in one job by default. Set `batch_step1=true` to fit each phenotype batch in a separate job. Every fit uses the whole-genome common-variant panel.
3. `SplitGenotypes` uses PLINK 2 to write one PGEN/PVAR/PSAM set per chromosome from each association input set. It applies `association_variant_ids` if supplied.
4. `CollectPredictions` collects the model outputs and groups the LOCO files by phenotype batch. `TestVariants` runs one job per chromosome and phenotype batch.

Model fitting and genotype splitting can run at the same time. Association testing waits for both.
You can supply the same whole-genome triplet to `step1_genotypes` and `association_genotypes`.
Each split task localizes its source triplet once and writes its chromosome files in sequence.
Each association task localizes its chromosome files and only its batch's LOCO files. The full phenotype table is retained; each job selects its batch by column name. Source tables are not modified.
The split files are also workflow outputs for reuse.

LOCO means leave one chromosome out. REGENIE omits the tested chromosome from each model prediction.
The pipeline performs single-variant association tests. Prepare genotype QC and the step 1 variant panel before use.

## Required inputs

| Input | Format | Purpose |
| --- | --- | --- |
| `step1_genotypes` | Struct with `pgen`, `pvar`, and `psam` File fields | Whole-genome model inputs |
| `step1_variant_ids` | Text file, one variant ID per line, no header | Selected common variants for step 1 |
| `association_genotypes` | Nonempty array of genotype structs | File sets to split by chromosome for step 2 |
| `phenotypes` | Text table with a header that starts with `FID IID` | Trait values |

Use uncompressed PVAR files. Each PGEN, PVAR, and PSAM triplet must describe the same genotype set.
The three input filenames do not need to share a prefix. The task creates a local prefix after localization.
The workflow keeps the `PgenFiles` input format and passes each member to a separate task `File` input.
This lets Cromwell map each genotype path before the command applies shell quoting.
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

- `phenotype_columns`: comma-separated trait names. Leave this field empty in Terra to use all columns after `FID IID` in the phenotype file. In an input JSON, omit the key or use `null`, `""`, or a whitespace-only string. Supply names such as `bmi,height` to run a subset. Both steps read the localized file header and use the same selection rules. All selected columns must use the chosen trait type and binary encoding.
- `phenotype_batch_size`: traits per association job. Default: `50`. Set to `0` to use one batch with all selected traits. Explicit selections retain their given order; automatic selection retains header order. Negative sizes fail before model fitting.
- `batch_step1`: default `false`, which shares one Step 1 fit across all batches. Set to `true` to run a Step 1 job for each batch. Separate jobs use less memory per job and can run in parallel, but repeat genotype localization and calculations. With missing phenotypes, separate fits can use different sample sets or cross-validation folds; do not assume exact agreement with a joint fit.
- `docker_image`: image string. Default: `ghcr.io/aou-multiomics-analysis/regenie:latest`. Omit the input to use that default. Supply another tag or a full `@sha256:...` reference to override it. `latest` can change; a digest fixes the software version for reproducible runs.
- `covariates`: table with a header that starts with `FID IID`.
- `covariate_columns`: selected column names. If absent, REGENIE uses all covariate columns.
- `categorical_covariates`: names of selected categorical columns.
- `keep_samples` and `remove_samples`: text files with `FID IID` pairs and no header. Both steps use these filters. If both are supplied, the task removes the listed pairs from the keep list and passes one keep file to REGENIE.
- `association_variant_ids`: variant IDs to retain during chromosome splitting and step 2. Omit this input to test all supplied variants that pass REGENIE filters. A source set with no retained variants causes a clear error.
- `inverse_normal_transform`: apply rank inverse normal transformation to quantitative traits in both steps. Default: `false`.
- `approximate_firth`: use approximate Firth correction for binary step 2. Default: `true`.
- `firth_threshold`: apply Firth correction below this p-value. Default: `0.01`.
- `min_mac`: minimum minor allele count in step 2. Default: `5`.
- `min_case_count`: minimum binary case count. Default: `10`.

If a selected phenotype does not produce step 1 predictions or a step 2 association file, the task fails with a clear error.
Do not combine quantitative traits and binary traits in one run.

For 500 traits, a batch size of 50, and 22 chromosomes, the workflow creates 220 association jobs.
It creates one Step 1 job by default, or ten when `batch_step1=true`.
Each association batch downloads and reads its chromosome genotypes again. More parallel jobs can reduce elapsed time while increasing total transfer and compute.
REGENIE already shares genotype calculations across phenotypes. Its native `--split-l0`, `--run-l0`, and `--run-l1` options parallelize Step 1 by genotype blocks; this workflow uses phenotype batches instead.
The batch size is a starting value and has not been benchmarked on All of Us data.

## Run in Terra

1. Omit `docker_image` to use `latest`, or supply the digest from a successful **Docker Image CI** run.
2. Make the image readable by Terra. A new GHCR package can be private. Set its visibility to public, or copy it to an accessible Google Artifact Registry.
3. Import `workflows/regenie.wdl` into Terra through Dockstore, or add it to the Terra method repository.
4. Start with [quantitative inputs](examples/quantitative.inputs.json) or [binary inputs](examples/binary.inputs.json).
5. Replace the example `gs://YOUR_BUCKET` paths, phenotype names, and covariate names. Set the batch size and Step 1 batching mode.
6. Use the input JSON in the Terra method configuration. Review the compute settings before submitting a run.

The example JSON files are workflow inputs. Tasks use typed WDL inputs and named CLI arguments.
The examples contain placeholder paths. They do not refer to a public test dataset.

| Resource input | Split default | Step 1 default | Step 2 default |
| --- | --- | --- | --- |
| CPU | `split_cpu=4` | 8 | 4 |
| Memory | `split_memory_gb=16` | 32 GB | 16 GB |
| Extra disk | `split_extra_disk_gb=20` | 100 GB | 50 GB |
| Block size | — | 1000 | 400 |

Each task adds the genotype and main input file sizes to its extra disk setting.
Step 2 also includes the prediction file sizes. Increase extra disk for large temporary files and outputs.
The split task reserves twice the genotype input size plus the optional variant list and extra disk.
This allows space for both the localized source and the chromosome files. These outputs also use bucket storage.
The splitter accepts human chromosomes 1–22, X, Y, XY, and MT, including `chr` prefixes and PLINK numeric aliases.
It detects chromosomes after the optional variant filter and does not make empty chromosome jobs.
Chromosome splitting does not perform QC or add support for other REGENIE chromosome models.
`lowmem=true` reduces step 1 memory use by writing temporary predictions to disk.
Adjust memory and disk for the sample count, variant count, and number of traits.
The defaults are starting values, not validated All of Us capacity estimates.
The default settings use no preemptible attempts and allow one task retry.
`PrepareBatches` uses one CPU, 1 GB RAM, and the phenotype input size plus 1 GB disk.
`CollectPredictions` uses one CPU, 2 GB RAM, and twice the LOCO input size plus the prediction lists and 5 GB disk. It writes copies with unique filenames so model jobs can safely reuse REGENIE's original basenames.

## Outputs

| Output | Contents |
| --- | --- |
| `associations_by_shard` | Association File arrays aligned with `association_jobs` |
| `associations` | Flat array of gzip-compressed `.regenie.gz` tables |
| `association_jobs` | Chromosome-set index, chromosome, phenotype-batch index, and trait names for each association job |
| `association_chromosomes` | Typed structs with chromosome label and PGEN/PVAR/PSAM File outputs |
| `phenotype_batches`, `model_batches` | Trait names and zero-based indices for the association and model batches |
| `phenotype_batch_manifest`, `model_batch_manifest` | JSON batch definitions |
| `prediction_batches`, `prediction_batch_manifest` | Batch-specific prediction lists and typed LOCO File outputs |
| `batch_preparation_log`, `prediction_collection_log` | Preparation and collection task logs |
| `split_manifests` | One JSON output manifest per original association input set |
| `split_plink_logs`, `split_task_logs` | PLINK logs grouped by source set, and split task logs |
| `prediction_list`, `prediction_files` | Combined trait-to-prediction list and LOCO files with unique basenames |
| `step1_logs`, `step1_task_logs` | Logs for every model job, aligned with `model_batches` |
| `step1_log`, `step1_task_log` | First model job's logs, retained for compatibility |
| `step2_logs`, `step2_task_logs` | Association logs aligned with `association_jobs` |

Results remain separate by chromosome and phenotype. `association_jobs`, `associations_by_shard`, and Step 2 logs use the same order: source input order, chromosome order (1–22, X, Y, XY, MT), then phenotype-batch order.
Each job's `chromosome_set_index` identifies its entry in `association_chromosomes`.
An input that already contains one chromosome produces one chromosome set. Separate inputs are not merged or deduplicated.
Files from different shards can have the same basename. Keep the Terra task directory paths when downloading results.
The prediction list contains task-local path metadata. To reuse predictions, supply the LOCO files as explicit File inputs and rebuild the list.
Files referenced inside a text or JSON file are not automatically localized by Terra.

## Container and tests

The default image is:

```text
ghcr.io/aou-multiomics-analysis/regenie:latest
```

The container uses `mambaorg/micromamba:2.3.3` and these direct dependency versions:

- REGENIE `4.1.2` from bioconda.
- PLINK 2 `2.0.0a.6.9` from bioconda, used to split association genotypes and create synthetic test data.
- Python `3.11.13` from conda-forge.

Direct versions are pinned. Transitive package builds are resolved during the image build.
Use the published image digest to retain the complete tested environment.
GitHub Actions builds the Linux AMD64 image and runs both REGENIE steps before it publishes the image.
The smoke test uses 240 synthetic samples, 200 variants, two chromosomes, two quantitative traits, and both binary encodings.
It tests explicit phenotype selection and selection of all columns with omitted, empty, and whitespace-only values.
It also tests covariates, sample filters, inverse normal transformation, Firth correction, and prediction path replacement.
The split test compares all synthetic genotype values, samples, variant IDs, and alleles with the source.
It checks optional variant selection, including a selection that retains only one chromosome.
It also compares joint and batched association results for complete quantitative and binary phenotypes, using both shared and separate Step 1 fits. The comparison checks variant identity, alleles, sample counts, effect estimates, standard errors, and test statistics.
PR builds test the image without publishing it. Push and manual builds publish a `sha-<commit>` tag after the tests pass.
Successful main-branch pushes also update `main` and `latest`. For a manual release, select **Promote the tested image to latest** when dispatching the image workflow. Feature-branch pushes and manual runs with promotion disabled do not change `latest`.

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
Command tests simulate cloud-to-local File mapping in all five tasks. They test safe quoting with special characters in localized paths.
GitHub Actions also runs the production workflow graph and command expressions with Cromwell 87 on the local backend.
This test checks that genotype files, optional files, and prediction files reach commands through localized task paths. It also checks that the split JSON output becomes typed chromosome File outputs and reaches the second scatter.
It runs both shared and separate model fits with two phenotype batches and two chromosomes. It checks that each association task receives only its batch's predictions and that identical model-output basenames do not mix traits.
It uses computation probes in place of PLINK and REGENIE and removes the Docker runtime for that test only.
The image smoke test runs real PLINK and REGENIE. It uses different task directories and cloud metadata in the prediction list.

**Terra validation status:** The complete workflow has not been run on Terra. Syntax checks and synthetic tests do not prove Terra execution.
No Terra job has been submitted as part of this implementation.

## References

- [REGENIE input formats and options](https://rgcgithub.github.io/regenie/options/)
- [REGENIE two-step method](https://rgcgithub.github.io/regenie/overview/)
- [PLINK 2 file formats](https://www.cog-genomics.org/plink/2.0/formats)
- [Bioconda REGENIE package](https://bioconda.github.io/recipes/regenie/README.html)

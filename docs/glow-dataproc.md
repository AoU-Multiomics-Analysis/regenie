# Glow REGENIE on Dataproc

This optional WDL 1.0 workflow runs Glow 2.0.3 on Dataproc Serverless.
For the default workflow that runs on one Terra VM, use [the Terra guide](glow-terra.md).
Its workflow name and input namespace are `GlowRegenieDataproc`.
Terra prepares the sample set and converts PGEN files to BGEN.
One controller task starts and monitors one Spark batch.
Spark fits models and runs association tests in phenotype batches.

**Validation status: Not tested on Terra.**
Local unit, file-localization, and WDL checks are available.
The GitHub image build, real PLINK/Spark tests, and Cromwell probes must pass
before the image can be used. A successful image run publishes its digest in
the GitHub Actions summary.
A complete Terra run is a separate acceptance test. No cloud jobs were submitted.

## Scope and scientific choices

Use diploid, biallelic variants on human autosomes 1–22.
The selected model panel must cover each chromosome in the association inputs.
Use the same reference build and allele definitions in all files.
Prepare genotype QC and a suitable common-variant model panel before this workflow.

This is a Glow implementation of the two-step whole-genome regression method.
It does not run the C++ REGENIE executable. Results can differ from native REGENIE.
The existing [native workflow](../workflows/regenie.wdl) remains available.
See the [Glow method documentation](https://glow.readthedocs.io/en/stable/tertiary/whole-genome-regression.html)
and [association tests](https://glow.readthedocs.io/en/stable/tertiary/regression-tests.html).

PLINK first writes a filtered PGEN set in the prepared sample order.
It then exports BGEN 1.2 with 16-bit probabilities and REF first.
Spark calculates ALT dosage from genotype probabilities.
It calculates observed MAC before it fills missing dosages with the variant mean.
All-missing and constant variants are removed. The association MAC threshold is 5.
The user-selected model panel has no additional MAC threshold.

Select one trait type per run. Binary values use 0=control and 1=case by default.
Set `binary_encoding=12` for 1=control and 2=case.
Each binary trait requires at least 10 cases and 10 controls.
Each model validation fold must also contain both classes.
Binary tests use approximate Firth correction below p=0.01 by default.
An uncorrected score test can have no effect estimate; keep the returned statistics.

The default missing-data policy is `error`.
The optional `complete_case` policy removes samples with a missing value in
**any** selected trait or covariate before phenotype batching.
This gives each batch the same sample set. It can remove more samples than a
trait-specific analysis. Review `preparation.json` before interpreting results.

The optional quantitative inverse-normal transform runs once before batching.
It uses average ranks and `norm.ppf((rank - 0.5) / n)`.
Categorical covariates use a sorted reference level and dummy columns.
Constant traits and rank-deficient covariate designs fail validation.
Prepared sample IDs preserve exact FID/IID pairs through a hash.
The workflow returns `sample_map.tsv` to map them back to the original strings.

## Files and inputs

Use [quantitative](../examples/glow.dataproc.quantitative.inputs.json) or
[binary](../examples/glow.dataproc.binary.inputs.json) as a template.
Replace the example paths, cloud resources, and column names. The examples pin
a tested Dataproc image digest at a placeholder Google registry path. Copy the
published image to your registry and use its digest reference.

| Input | Required contents |
| --- | --- |
| `step1_genotypes` | Struct with explicit PGEN, uncompressed PVAR, and PSAM File fields |
| `step1_variant_ids` | One model variant ID per line, without a header |
| `association_genotypes` | Nonempty array of PGEN/PVAR/PSAM structs |
| `phenotypes` | Whitespace-delimited table with FID, IID, and selected traits |
| `covariates` | Optional table with FID, IID, and selected covariates |
| `keep_samples`, `remove_samples` | Optional sample-ID tables |
| `association_variant_ids` | Optional variant-ID list applied to each association set |

Preserve leading zeros and use identical FID/IID strings in each table.
If PSAM omits FID, use FID=0 in the other tables.
Use `NA` for missing values. Column selections are comma-separated names.
Leave `phenotype_columns` absent to use all non-ID phenotype columns.
The same rule applies to covariate columns.

WDL passes each input file through a typed File input.
A task checks that its local input files are readable before computation.
Files listed inside a user JSON or text file are not automatically localized.
Do not replace the explicit genotype File fields with a list of cloud paths.

## Project and service account defaults

Leave `project` and `service_account` unset to use the Google project and
service account attached to the Terra controller task VM. The controller reads
these values from the VM metadata service and logs the resolved names.
You do not need a service-account key or a new account for this detection.
Explicit values remain available as overrides.
If metadata is unavailable, the controller stops and asks for the missing override.
It does not switch to the Compute Engine default service account.

Terra configures the identity for a Hail notebook or managed Spark environment.
This WDL starts a separate Dataproc Serverless batch. Its attached account must
have permission to submit that batch and act as the workload account. The
workload account must have access to its inputs, outputs, network, and image.
Automatic detection does not grant those permissions or create cloud resources.
See [Terra service accounts](https://support.terra.bio/hc/en-us/articles/7448594459931-How-to-use-a-service-account-in-Terra)
and [Google VM metadata](https://docs.cloud.google.com/compute/docs/metadata/predefined-metadata-keys).

## Google Cloud setup

Ask the workspace administrator to prepare these resources.
This guide does not grant roles or create resources.

1. Use the Terra task's Google project by default, or set `project` to an approved
   alternative. It must have sufficient regional CPU, memory, disk,
   and Dataproc batch quotas. Enable the Dataproc, Compute Engine, Cloud Storage,
   and Artifact Registry APIs.
2. Select a subnet in the batch region. Configure internal Spark traffic,
   Private Google Access, and access to Google APIs and the image registry.
   Use the [serverless network requirements](https://docs.cloud.google.com/managed-spark/docs/concepts/network).
3. Select a regional staging bucket and a protected run prefix.
   The prefix contains genotype, phenotype, prediction, and association data.
   Apply the same access rules as the source study.
4. Identify the **Terra task service account** that runs the controller.
   It must be able to submit, inspect, and cancel Dataproc batches; act as the
   selected workload service account; and read, create, and list objects in the
   run prefix. For runtime 1.2, `roles/dataproc.editor` and
   `roles/iam.serviceAccountUser` on the workload account are documented role
   choices. A custom role can restrict permissions.
5. The **Spark workload service account** defaults to the account attached to
   the Terra task VM. Set `service_account` only to use another approved account.
   For runtime 1.2 the workload account needs Dataproc Worker permissions,
   access to staged inputs,
   permission to write/list result objects, and Artifact Registry Reader on the
   image repository. The Dataproc service agent also needs its documented
   subnet and image access, including Shared VPC permissions where applicable.
   Check the [runtime-specific IAM rules](https://docs.cloud.google.com/managed-spark/docs/concepts/iam-serverless).
6. Confirm that Terra's task VMs can pull the chosen task image.
   Registry permissions for Terra and Dataproc are separate.

A standard Terra workspace does not by itself grant the controller permission
to start Dataproc batches. An administrator must supply that permission.

## Build and copy the image

The [Glow image workflow](../.github/workflows/glow-image.yml) builds the source
artifacts and Docker images on GitHub. It runs real PLINK export, dosage checks,
independent numerical reference tests, both trait types, and comparisons between
one-worker and two-worker Spark clusters. It publishes only after those checks pass.

Automatic builds run only when included files under `scripts/glow/` change.
For Dockerfile or environment changes, use the workflow's manual dispatch.
There is no local Docker build step.

Glow Python and JAR artifacts are built from commit
`72637ad17ff267ca6bbc63cfc77c78db619f38b5`.
The image records the source and artifact hashes in `/opt/glow/build.json`.
The production image uses micromamba with pinned conda-forge/bioconda packages.
Dataproc supplies Spark and Java. Use the published `sha-<commit>` image for this workflow. The `local-sha-<commit>` image is for the separate single-VM workflow.

After a successful build, copy the published production image to Google Artifact
Registry. Use the digest from that build, rather than the example digest:

```bash
source_image='ghcr.io/aou-multiomics-analysis/regenie-glow@sha256:639d9340dfcb9c20b515389b59c62d4301306cdc41d03b80858fff9014cb8462'
target_image='us-central1-docker.pkg.dev/PROJECT/REPOSITORY/regenie-glow:RELEASE'
crane copy "$source_image" "$target_image"
crane digest "$source_image"
crane digest "$target_image"
```

The two digests must match. [crane copy](https://github.com/google/go-containerregistry/blob/main/cmd/crane/doc/crane_copy.md)
retains the image digest.
Use a full `@sha256:...` reference for both `docker_image` and `spark_image`.
Both must identify the same content. The Spark reference must use a supported
Google registry. You can use the same Artifact Registry reference for both.

The workflow selects Dataproc runtime 1.2, which supplies Spark 3.5 and Scala 2.12.
Google can update the runtime subminor version. The application checks these
versions and records the resolved runtime. Test a runtime change before study use.

## Import and run in Terra

1. Import Dockstore workflow `glow-regenie-dataproc`, or import
   [workflows/glow_regenie_dataproc.wdl](../workflows/glow_regenie_dataproc.wdl) directly.
2. Load the matching example JSON and replace every placeholder.
3. Set the region, subnet, staging bucket, and prefix. Leave `project` and
   `service_account` unset to detect them from the Terra task VM, or supply overrides.
   Set a new `run_id` with 4–40 lowercase letters, digits, or hyphens.
   Start it with a letter.
4. Disable Terra call caching for initial use. Review input files and resources.
5. Submit the workflow only after the workspace owner approves the cloud run.

The controller uploads localized files to `staging_prefix/run_id`.
It uses stable batch ID `run_id-glow`.
A retry attaches to that batch only when the input hashes and request match.
A terminal failed or cancelled batch needs a new run ID.

## Resources and costs

These are adjustable starting values. They have not been benchmarked on a study.

| Setting | Default |
| --- | --- |
| Traits per batch | 50 |
| Variants per model block | 1,000 |
| Sample blocks | 10 |
| Driver and executor cores | 4 each |
| Driver and executor JVM memory | 8 GiB each |
| Driver and executor overhead | 4 GiB each |
| Executors | Minimum 2, maximum 4 |
| Driver and executor disk | 250 GiB each |
| Remote lifetime | 6 hours |
| Controller polling interval | 30 seconds |
| Terra conversion VM | 4 CPUs, 16 GiB memory |
| Conversion extra disk allowance | 100 GiB per source |
| Controller extra disk allowance | 100 GiB |

Increase conversion disk for the temporary sorted PGEN set, BGEN export,
and other temporary files.
Increase controller disk for downloaded results.
Phenotype batches run in sequence. Spark distributes calculations within each batch.
More phenotypes increase the number of model fits and output files.

The Terra controller VM remains active while Spark runs.
Charges include Terra conversion/controller VMs, Spark driver/executors,
staging objects, Parquet work files, result downloads, and storage retention.
No cost or speed improvement over native REGENIE is established yet.

## Outputs and recovery

WDL returns explicit association parts, LOCO predictions, logs, reports, and
the sample map. Association parts have unique source/chromosome/batch names.
Predictions use the prepared sample IDs.
Use the returned manifest to identify remote partitions and verify checksums.
Cloud paths in that manifest are metadata; WDL collects local files through
declared File outputs and directory globs.

The worker writes `outputs.json`, then writes `COMPLETE` last.
The controller verifies identity, size, paths, and SHA-256 before returning results.
Remote success without matching, complete output files fails the workflow.

Use the execution report to find the batch and its operation:

```bash
gcloud dataproc batches describe RUN_ID-glow --region REGION --project PROJECT
gcloud dataproc batches cancel RUN_ID-glow --region REGION --project PROJECT
```

[Cancel](https://docs.cloud.google.com/sdk/gcloud/reference/dataproc/batches/cancel)
stops the workload and retains its metadata.
The controller also requests cancellation on a catchable termination signal or
its monitoring deadline. Forced VM termination can prevent that request.
The remote lifetime remains the final compute limit.

Staged files and work files are retained. The workflow does not delete them.
Set a bucket lifecycle policy after choosing the required evidence-retention period.
Do not delete an active run prefix.

## Validation layers

- Local tests: sample identity, scientific options, controller states, checksums,
  safe command rendering, cloud-to-local File substitution, and example schemas.
- Static WDL checks: strict miniwdl checks and a regression check that rejects
  workflow-scope file-writing functions.
- GitHub Actions: Cromwell localization probes and real image/Spark tests.
  These must succeed on the implementation commit before use.
- Terra: a complete approved run is still required. Local and CI checks do not
  establish Terra validation.

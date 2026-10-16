# Glow Spark workflow for Terra

## Purpose and agreed scope

Add a separate WDL 1.0 workflow for a standard Terra workspace on Google Cloud.
Terra will run preparation tasks and a controller task.
The controller will submit a Glow analysis to Dataproc Serverless, wait for completion, and collect the outputs.
Spark will distribute the analysis across worker machines.

The user approved this execution method on October 15, 2026.
This specification defines the proposed implementation for review.
The existing C++ REGENIE workflow remains available as a separate analysis method.
Do not claim numerical equivalence between Glow and C++ REGENIE without a comparison that establishes it.

Support quantitative traits and binary traits, with one trait type per run.
Accept the repository's PGEN, PVAR, and PSAM input format.
Perform whole-genome model fitting, leave-one-chromosome-out predictions, and single-variant association tests.
LOCO means leave one chromosome out.
Each association test must use the prediction that excludes its chromosome.

The first release supports biallelic, diploid variants on human autosomes 1–22.
Reject other chromosomes and unsupported ploidy before model fitting.
Genotype QC, ancestry selection, gene tests, interaction tests, and sex-chromosome models are outside this release.
The user must supply a suitable common-variant panel for model fitting.

## Execution method and software

Use Glow 2.0.3 for distributed model fitting and association testing.
Use Dataproc runtime `1.2`, whose documented components include Spark 3.5.1 and Scala 2.12.
This matches Glow's documented Spark requirement and Scala artifact.
Build the Glow Python package and matching Scala 2.12 JAR from the official 2.0.3 source release.
Pin the source to commit `72637ad17ff267ca6bbc63cfc77c78db619f38b5`.
[Glow installation](https://glow.readthedocs.io/en/stable/getting-started.html),
[Dataproc runtime components](https://docs.cloud.google.com/managed-spark/docs/concepts/versions/spark-runtime-1.2).

Specify runtime `1.2` in every submission.
Record the resolved runtime and image digest in the run metadata.
Google does not support pinning a serverless runtime's subminor version.
Do not describe the runtime as fully fixed by an image digest.
[Runtime version policy](https://docs.cloud.google.com/managed-spark/docs/concepts/versions/serverless-versions).

Use one remote Spark batch per workflow run.
Within that batch, process phenotype batches in sequence.
Distribute genotype calculations across Spark workers within each phenotype batch.
This stages the genotype data once and bounds the number of remote jobs.
The first release does not create a persistent cluster or require a notebook session.

## WDL tasks and data flow

1. `PrepareGlowInputs` validates tables, selects traits and samples, and creates an ordered sample map.
   It writes cleaned phenotype and covariate tables, phenotype batches, and a preparation report.
2. `ConvertGlowModelGenotypes` reads the whole-genome model PGEN triplet and the common-variant list.
   It applies the sample selection and exports BGEN and sample files.
3. `ConvertGlowAssociationGenotypes` runs in a WDL scatter over association PGEN triplets.
   Each task applies the same sample selection and exports BGEN and sample files.
4. `RunGlowAnalysis` receives all converted files as explicit WDL `File` inputs.
   It stages them in Cloud Storage, submits the Spark batch, monitors its state, and downloads the results.
5. Workflow outputs expose the collected files, conversion logs, preparation report, and remote job metadata.

Preparation must finish before genotype conversion.
Model conversion and association conversion can then run at the same time.
The controller must wait for all converted inputs.

Use typed structs for genotype sets.
Keep file members as `File` until command rendering.
An association set retains its source index in filenames and result metadata.
Do not merge or silently deduplicate separate source sets.

## Input interface

Required scientific inputs are:

| Input | WDL type | Meaning |
| --- | --- | --- |
| `step1_genotypes` | PGEN struct with three `File` members | Whole-genome model data |
| `step1_variant_ids` | `File` | Common variants selected for model fitting |
| `association_genotypes` | Nonempty array of PGEN structs | Data for association tests |
| `phenotypes` | `File` | Table beginning with `FID IID` |
| `trait_type` | `String` | `quantitative` or `binary` |

Required execution inputs are the Google project, region, regional subnet, workload service account,
staging bucket, staging prefix, unique `run_id`, and image digest references.
Keep the Dataproc staging bucket name separate from the application staging `gs://` prefix.
Require a nonempty application prefix within the bucket.
Validate project, region, subnet, and identifier formats before submission.

Use a required run identifier of 4–40 lowercase letters, digits, and hyphens, starting with a letter.
Derive the remote batch ID from this identifier.
Require a new identifier for a new analysis.
A retry of the same analysis must retain the identifier.

Support optional phenotype selection, numerical covariates, categorical covariates, keep/remove sample files,
association variant selection, and binary encoding `01` or `12`.
An omitted phenotype selection means all columns after `FID IID`.
Use named CLI arguments for task scripts.
Add an optional CLI argument only when its WDL input is present.

Use these scientific defaults: 50 traits per phenotype batch, 1,000 model variants per block,
10 sample blocks, minimum minor allele count 5, and minimum binary case and control counts of 10 each.
Allow a phenotype batch size of zero to select all traits in one batch.
Apply approximate Firth correction below p-value 0.01 by default for binary association tests.
Expose the random seed and quantitative inverse normal transformation as explicit inputs.
These defaults are starting values; they are not capacity or performance estimates.

Expose preparation CPU, memory, and extra disk separately from Spark driver and executor resources.
Use explicit executor bounds and a remote job lifetime, with a six-hour lifetime as the starting default.
Validate the resource values against the selected runtime's supported settings.

## Sample identity and missing data

Treat `FID` and `IID` as strings and use the pair as the sample identity.
For PSAM inputs without FID, use `FID=0`, consistent with the existing workflow.
Reject duplicate pairs and sample IDs that cannot be represented safely in the input tables.

Create stable internal sample IDs and an explicit mapping back to `FID IID`.
Do not form IDs by joining the two values with an ambiguous separator.
Apply the same retained sample order to model and association inputs.
Fail if an association set lacks a retained model sample.
Extra genotype samples can be excluded through the prepared sample selection.

Glow documents complete phenotype and covariate tables as model inputs.
Use `missing_data_policy=error` by default.
Report missing values among the selected samples and columns before conversion.
Offer `complete_case` as an explicit alternative.
That policy retains samples complete across all selected traits and covariates in the run.
Report the number removed and the final case/control counts for every binary trait.
Do not impute phenotype or covariate values.
[Glow model inputs](https://glow.readthedocs.io/en/stable/tertiary/whole-genome-regression.html).

Apply keep/remove filters before the missing-data policy.
Create phenotype batches after the final sample set is fixed.
Use the same sample set for all batches.
Reject a constant phenotype, an invalid binary value, or insufficient cases or controls.
Use deterministic sample blocks and check that binary model training subsets contain both classes.

Encode categorical covariates with a recorded reference category and fixed column order.
Reject missing categories and a rank-deficient covariate design.
Apply an optional quantitative transformation once, after sample selection.
Use that prepared table for model fitting and association testing.

## Genotype conversion and allele handling

Export BGEN 1.2 with 16-bit probabilities and the reference allele first.
Export the associated sample file and verify its identifiers against the prepared sample map.
Preserve variant IDs, chromosome labels, positions, reference alleles, and alternate alleles.
Check for duplicate variant IDs within each source set.
PLINK supports these BGEN export settings.
[PLINK export options](https://www.cog-genomics.org/plink/2.0/data#export).

Read BGEN probabilities in Glow and calculate alternate-allele dosage from them.
Do not replace dosages with thresholded genotype calls.
For unphased diploid biallelic records, dosage is `P(REF/ALT) + 2 * P(ALT/ALT)`.
Check phasing and ploidy explicitly before applying this formula.
Reject unsupported records rather than apply a diploid formula to them.
Glow assumes that the first BGEN allele is the reference allele.
[Glow BGEN reader](https://glow.readthedocs.io/en/stable/etl/variant-data.html).

Use the same dosage definition in both analysis steps.
Mean-impute missing genotype dosages using the retained analysis samples for each variant.
Exclude variants with no observed dosage or zero variance and report their counts.
Calculate the minor allele count from observed dosages before genotype imputation.
Keep the effect allele explicit in every association output.

Conversion tests must include non-integer dosages, missing genotypes, reversed allele frequencies, and reordered samples.
Require an absolute dosage difference no greater than `1e-4` for nonmissing synthetic values after BGEN conversion.
Hardcall-only fixtures do not establish dosage preservation.

## Spark analysis and outputs

Persist prepared genotype data as Parquet within the run's remote working prefix.
Reuse this data across phenotype batches.
Use Glow's model reduction and ridge models for each selected trait type.
Generate LOCO predictions one chromosome at a time to limit driver memory use.
Pass the matching predictions to Glow association tests.
For binary traits, use the linear genetic prediction as the offset and refit covariate effects in the test.
[Glow predictions](https://glow.readthedocs.io/en/stable/tertiary/whole-genome-regression.html),
[Glow association tests](https://glow.readthedocs.io/en/stable/tertiary/regression-tests.html).

Write association results as partitioned, gzip-compressed TSV files.
Partition by source set, chromosome, and phenotype batch.
Include variant identity, effect allele, phenotype, sample count, effect estimate, standard error,
and p-value when the selected test provides them.
Record any additional Glow test statistics with their documented names.
Do not label a Glow output as a native `.regenie` file.

Write prediction files, model validation summaries, sample counts, software versions, and an output manifest.
The manifest describes generated objects and records their sizes and checksums.
Publish a completion marker only after all required outputs are present.
The controller must check the marker and output manifest after remote success.

Download final outputs into distinct local directories in the controller task.
Require manifest objects to remain within the run prefix and reject unsafe local output paths.
Expose association and prediction files with explicit WDL `File` declarations or `glob()` outputs.
Expose the manifest and execution report as separate `File` outputs.
Keep remote object locations as string metadata.
Do not derive WDL `File` arrays from paths embedded in a manifest.

## Terra localization and submission boundary

Create files only in task scope or task commands.
Do not use workflow-scope `write_lines`, `write_tsv`, `write_map`, or `write_json` calls.
This restriction also applies inside workflow scatters, conditions, call expressions, and outputs.

Every local file opened by a preparation, conversion, or controller script must be an explicit WDL `File` input.
Check readability before computation and identify an unresolved cloud URI as a localization error.
Render named CLI arguments with safe shell quoting.
For arrays of files, create task-local newline-delimited lists after localization.
Reject newline characters in input paths used by these lists.

The controller uploads localized files before it submits the remote job.
Spark workers use the resulting Cloud Storage objects as remote inputs.
Those URI strings describe explicitly staged objects; they are not local filesystem paths.
Files referenced by remote arguments must be staged from the declared WDL inputs.
No downstream WDL task may treat those URI strings as automatically localized files.

Use the Dataproc Python client to build a typed batch request.
Pass scientific options to the Spark application as named CLI arguments.
JSON is allowed for the Dataproc API request and structured output metadata.
It must not become a generic argument wrapper between WDL tasks and analysis scripts.

## Remote job state and recovery

Store a run identity containing input checksums, scientific settings, image digests, and the application version.
Write that identity to the run prefix with an object creation precondition.
Fail if an existing identifier refers to different inputs or settings.
Stage inputs with fixed object names and verify checksums before reusing them.
Do not overwrite objects from another analysis.

A controller retry must attach to the same remote batch.
It must not create a second batch when the first submission has an uncertain response.
For an existing successful batch, verify the output identity and collect the results.
For a failed or cancelled batch, report its terminal state.
A new compute attempt after terminal failure requires a new run identifier.

Log batch creation, state changes, transfer counts, and output verification.
Log references to remote diagnostics without printing credentials or phenotype values.
Use bounded API retries and a polling interval of at least 30 seconds.
Fail on remote failure, an invalid output manifest, missing outputs, or checksum errors.

Set `environmentConfig.executionConfig.ttl` on the batch.
On controller timeout or a catchable termination signal, request cancellation and report its result.
Forced VM termination can prevent signal handling; do not promise immediate cancellation in that case.
The remote lifetime remains the independent limit.
Include the remote batch ID and a manual cancellation command in the run instructions.
[Dataproc execution configuration](https://docs.cloud.google.com/managed-spark/docs/reference/rest/v1/EnvironmentConfig).

Disable Terra call caching for the initial release's workflow submissions.
Controller retries use the run identity described above.
Do not automatically delete staging data or failed outputs.
Document their location and the separate storage lifecycle settings that an operator can apply.

## Container and access requirements

Build a separate image from `mambaorg/micromamba:2.3.3`.
Use Python 3.11.13 and PLINK 2 `2.0.0a.6.9`, consistent with the existing environment.
Pin all direct Python and conda dependencies in committed environment files.
Use conda-forge and bioconda where packages are available.
Build and install the Glow Python package at version 2.0.3 from the pinned source commit.
Build the matching JAR and package its dependencies during the GitHub Actions image build.
Use sbt 1.9.9 and Scala 2.12.19, as specified by that source release.
[Glow source release](https://github.com/projectglow/glow/releases/tag/v2.0.3).
The production job must not depend on downloading Maven packages when it starts.

The image contains the preparation scripts, controller, Spark application, and required Python libraries.
Configure the micromamba Python path for both Spark driver and executors.
Include the required system utilities and permissions for Dataproc's UID/GID 1099.
Do not install Spark, PySpark, or a Java runtime in the production image.
Dataproc supplies those components when it runs the batch.
Use a separate CI test target to supply Spark and Java for synthetic tests.

Publish the tested image by commit tag and digest.
Terra can use a readable image registry for its tasks.
Dataproc must use a copy in Google Artifact Registry or another Google registry supported by the service.
Verify that both image references identify the same tested image content.
[Dataproc custom containers](https://docs.cloud.google.com/managed-spark/docs/guides/custom-containers).

Document the permissions for the Terra task identity and the Dataproc workload identity separately.
The submitter needs batch submission, monitoring, cancellation, and permission to act as the workload service account.
The workload identity needs Dataproc Worker permissions plus access to its staging data and container image.
Both identities need the relevant Cloud Storage access for staging and result collection.
The project needs the required APIs, network access, and available compute quota.
The workflow will not create service accounts, grant IAM roles, or change network policies.
[Serverless permissions](https://docs.cloud.google.com/managed-spark/docs/concepts/iam-serverless).

Input localization, conversion, staging uploads, and result downloads require additional transfer and disk space.
Size conversion disks for the source files, exported BGEN files, and temporary output.
Size the controller disk for its localized inputs and collected final outputs.
Document the continued Terra controller VM charge while Spark runs.
Measure complete-workflow cost and elapsed time before recommending production resource settings.

## Verification and delivery

Add `workflows/glow_regenie.wdl`, separate scripts under `scripts/glow/`, a separate image definition,
Terra input examples, a Dockstore entry, and setup instructions.
Keep the new image name and its build workflow distinct from the C++ REGENIE image.

Automatic image builds must trigger only when a script included in that image changes.
Keep manual dispatch available for Dockerfile or dependency changes.
Run image builds and synthetic Spark smoke tests in GitHub Actions.
Do not build a local Docker image for development checks.

Required verification includes:

- Strict WDL 1.0 validation and the existing Terra static checker applied to the new workflow.
- Regression cases for workflow-scope file writing and manifest-derived `File` outputs.
- Cloud-to-local substitution tests for every new WDL task, including changed paths and unusual filenames.
- Dosage, sample-order, chromosome, and effect-allele checks across conversion and Spark loading.
- Quantitative and binary smoke tests that execute model fitting, LOCO prediction, and association testing.
- Independent reference tests for association calculations with fixed offsets and covariates.
- Tests that identify incorrect chromosome offsets and shuffled phenotype/sample order.
- Mocked cloud tests for submission uncertainty, retry attachment, timeout, cancellation, and incomplete outputs.
- Numerical comparison of one-worker and multiple-worker synthetic results within an explicit test tolerance.
- Preservation of all existing C++ workflow tests.

Record software versions, test settings, sample counts, and retained variant counts with each smoke test.
Publish the image only after the smoke tests pass.
Use local Python and WDL checks where available; use CI for the image and Spark execution tests.

Complete-workflow validation requires an approved Terra submission with the required cloud setup.
No cloud submission is authorized during implementation.
Report the workflow as not tested on Terra until a complete Terra run succeeds.
Local Spark tests, mocked cloud tests, and WDL syntax checks cannot establish that result.

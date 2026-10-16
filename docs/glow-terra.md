# Glow REGENIE on one Terra VM

The default [GlowRegenie workflow](../workflows/glow_regenie.wdl) uses WDL 1.0.
Terra prepares the sample set and converts PGEN files to BGEN.
One analysis task runs Spark with `local[N]` on its Terra VM.
`N` is `analysis_cpu`. This follows the single-VM execution style in
[MTtoVCF](https://github.com/AoU-Multiomics-Analysis/MTtoVCF/blob/main/scripts/ExportVCF.py).
Spark uses several CPU threads on this VM. It does not start other VMs.

Use Terra's normal workspace access and billing setup.
There are no project, service account, region, subnet, or staging-bucket inputs.
Cromwell localizes input files and uploads declared output files.
The optional [Dataproc workflow](glow-dataproc.md) remains available for a cluster.

**Validation status: Not tested on Terra.**
GitHub Actions checks real PLINK and Spark calculations and Cromwell file handling.
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

Use [quantitative](../examples/glow.quantitative.inputs.json) or
[binary](../examples/glow.binary.inputs.json) as a template.
Every example value is a placeholder. The repeated `a` image digest is invalid
for use and must be replaced with the digest from a successful image build.

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

## Import and run in Terra

1. Import Dockstore workflow `glow-regenie`, with
   [workflows/glow_regenie.wdl](../workflows/glow_regenie.wdl) as its descriptor.
   Include its `glow_regenie_dataproc.wdl` import; the default workflow calls only
   the shared preparation and conversion tasks from that file.
2. Load the quantitative or binary example. Replace the input paths and column names.
3. Set `docker_image` to the **Terra local image** digest from a successful
   [Glow Spark Image build](../.github/workflows/glow-image.yml).
   Its tag starts with `local-sha-`. This image includes Spark and Java.
   The `sha-` image is for Dataproc, which supplies its own Spark runtime.
4. Make the image readable by Terra. If GHCR access is restricted, copy the same
   image to an accessible Google Artifact Registry and use its digest reference.
5. Review the resource settings before an approved Terra submission.

## Resources

These are adjustable starting values. They have not been benchmarked on study data.
The analysis defaults follow the VM sizes used by MTtoVCF.

| Input | Default | Use |
| --- | --- | --- |
| `analysis_cpu` | 64 | VM CPUs and Spark local threads |
| `analysis_memory_gb` | 256 | Total VM memory in GB |
| `analysis_extra_disk_gb` | 1000 | Disk allowance beyond localized input sizes |
| `spark_driver_memory_gb` | 64 | Java heap on this VM |
| `spark_parallelism` | 100 | Default Spark partition count |
| `spark_shuffle_partitions` | 100 | SQL shuffle partition count |
| `prepare_cpu`, `prepare_memory_gb` | 1, 8 | Cohort preparation resources |
| `convert_cpu`, `convert_memory_gb` | 4, 16 | Resources per PGEN conversion task |
| `convert_extra_disk_gb` | 100 | Conversion disk allowance |
| `batch_size` | 50 | Traits per batch |
| `variants_per_block` | 1000 | Variants per model block |
| `sample_block_count` | 10 | Sample blocks for model fitting |

Keep driver heap below VM memory. Python workers, Arrow buffers, native libraries,
and the operating system also need memory. Local Spark shares one Java process;
executor counts do not create extra machines.
Increase analysis disk for Parquet work files, predictions, result parts, and the
checked output copies. Increase conversion disk for sorted PGEN and BGEN files.
Phenotype batches run in sequence. These defaults do not establish capacity for
All of Us data or a speed or cost benefit over native REGENIE.

## Outputs

WDL returns association parts, LOCO predictions, metrics, application logs,
`analysis_log`, preparation and conversion logs, the sample map, converted BGEN
files, an execution report, and an output manifest.
Association parts have distinct source, chromosome, and batch names.
Predictions use the prepared sample IDs. Use the sample map to recover FID/IID.

The runner checks the completion marker, run identity, file sizes, and SHA-256
before it collects outputs. The local output manifest records relative paths,
sizes, and checksums. Use WDL File outputs to reuse results in another task;
manifest paths are metadata. Temporary Spark and Parquet files stay on the task
disk and are removed after output collection.

## Image and tests

The image uses micromamba 2.3.3, pinned conda-forge and bioconda packages,
Spark 3.5.1, Java 17.0.15, and Glow 2.0.3 built from source commit
`72637ad17ff267ca6bbc63cfc77c78db619f38b5`.
The build records the source and artifact hashes in `/opt/glow/build.json`.
Use a full image digest to fix the complete environment.

GitHub Actions builds and tests both image variants before publication.
The tests use 240 synthetic samples, 200 variants, two chromosomes,
two quantitative traits, and both binary encodings.
They run the exact single-VM command used by WDL and compare its results with
one-worker and two-worker Spark runs. They also check dosages and numerical
reference results. These tests do not establish agreement with native REGENIE.
Cromwell 87 probes check file localization and output collection for both WDLs.
Strict WDL checks reject workflow-scope file-writing functions and
manifest-derived File outputs.

Automatic image builds run only when included `scripts/glow/**/*.py` files change.
Use manual dispatch after a Dockerfile or dependency-only change.
No local Docker build is required.

# REGENIE pipeline design

## Purpose

Build a WDL 1.0 workflow for Terra. Use PLINK 2 PGEN, PVAR, and PSAM files.
Support quantitative traits and binary traits. Push the work to the supplied GitHub repository.

## Data flow

1. Fit REGENIE step 1 with a whole-genome PGEN file set and a selected variant list.
2. Return the prediction list and every LOCO prediction file as typed File outputs.
3. Run step 2 for each association PGEN file set. Each set can contain one chromosome or the whole genome.
4. Return association tables, prediction files, REGENIE logs, and task logs.

Use separate file sets for the common variants used in step 1 and the variants tested in step 2.
Require a step 1 variant list. The user must select suitable common variants before the run.
This version performs single-variant tests. It does not perform genotype QC, gene tests, or ancestry selection.

## Inputs and models

Use typed WDL inputs. Preserve File types inside the genotype struct.
Use named command-line arguments. Use a task-local newline file list for prediction files.
Use one trait type per run. Multiple phenotype columns can use the same model.
Allow a covariate file, selected covariates, categorical covariates, and sample keep/remove files.
Use 0/1 binary encoding by default. Allow 1/2 encoding.
Use approximate Firth correction for binary step 2 by default.
Allow inverse normal transformation for quantitative traits, applied in both steps.
Expose CPU, memory, disk, block size, cross-validation folds, and retry settings.

## Terra file handling

Create files only within tasks. Do not call file-writing functions at workflow scope.
Check that each local input is readable before computation. Reject unresolved cloud URIs.
Create a fixed local PGEN prefix with links to localized files.
Read the step 1 prediction list as structured data. Rebuild its paths from explicitly localized prediction files.
Use safe shell quoting at command rendering. Do not put input paths in JSON argument wrappers.

## Container and CI

Use a pinned micromamba image. Install pinned REGENIE, PLINK 2, and Python versions from conda-forge and bioconda.
Build the image in GitHub Actions. Run synthetic quantitative and binary tests before publishing the image.
Automatic image builds start only when a script copied into the image changes.
Permit manual builds for changes to the Dockerfile or dependency file.
Use a SHA image tag and a digest for reproducible Terra runs.
Run WDL checks, Python tests, and a workflow-scope file-writing regression check in separate CI.
The tests must cover cloud-to-local path replacement, including both WDL tasks.
Do not submit a Terra job. Do not claim that this workflow has run on Terra.

## Delivery

Provide Terra input examples, a Dockstore descriptor, usage instructions, and test results.
Push a feature branch and open a reviewable pull request. Do not merge the pull request.

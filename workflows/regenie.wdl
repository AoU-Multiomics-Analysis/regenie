version 1.0

struct PgenFiles {
  File pgen
  File pvar
  File psam
}

struct ChromosomePgen {
  String chromosome
  File pgen
  File pvar
  File psam
}

struct PhenotypeBatch {
  Int batch_index
  String phenotype_columns
}

struct PredictionBatch {
  Int batch_index
  String phenotype_columns
  File prediction_list
  Array[File] prediction_files
}

struct PredictionBatchIndices {
  Int batch_index
  String phenotype_columns
  Int prediction_list_index
  Array[Int] prediction_file_indices
}

struct AssociationJob {
  Int chromosome_set_index
  String chromosome
  Int phenotype_batch_index
  String phenotype_columns
}

workflow Regenie {
  input {
    PgenFiles step1_genotypes
    File step1_variant_ids
    Array[PgenFiles]+ association_genotypes
    File phenotypes
    String? phenotype_columns
    String trait_type = "quantitative"
    String docker_image = "ghcr.io/aou-multiomics-analysis/regenie:latest"
    Int phenotype_batch_size = 50
    Boolean batch_step1 = false

    File? covariates
    String? covariate_columns
    String? categorical_covariates
    File? keep_samples
    File? remove_samples
    File? association_variant_ids

    String binary_encoding = "01"
    Boolean inverse_normal_transform = false
    Boolean approximate_firth = true
    Float firth_threshold = 0.01
    Float min_mac = 5.0
    Int min_case_count = 10
    Int cv_folds = 5
    Boolean lowmem = true
    Int step1_block_size = 1000
    Int step2_block_size = 400

    Int step1_cpu = 8
    Int step2_cpu = 4
    Int split_cpu = 4
    Int step1_memory_gb = 32
    Int step2_memory_gb = 16
    Int split_memory_gb = 16
    Int step1_extra_disk_gb = 100
    Int step2_extra_disk_gb = 50
    Int split_extra_disk_gb = 20
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  call PrepareBatches {
    input:
      phenotypes = phenotypes,
      phenotype_columns = phenotype_columns,
      batch_size = phenotype_batch_size,
      batch_step1 = batch_step1,
      docker_image = docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  scatter (model_batch in PrepareBatches.model_batches) {
    call FitModel {
      input:
        pgen = step1_genotypes.pgen,
        pvar = step1_genotypes.pvar,
        psam = step1_genotypes.psam,
        variant_ids = step1_variant_ids,
        phenotypes = phenotypes,
        phenotype_columns = model_batch.phenotype_columns,
        trait_type = trait_type,
        docker_image = docker_image,
        covariates = covariates,
        covariate_columns = covariate_columns,
        categorical_covariates = categorical_covariates,
        keep_samples = keep_samples,
        remove_samples = remove_samples,
        binary_encoding = binary_encoding,
        inverse_normal_transform = inverse_normal_transform,
        min_case_count = min_case_count,
        cv_folds = cv_folds,
        lowmem = lowmem,
        block_size = step1_block_size,
        cpu = step1_cpu,
        memory_gb = step1_memory_gb,
        extra_disk_gb = step1_extra_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  call CollectPredictions {
    input:
      prediction_lists = FitModel.prediction_list,
      prediction_files = flatten(FitModel.prediction_files),
      prediction_file_counts = FitModel.prediction_file_count,
      phenotype_columns = PrepareBatches.selected_phenotypes,
      batch_size = phenotype_batch_size,
      docker_image = docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  scatter (batch_indices in CollectPredictions.batch_indices) {
    scatter (prediction_index in batch_indices.prediction_file_indices) {
      File batch_prediction_file = CollectPredictions.collected_prediction_files[prediction_index]
    }
    PredictionBatch prediction_batch = object {
      batch_index: batch_indices.batch_index,
      phenotype_columns: batch_indices.phenotype_columns,
      prediction_list: CollectPredictions.batch_prediction_lists[batch_indices.prediction_list_index],
      prediction_files: batch_prediction_file
    }
  }
  Array[PredictionBatch] collected_prediction_batches = prediction_batch

  scatter (genotype_set in association_genotypes) {
    call SplitGenotypes {
      input:
        pgen = genotype_set.pgen,
        pvar = genotype_set.pvar,
        psam = genotype_set.psam,
        variant_ids = association_variant_ids,
        docker_image = docker_image,
        cpu = split_cpu,
        memory_gb = split_memory_gb,
        extra_disk_gb = split_extra_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  Array[File] chromosome_pgen_files = flatten(SplitGenotypes.pgen_files)
  Array[File] chromosome_pvar_files = flatten(SplitGenotypes.pvar_files)
  Array[File] chromosome_psam_files = flatten(SplitGenotypes.psam_files)
  scatter (chromosome_index in range(length(chromosome_pgen_files))) {
    ChromosomePgen chromosome_file_set = object {
      chromosome: sub(basename(chromosome_pgen_files[chromosome_index], ".pgen"), "^chr", ""),
      pgen: chromosome_pgen_files[chromosome_index],
      pvar: chromosome_pvar_files[chromosome_index],
      psam: chromosome_psam_files[chromosome_index]
    }
  }
  Array[ChromosomePgen] chromosome_genotypes = chromosome_file_set

  Array[Pair[Int, Int]] job_indices = cross(range(length(chromosome_genotypes)), range(length(collected_prediction_batches)))

  scatter (job_index in job_indices) {
    ChromosomePgen chromosome_set = chromosome_genotypes[job_index.left]
    PredictionBatch phenotype_batch = collected_prediction_batches[job_index.right]
    AssociationJob association_job = object {
      chromosome_set_index: job_index.left,
      chromosome: chromosome_set.chromosome,
      phenotype_batch_index: phenotype_batch.batch_index,
      phenotype_columns: phenotype_batch.phenotype_columns
    }
    call TestVariants {
      input:
        pgen = chromosome_set.pgen,
        pvar = chromosome_set.pvar,
        psam = chromosome_set.psam,
        phenotypes = phenotypes,
        phenotype_columns = phenotype_batch.phenotype_columns,
        trait_type = trait_type,
        docker_image = docker_image,
        prediction_list = phenotype_batch.prediction_list,
        prediction_files = phenotype_batch.prediction_files,
        covariates = covariates,
        covariate_columns = covariate_columns,
        categorical_covariates = categorical_covariates,
        keep_samples = keep_samples,
        remove_samples = remove_samples,
        variant_ids = association_variant_ids,
        binary_encoding = binary_encoding,
        inverse_normal_transform = inverse_normal_transform,
        approximate_firth = approximate_firth,
        firth_threshold = firth_threshold,
        min_mac = min_mac,
        min_case_count = min_case_count,
        block_size = step2_block_size,
        cpu = step2_cpu,
        memory_gb = step2_memory_gb,
        extra_disk_gb = step2_extra_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  output {
    File prediction_list = CollectPredictions.prediction_list
    Array[File] prediction_files = CollectPredictions.collected_prediction_files
    File step1_log = FitModel.regenie_log[0]
    File step1_task_log = FitModel.task_log[0]
    Array[File] step1_logs = FitModel.regenie_log
    Array[File] step1_task_logs = FitModel.task_log
    Array[PhenotypeBatch] model_batches = PrepareBatches.model_batches
    Array[PhenotypeBatch] phenotype_batches = PrepareBatches.phenotype_batches
    File phenotype_batch_manifest = PrepareBatches.phenotype_batch_manifest
    File model_batch_manifest = PrepareBatches.model_batch_manifest
    File batch_preparation_log = PrepareBatches.task_log
    Array[PredictionBatch] prediction_batches = prediction_batch
    File prediction_batch_manifest = CollectPredictions.manifest
    File prediction_collection_log = CollectPredictions.task_log
    Array[AssociationJob] association_jobs = association_job
    Array[Array[File]] associations_by_shard = TestVariants.associations
    Array[File] associations = flatten(TestVariants.associations)
    Array[File] step2_logs = TestVariants.regenie_log
    Array[File] step2_task_logs = TestVariants.task_log
    Array[ChromosomePgen] association_chromosomes = chromosome_genotypes
    Array[File] split_manifests = SplitGenotypes.manifest
    Array[Array[File]] split_plink_logs = SplitGenotypes.plink_logs
    Array[File] split_task_logs = SplitGenotypes.task_log
  }

  meta {
    description: "Fit a whole-genome REGENIE model, split association genotypes by chromosome, and test variants in Terra."
  }
  parameter_meta {
    step1_genotypes: "Whole-genome PGEN, uncompressed PVAR, and PSAM files."
    step1_variant_ids: "One selected common-variant ID per line, with no header. Select variants before this workflow."
    association_genotypes: "One or more PGEN file sets. Each set is split by chromosome before TestVariants. The same whole-genome files can be supplied for both steps."
    phenotypes: "Whitespace-delimited table with a header that starts with FID IID. Use NA for missing traits."
    phenotype_columns: "Optional comma-separated phenotype names. Omit or leave blank to use all columns after FID IID. All selected traits must use the same trait type."
    phenotype_batch_size: "Traits per association job. Default 50. Set to 0 to use all selected traits in one batch."
    batch_step1: "False: fit all selected traits together once. True: fit one model job per phenotype batch. Each job still uses whole-genome common SNPs."
    trait_type: "quantitative or binary."
    docker_image: "Container image string. Defaults to ghcr.io/aou-multiomics-analysis/regenie:latest. Supply a tag or digest to override it."
    covariates: "Optional table with a header that starts with FID IID."
    covariate_columns: "Optional comma-separated covariate names. If absent, use all covariate columns."
    categorical_covariates: "Optional comma-separated names of selected categorical covariates."
    binary_encoding: "01: 0=control and 1=case. 12: 1=control and 2=case."
    step1_extra_disk_gb: "Disk space in addition to localized input size. Include low-memory temporary files and predictions."
    step2_extra_disk_gb: "Disk space in addition to genotype and prediction input size. Include association outputs."
    split_extra_disk_gb: "Disk space in addition to twice the genotype input size. The split task retains localized inputs and chromosome outputs."
  }
}

task PrepareBatches {
  input {
    File phenotypes
    String? phenotype_columns
    Int batch_size = 50
    Boolean batch_step1 = false
    String docker_image
    Int preemptible_attempts = 0
    Int max_retries = 1
  }
  Int disk_gb = ceil(size(phenotypes, "GB")) + 1
  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start phenotype batch task"
    python3 /opt/regenie/scripts/prepare_batches.py \
      --phenotypes '~{sub(phenotypes, "'", "'\"'\"'")}' \
      ~{if defined(phenotype_columns) then "--phenotype-columns '" + sub(select_first([phenotype_columns]), "'", "'\"'\"'") + "'" else ""} \
      --batch-size ~{batch_size} \
      ~{if batch_step1 then "--batch-step1" else ""} \
      2>&1 | tee task.log
    echo "[$(date -u +%FT%TZ)] Phenotype batch task completed"
  >>>
  output {
    String selected_phenotypes = read_string("selected_phenotypes.txt")
    Array[PhenotypeBatch] phenotype_batches = read_json("phenotype_batches.json")
    Array[PhenotypeBatch] model_batches = read_json("model_batches.json")
    File phenotype_batch_manifest = "phenotype_batches.json"
    File model_batch_manifest = "model_batches.json"
    File task_log = "task.log"
  }
  runtime {
    docker: docker_image
    cpu: 1
    memory: "1 GB"
    disks: "local-disk ${disk_gb} SSD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}

task CollectPredictions {
  input {
    Array[File] prediction_lists
    Array[File] prediction_files
    Array[Int] prediction_file_counts
    String phenotype_columns
    Int batch_size = 50
    String docker_image
    Int preemptible_attempts = 0
    Int max_retries = 1
  }
  Int disk_gb = ceil(2 * size(prediction_files, "GB") + size(prediction_lists, "GB")) + 5
  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start prediction collection task"
    # Render localized File arrays at command time. Model outputs use fixed names.
    cat > prediction_lists.txt <<'REGENIE_LISTS_4c76c013733743c0a53d71e241f942c4'
    ~{sep="\n" prediction_lists}
    REGENIE_LISTS_4c76c013733743c0a53d71e241f942c4
    cat > prediction_files.txt <<'REGENIE_LOCO_a0e9c96ca4a44308ab087c14c09d9638'
    ~{sep="\n" prediction_files}
    REGENIE_LOCO_a0e9c96ca4a44308ab087c14c09d9638
    python3 /opt/regenie/scripts/collect_predictions.py \
      --prediction-lists prediction_lists.txt --prediction-files prediction_files.txt \
      --prediction-file-counts '~{sep="," prediction_file_counts}' \
      --phenotype-columns '~{sub(phenotype_columns, "'", "'\"'\"'")}' \
      --batch-size ~{batch_size} \
      2>&1 | tee task.log
    echo "[$(date -u +%FT%TZ)] Prediction collection task completed"
  >>>
  output {
    File prediction_list = "results/all_pred.list"
    Array[File] collected_prediction_files = glob("results/predictions/*.loco")
    Array[File] batch_prediction_lists = glob("results/batch_*_pred.list")
    Array[PredictionBatchIndices] batch_indices = read_json("results/prediction_batch_indices.json")
    File manifest = "results/prediction_batches.json"
    File task_log = "task.log"
  }
  runtime {
    docker: docker_image
    cpu: 1
    memory: "2 GB"
    disks: "local-disk ${disk_gb} SSD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}

task SplitGenotypes {
  input {
    File pgen
    File pvar
    File psam
    File? variant_ids
    String docker_image
    Int cpu = 4
    Int memory_gb = 16
    Int extra_disk_gb = 20
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  Int disk_gb = ceil(2 * size([pgen, pvar, psam], "GB") + size(variant_ids, "GB")) + extra_disk_gb
  Int plink_memory_mb = floor(memory_gb * 1024 * 0.8)

  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start genotype split task"
    python3 /opt/regenie/scripts/split_genotypes.py \
      --pgen '~{sub(pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(psam, "'", "'\"'\"'")}' \
      ~{if defined(variant_ids) then "--extract '" + sub(select_first([variant_ids]), "'", "'\"'\"'") + "'" else ""} \
      --threads ~{cpu} --memory-mb ~{plink_memory_mb} \
      2>&1 | tee task.log
    echo "[$(date -u +%FT%TZ)] Genotype split task completed"
  >>>

  output {
    # Explicit glob outputs let the Google backend plan file uploads before
    # the task runs. Groups preserve numeric chromosome order, then X/Y/XY/MT.
    Array[File] pgen_files = flatten([
      glob("chromosomes/chr[1-9].pgen"), glob("chromosomes/chr[12][0-9].pgen"),
      glob("chromosomes/chrX.pgen"), glob("chromosomes/chrY.pgen"),
      glob("chromosomes/chrXY.pgen"), glob("chromosomes/chrMT.pgen")
    ])
    Array[File] pvar_files = flatten([
      glob("chromosomes/chr[1-9].pvar"), glob("chromosomes/chr[12][0-9].pvar"),
      glob("chromosomes/chrX.pvar"), glob("chromosomes/chrY.pvar"),
      glob("chromosomes/chrXY.pvar"), glob("chromosomes/chrMT.pvar")
    ])
    Array[File] psam_files = flatten([
      glob("chromosomes/chr[1-9].psam"), glob("chromosomes/chr[12][0-9].psam"),
      glob("chromosomes/chrX.psam"), glob("chromosomes/chrY.psam"),
      glob("chromosomes/chrXY.psam"), glob("chromosomes/chrMT.psam")
    ])
    File manifest = "chromosome_genotypes.json"
    Array[File] plink_logs = glob("chromosomes/*.log")
    File task_log = "task.log"
  }
  runtime {
    docker: docker_image
    cpu: cpu
    memory: "${memory_gb} GB"
    disks: "local-disk ${disk_gb} SSD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}

task FitModel {
  input {
    File pgen
    File pvar
    File psam
    File variant_ids
    File phenotypes
    String? phenotype_columns
    String trait_type
    String docker_image
    File? covariates
    String? covariate_columns
    String? categorical_covariates
    File? keep_samples
    File? remove_samples
    String binary_encoding = "01"
    Boolean inverse_normal_transform = false
    Int min_case_count = 10
    Int cv_folds = 5
    Boolean lowmem = true
    Int block_size = 1000
    Int cpu = 8
    Int memory_gb = 32
    Int extra_disk_gb = 100
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  Int disk_gb = ceil(size([pgen, pvar, psam, phenotypes, variant_ids], "GB") + size(covariates, "GB")) + extra_disk_gb

  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start REGENIE model task"
    python3 /opt/regenie/scripts/run_regenie.py step1 \
      --pgen '~{sub(pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(psam, "'", "'\"'\"'")}' \
      --phenotypes '~{sub(phenotypes, "'", "'\"'\"'")}' \
      ~{if defined(phenotype_columns) then "--phenotype-columns '" + sub(select_first([phenotype_columns]), "'", "'\"'\"'") + "'" else ""} \
      --trait-type '~{sub(trait_type, "'", "'\"'\"'")}' \
      --extract '~{sub(variant_ids, "'", "'\"'\"'")}' \
      --binary-encoding '~{sub(binary_encoding, "'", "'\"'\"'")}' \
      --threads ~{cpu} --block-size ~{block_size} --cv-folds ~{cv_folds} \
      --min-case-count ~{min_case_count} \
      ~{if lowmem then "--lowmem" else ""} \
      ~{if inverse_normal_transform then "--rint" else ""} \
      ~{if defined(covariates) then "--covariates '" + sub(select_first([covariates]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(covariate_columns) then "--covariate-columns '" + sub(select_first([covariate_columns]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(categorical_covariates) then "--categorical-covariates '" + sub(select_first([categorical_covariates]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(keep_samples) then "--keep '" + sub(select_first([keep_samples]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(remove_samples) then "--remove '" + sub(select_first([remove_samples]), "'", "'\"'\"'") + "'" else ""} \
      2>&1 | tee task.log
    echo "[$(date -u +%FT%TZ)] REGENIE model task completed"
  >>>

  output {
    File prediction_list = "results/fit_pred.list"
    Array[File] prediction_files = glob("results/fit_*.loco")
    Int prediction_file_count = length(glob("results/fit_*.loco"))
    File regenie_log = "results/fit.log"
    File task_log = "task.log"
  }
  runtime {
    docker: docker_image
    cpu: cpu
    memory: "${memory_gb} GB"
    disks: "local-disk ${disk_gb} SSD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}

task TestVariants {
  input {
    File pgen
    File pvar
    File psam
    File phenotypes
    String? phenotype_columns
    String trait_type
    String docker_image
    File prediction_list
    Array[File] prediction_files
    File? covariates
    String? covariate_columns
    String? categorical_covariates
    File? keep_samples
    File? remove_samples
    File? variant_ids
    String binary_encoding = "01"
    Boolean inverse_normal_transform = false
    Boolean approximate_firth = true
    Float firth_threshold = 0.01
    Float min_mac = 5.0
    Int min_case_count = 10
    Int block_size = 400
    Int cpu = 4
    Int memory_gb = 16
    Int extra_disk_gb = 50
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  Int disk_gb = ceil(size([pgen, pvar, psam, phenotypes, prediction_list], "GB") + size(prediction_files, "GB") + size(covariates, "GB")) + extra_disk_gb

  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start REGENIE association task"
    # Render File paths after localization. The quoted delimiter prevents shell expansion.
    # Files come from FitModel with fixed fit_<number>.loco names.
    cat > prediction_files.txt <<'REGENIE_FILES_70ce88a66f5a40a4914559e3c5b3f6cd'
    ~{sep="\n" prediction_files}
    REGENIE_FILES_70ce88a66f5a40a4914559e3c5b3f6cd
    python3 /opt/regenie/scripts/run_regenie.py step2 \
      --pgen '~{sub(pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(psam, "'", "'\"'\"'")}' \
      --phenotypes '~{sub(phenotypes, "'", "'\"'\"'")}' \
      ~{if defined(phenotype_columns) then "--phenotype-columns '" + sub(select_first([phenotype_columns]), "'", "'\"'\"'") + "'" else ""} \
      --trait-type '~{sub(trait_type, "'", "'\"'\"'")}' \
      --prediction-list '~{sub(prediction_list, "'", "'\"'\"'")}' \
      --prediction-files prediction_files.txt \
      --binary-encoding '~{sub(binary_encoding, "'", "'\"'\"'")}' \
      --threads ~{cpu} --block-size ~{block_size} --min-mac ~{min_mac} \
      --min-case-count ~{min_case_count} --firth-threshold ~{firth_threshold} \
      ~{if inverse_normal_transform then "--rint" else ""} \
      ~{if approximate_firth then "" else "--no-firth"} \
      ~{if defined(covariates) then "--covariates '" + sub(select_first([covariates]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(covariate_columns) then "--covariate-columns '" + sub(select_first([covariate_columns]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(categorical_covariates) then "--categorical-covariates '" + sub(select_first([categorical_covariates]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(keep_samples) then "--keep '" + sub(select_first([keep_samples]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(remove_samples) then "--remove '" + sub(select_first([remove_samples]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(variant_ids) then "--extract '" + sub(select_first([variant_ids]), "'", "'\"'\"'") + "'" else ""} \
      2>&1 | tee task.log
    echo "[$(date -u +%FT%TZ)] REGENIE association task completed"
  >>>

  output {
    Array[File] associations = glob("results/assoc*.regenie.gz")
    File regenie_log = "results/assoc.log"
    File task_log = "task.log"
  }
  runtime {
    docker: docker_image
    cpu: cpu
    memory: "${memory_gb} GB"
    disks: "local-disk ${disk_gb} SSD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}

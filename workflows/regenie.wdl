version 1.0

struct PgenFiles {
  File pgen
  File pvar
  File psam
}

workflow Regenie {
  input {
    PgenFiles step1_genotypes
    File step1_variant_ids
    Array[PgenFiles]+ association_genotypes
    File phenotypes
    String phenotype_columns
    String trait_type = "quantitative"
    String docker_image

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
    Int step1_memory_gb = 32
    Int step2_memory_gb = 16
    Int step1_extra_disk_gb = 100
    Int step2_extra_disk_gb = 50
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  call FitModel {
    input:
      genotypes = step1_genotypes,
      variant_ids = step1_variant_ids,
      phenotypes = phenotypes,
      phenotype_columns = phenotype_columns,
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

  scatter (genotype_set in association_genotypes) {
    call TestVariants {
      input:
        genotypes = genotype_set,
        phenotypes = phenotypes,
        phenotype_columns = phenotype_columns,
        trait_type = trait_type,
        docker_image = docker_image,
        prediction_list = FitModel.prediction_list,
        prediction_files = FitModel.prediction_files,
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
    File prediction_list = FitModel.prediction_list
    Array[File] prediction_files = FitModel.prediction_files
    File step1_log = FitModel.regenie_log
    File step1_task_log = FitModel.task_log
    Array[Array[File]] associations_by_shard = TestVariants.associations
    Array[File] associations = flatten(TestVariants.associations)
    Array[File] step2_logs = TestVariants.regenie_log
    Array[File] step2_task_logs = TestVariants.task_log
  }

  meta {
    description: "Fit a whole-genome REGENIE model and test PGEN variants in Terra. Use one trait type per run."
  }
  parameter_meta {
    step1_genotypes: "Whole-genome PGEN, uncompressed PVAR, and PSAM files."
    step1_variant_ids: "One selected common-variant ID per line, with no header. Select variants before this workflow."
    association_genotypes: "One or more PGEN file sets. Each set can contain one chromosome or the whole genome."
    phenotypes: "Whitespace-delimited table with a header that starts with FID IID. Use NA for missing traits."
    phenotype_columns: "Comma-separated phenotype names. All selected traits must use the same trait type."
    trait_type: "quantitative or binary."
    docker_image: "Publicly readable container image. Use the digest from the GitHub Actions image build."
    covariates: "Optional table with a header that starts with FID IID."
    covariate_columns: "Optional comma-separated covariate names. If absent, use all covariate columns."
    categorical_covariates: "Optional comma-separated names of selected categorical covariates."
    binary_encoding: "01: 0=control and 1=case. 12: 1=control and 2=case."
    step1_extra_disk_gb: "Disk space in addition to localized input size. Include low-memory temporary files and predictions."
    step2_extra_disk_gb: "Disk space in addition to genotype and prediction input size. Include association outputs."
  }
}

task FitModel {
  input {
    PgenFiles genotypes
    File variant_ids
    File phenotypes
    String phenotype_columns
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

  Int disk_gb = ceil(size([genotypes.pgen, genotypes.pvar, genotypes.psam, phenotypes, variant_ids], "GB") + size(covariates, "GB")) + extra_disk_gb

  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start REGENIE model task"
    python3 /opt/regenie/scripts/run_regenie.py step1 \
      --pgen '~{sub(genotypes.pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(genotypes.pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(genotypes.psam, "'", "'\"'\"'")}' \
      --phenotypes '~{sub(phenotypes, "'", "'\"'\"'")}' \
      --phenotype-columns '~{sub(phenotype_columns, "'", "'\"'\"'")}' \
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
    PgenFiles genotypes
    File phenotypes
    String phenotype_columns
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

  Int disk_gb = ceil(size([genotypes.pgen, genotypes.pvar, genotypes.psam, phenotypes, prediction_list], "GB") + size(prediction_files, "GB") + size(covariates, "GB")) + extra_disk_gb

  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start REGENIE association task"
    # Render File paths after localization. The quoted delimiter prevents shell expansion.
    # Files come from FitModel with fixed fit_<number>.loco names.
    cat > prediction_files.txt <<'REGENIE_FILES_70ce88a66f5a40a4914559e3c5b3f6cd'
    ~{sep="\n" prediction_files}
    REGENIE_FILES_70ce88a66f5a40a4914559e3c5b3f6cd
    python3 /opt/regenie/scripts/run_regenie.py step2 \
      --pgen '~{sub(genotypes.pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(genotypes.pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(genotypes.psam, "'", "'\"'\"'")}' \
      --phenotypes '~{sub(phenotypes, "'", "'\"'\"'")}' \
      --phenotype-columns '~{sub(phenotype_columns, "'", "'\"'\"'")}' \
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

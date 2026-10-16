version 1.0

import "glow_regenie_dataproc.wdl" as shared

struct GlowLocalPgenFiles { File pgen File pvar File psam }
struct GlowLocalBgenFiles { File bgen File sample File variants Int source_index }

workflow GlowRegenie {
  input {
    GlowLocalPgenFiles step1_genotypes
    File step1_variant_ids
    Array[GlowLocalPgenFiles]+ association_genotypes
    File phenotypes
    File? covariates
    File? keep_samples
    File? remove_samples
    File? association_variant_ids
    String docker_image = "ghcr.io/aou-multiomics-analysis/regenie-glow@sha256:c2d18e1db4049aa2c3fe9a9833037cfddad441d3c2074e04b9c9fb6fb00956b2"
    String trait_type
    String? phenotype_columns
    String? covariate_columns
    String? categorical_covariates
    String binary_encoding = "01"
    String missing_data_policy = "error"
    Int batch_size = 50
    Int variants_per_block = 1000
    Int sample_block_count = 10
    Float min_mac = 5.0
    Int min_cases = 10
    Int min_controls = 10
    Float firth_threshold = 0.01
    Int seed = 1
    Boolean inverse_normal_transform = false
    Int analysis_cpu = 64
    Int analysis_memory_gb = 256
    Int analysis_extra_disk_gb = 1000
    Int spark_driver_memory_gb = 64
    Int spark_parallelism = 100
    Int spark_shuffle_partitions = 100
    Int prepare_cpu = 1
    Int prepare_memory_gb = 8
    Int prepare_extra_disk_gb = 10
    Int convert_cpu = 4
    Int convert_memory_gb = 16
    Int convert_extra_disk_gb = 100
  }

  call shared.PrepareGlowInputs {
    input:
      model_psam = step1_genotypes.psam, phenotypes = phenotypes, covariates = covariates,
      keep_samples = keep_samples, remove_samples = remove_samples,
      trait_type = trait_type,
      phenotype_columns = phenotype_columns,
      covariate_columns = covariate_columns,
      categorical_covariates = categorical_covariates,
      binary_encoding = binary_encoding,
      missing_data_policy = missing_data_policy,
      batch_size = batch_size,
      variants_per_block = variants_per_block,
      sample_block_count = sample_block_count,
      min_mac = min_mac,
      min_cases = min_cases,
      min_controls = min_controls,
      firth_threshold = firth_threshold,
      seed = seed,
      inverse_normal_transform = inverse_normal_transform,
      docker_image = docker_image, cpu = prepare_cpu, memory_gb = prepare_memory_gb,
      extra_disk_gb = prepare_extra_disk_gb
  }

  call shared.ConvertGlowGenotypes as ConvertGlowModelGenotypes {
    input:
      pgen = step1_genotypes.pgen, pvar = step1_genotypes.pvar, psam = step1_genotypes.psam,
      sample_map = PrepareGlowInputs.sample_map, variant_ids = step1_variant_ids,
      source_index = 0, role = "model", docker_image = docker_image,
      cpu = convert_cpu, memory_gb = convert_memory_gb, extra_disk_gb = convert_extra_disk_gb
  }

  scatter (source_index in range(length(association_genotypes))) {
    call shared.ConvertGlowGenotypes as ConvertGlowAssociationGenotypes {
      input:
        pgen = association_genotypes[source_index].pgen, pvar = association_genotypes[source_index].pvar,
        psam = association_genotypes[source_index].psam, sample_map = PrepareGlowInputs.sample_map,
        variant_ids = association_variant_ids, source_index = source_index, role = "association",
        docker_image = docker_image, cpu = convert_cpu, memory_gb = convert_memory_gb, extra_disk_gb = convert_extra_disk_gb
    }
    GlowLocalBgenFiles association_set = object {
      bgen: ConvertGlowAssociationGenotypes.bgen, sample: ConvertGlowAssociationGenotypes.sample,
      variants: ConvertGlowAssociationGenotypes.variants, source_index: source_index
    }
  }

  call RunGlowLocal {
    input:
      model_bgen = ConvertGlowModelGenotypes.bgen, model_sample = ConvertGlowModelGenotypes.sample,
      model_variants = ConvertGlowModelGenotypes.variants,
      association_bgens = ConvertGlowAssociationGenotypes.bgen,
      association_samples = ConvertGlowAssociationGenotypes.sample,
      association_variants = ConvertGlowAssociationGenotypes.variants,
      sample_map = PrepareGlowInputs.sample_map, phenotypes = PrepareGlowInputs.prepared_phenotypes,
      covariates = PrepareGlowInputs.prepared_covariates, phenotype_batches = PrepareGlowInputs.phenotype_batches,
      trait_type = trait_type,
      phenotype_columns = phenotype_columns,
      covariate_columns = covariate_columns,
      categorical_covariates = categorical_covariates,
      binary_encoding = binary_encoding,
      missing_data_policy = missing_data_policy,
      batch_size = batch_size,
      variants_per_block = variants_per_block,
      sample_block_count = sample_block_count,
      min_mac = min_mac,
      min_cases = min_cases,
      min_controls = min_controls,
      firth_threshold = firth_threshold,
      seed = seed,
      inverse_normal_transform = inverse_normal_transform,
      driver_memory_gb = spark_driver_memory_gb,
      parallelism = spark_parallelism,
      shuffle_partitions = spark_shuffle_partitions,
      docker_image = docker_image, cpu = analysis_cpu, memory_gb = analysis_memory_gb,
      extra_disk_gb = analysis_extra_disk_gb
  }

  output {
    Array[File] associations = RunGlowLocal.associations
    Array[File] predictions = RunGlowLocal.predictions
    Array[File] metrics = RunGlowLocal.metrics
    Array[File] application_logs = RunGlowLocal.application_logs
    File output_manifest = RunGlowLocal.output_manifest
    File execution_report = RunGlowLocal.execution_report
    File analysis_log = RunGlowLocal.task_log
    File preparation_report = PrepareGlowInputs.report
    File preparation_log = PrepareGlowInputs.task_log
    File sample_map = PrepareGlowInputs.sample_map
    Array[GlowLocalBgenFiles] association_sets = association_set
    File model_bgen = ConvertGlowModelGenotypes.bgen
    File model_sample = ConvertGlowModelGenotypes.sample
    File model_variant_metadata = ConvertGlowModelGenotypes.variants
    Array[File] conversion_logs = flatten([[ConvertGlowModelGenotypes.plink_log], ConvertGlowAssociationGenotypes.plink_log])
    Array[File] conversion_task_logs = flatten([[ConvertGlowModelGenotypes.task_log], ConvertGlowAssociationGenotypes.task_log])
  }

  meta {
    description: "Glow with Spark on one Terra VM. WDL 1.0. Not tested on Terra."
  }
}

task RunGlowLocal {
  input {
    File model_bgen
    File model_sample
    File model_variants
    Array[File]+ association_bgens
    Array[File]+ association_samples
    Array[File]+ association_variants
    File sample_map
    File phenotypes
    File covariates
    File phenotype_batches
    String trait_type
    String? phenotype_columns
    String? covariate_columns
    String? categorical_covariates
    String binary_encoding = "01"
    String missing_data_policy = "error"
    Int batch_size = 50
    Int variants_per_block = 1000
    Int sample_block_count = 10
    Float min_mac = 5.0
    Int min_cases = 10
    Int min_controls = 10
    Float firth_threshold = 0.01
    Int seed = 1
    Boolean inverse_normal_transform = false
    Int driver_memory_gb = 64
    Int parallelism = 100
    Int shuffle_partitions = 100
    String docker_image
    Int cpu = 64
    Int memory_gb = 256
    Int extra_disk_gb = 1000
  }
  Int disk_gb = ceil(size([model_bgen, model_sample, model_variants, sample_map, phenotypes, covariates, phenotype_batches], "GB") + size(association_bgens, "GB") + size(association_samples, "GB") + size(association_variants, "GB")) + extra_disk_gb
  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start Glow local analysis"
    cat > association_bgens.txt <<'GLOW_BGEN_8de2c5e73b704d948777'
    ~{sep="\n" association_bgens}
    GLOW_BGEN_8de2c5e73b704d948777
    cat > association_samples.txt <<'GLOW_SAMPLE_07536efbe14c4cae982e'
    ~{sep="\n" association_samples}
    GLOW_SAMPLE_07536efbe14c4cae982e
    cat > association_variants.txt <<'GLOW_VARIANTS_c9a7b2de62fe427e9f56'
    ~{sep="\n" association_variants}
    GLOW_VARIANTS_c9a7b2de62fe427e9f56
    python3 -m scripts.glow.local_runner \
      --model-bgen '~{sub(model_bgen, "'", "'\"'\"'")}' \
      --model-sample '~{sub(model_sample, "'", "'\"'\"'")}' \
      --model-variants '~{sub(model_variants, "'", "'\"'\"'")}' \
      --sample-map '~{sub(sample_map, "'", "'\"'\"'")}' \
      --phenotypes '~{sub(phenotypes, "'", "'\"'\"'")}' \
      --covariates '~{sub(covariates, "'", "'\"'\"'")}' \
      --phenotype-batches '~{sub(phenotype_batches, "'", "'\"'\"'")}' \
      --association-bgen-list association_bgens.txt \
      --association-sample-list association_samples.txt \
      --association-variants-list association_variants.txt \
      --cpu ~{cpu} \
      --memory-gb ~{memory_gb} \
      --driver-memory-gb ~{driver_memory_gb} \
      --parallelism ~{parallelism} \
      --shuffle-partitions ~{shuffle_partitions} \
      --trait-type '~{sub(trait_type, "'", "'\"'\"'")}' \
      ~{if defined(phenotype_columns) then "--phenotype-columns '" + sub(select_first([phenotype_columns]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(covariate_columns) then "--covariate-columns '" + sub(select_first([covariate_columns]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(categorical_covariates) then "--categorical-covariates '" + sub(select_first([categorical_covariates]), "'", "'\"'\"'") + "'" else ""} \
      --binary-encoding '~{sub(binary_encoding, "'", "'\"'\"'")}' \
      --missing-data-policy '~{sub(missing_data_policy, "'", "'\"'\"'")}' \
      --batch-size ~{batch_size} \
      --variants-per-block ~{variants_per_block} \
      --sample-block-count ~{sample_block_count} \
      --min-mac ~{min_mac} \
      --min-cases ~{min_cases} \
      --min-controls ~{min_controls} \
      --firth-threshold ~{firth_threshold} \
      --seed ~{seed} \
      ~{if inverse_normal_transform then "--inverse-normal-transform" else ""} \
      --out-dir results \
      2>&1 | tee task.log
    echo "[$(date -u +%FT%TZ)] Glow local analysis completed"
  >>>
  output {
    Array[File] associations = glob("results/associations/*")
    Array[File] predictions = glob("results/predictions/*")
    Array[File] metrics = glob("results/metrics/*")
    Array[File] application_logs = glob("results/logs/*")
    File output_manifest = "results/output_manifest.json"
    File execution_report = "results/execution.json"
    File task_log = "task.log"
  }
  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: 0
    maxRetries: 1
  }
}

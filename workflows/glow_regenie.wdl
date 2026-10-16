version 1.0

struct GlowPgenFiles { File pgen File pvar File psam }
struct GlowBgenFiles { File bgen File sample File variants Int source_index }

workflow GlowRegenie {
  input {
    GlowPgenFiles step1_genotypes
    File step1_variant_ids
    Array[GlowPgenFiles]+ association_genotypes
    File phenotypes
    File? covariates
    File? keep_samples
    File? remove_samples
    File? association_variant_ids
    String docker_image
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
    String project
    String region
    String subnetwork_uri
    String service_account
    String staging_bucket
    String staging_prefix
    String run_id
    String spark_image
    Int driver_cores = 4
    Int driver_memory_mb = 8192
    Int driver_overhead_mb = 4096
    Int executor_cores = 4
    Int executor_memory_mb = 8192
    Int executor_overhead_mb = 4096
    Int min_executors = 2
    Int max_executors = 4
    Int driver_disk_gb = 250
    Int executor_disk_gb = 250
    Int ttl_seconds = 21600
    Int poll_seconds = 30
    Int prepare_cpu = 1
    Int prepare_memory_gb = 8
    Int prepare_extra_disk_gb = 10
    Int convert_cpu = 4
    Int convert_memory_gb = 16
    Int convert_extra_disk_gb = 100
    Int controller_cpu = 1
    Int controller_memory_gb = 4
    Int controller_extra_disk_gb = 100
  }

  call PrepareGlowInputs {
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

  call ConvertGlowGenotypes as ConvertGlowModelGenotypes {
    input:
      pgen = step1_genotypes.pgen, pvar = step1_genotypes.pvar, psam = step1_genotypes.psam,
      sample_map = PrepareGlowInputs.sample_map, variant_ids = step1_variant_ids,
      source_index = 0, role = "model", docker_image = docker_image,
      cpu = convert_cpu, memory_gb = convert_memory_gb, extra_disk_gb = convert_extra_disk_gb
  }

  scatter (source_index in range(length(association_genotypes))) {
    call ConvertGlowGenotypes as ConvertGlowAssociationGenotypes {
      input:
        pgen = association_genotypes[source_index].pgen, pvar = association_genotypes[source_index].pvar,
        psam = association_genotypes[source_index].psam, sample_map = PrepareGlowInputs.sample_map,
        variant_ids = association_variant_ids, source_index = source_index, role = "association",
        docker_image = docker_image, cpu = convert_cpu, memory_gb = convert_memory_gb, extra_disk_gb = convert_extra_disk_gb
    }
    GlowBgenFiles association_set = object {
      bgen: ConvertGlowAssociationGenotypes.bgen, sample: ConvertGlowAssociationGenotypes.sample,
      variants: ConvertGlowAssociationGenotypes.variants, source_index: source_index
    }
  }

  call RunGlowAnalysis {
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
      project = project,
      region = region,
      subnetwork_uri = subnetwork_uri,
      service_account = service_account,
      staging_bucket = staging_bucket,
      staging_prefix = staging_prefix,
      run_id = run_id,
      spark_image = spark_image,
      driver_cores = driver_cores,
      driver_memory_mb = driver_memory_mb,
      driver_overhead_mb = driver_overhead_mb,
      executor_cores = executor_cores,
      executor_memory_mb = executor_memory_mb,
      executor_overhead_mb = executor_overhead_mb,
      min_executors = min_executors,
      max_executors = max_executors,
      driver_disk_gb = driver_disk_gb,
      executor_disk_gb = executor_disk_gb,
      ttl_seconds = ttl_seconds,
      poll_seconds = poll_seconds,
      docker_image = docker_image, cpu = controller_cpu, memory_gb = controller_memory_gb,
      extra_disk_gb = controller_extra_disk_gb
  }

  output {
    Array[File] associations = RunGlowAnalysis.associations
    Array[File] predictions = RunGlowAnalysis.predictions
    Array[File] metrics = RunGlowAnalysis.metrics
    Array[File] application_logs = RunGlowAnalysis.application_logs
    File output_manifest = RunGlowAnalysis.output_manifest
    File execution_report = RunGlowAnalysis.execution_report
    File controller_log = RunGlowAnalysis.task_log
    File preparation_report = PrepareGlowInputs.report
    File preparation_log = PrepareGlowInputs.task_log
    File sample_map = PrepareGlowInputs.sample_map
    Array[GlowBgenFiles] association_sets = association_set
    File model_bgen = ConvertGlowModelGenotypes.bgen
    File model_sample = ConvertGlowModelGenotypes.sample
    File model_variant_metadata = ConvertGlowModelGenotypes.variants
    Array[File] conversion_logs = flatten([[ConvertGlowModelGenotypes.task_log], ConvertGlowAssociationGenotypes.task_log])
  }

  meta {
    description: "WDL 1.0 controller for Glow on Dataproc Serverless. Not tested on Terra. Disable call caching."
  }
}

task PrepareGlowInputs {
  input {
    File model_psam
    File phenotypes
    File? covariates
    File? keep_samples
    File? remove_samples
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
    String docker_image
    Int cpu = 1
    Int memory_gb = 8
    Int extra_disk_gb = 10
  }
  Int disk_gb = ceil(size([model_psam, phenotypes], "GB") + size(covariates, "GB") + size(keep_samples, "GB") + size(remove_samples, "GB")) + extra_disk_gb
  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start Glow cohort preparation"
    python3 -m scripts.glow.prepare \
      --psam '~{sub(model_psam, "'", "'\"'\"'")}' \
      --phenotypes '~{sub(phenotypes, "'", "'\"'\"'")}' \
      ~{if defined(covariates) then "--covariates '" + sub(select_first([covariates]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(keep_samples) then "--keep '" + sub(select_first([keep_samples]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(remove_samples) then "--remove '" + sub(select_first([remove_samples]), "'", "'\"'\"'") + "'" else ""} \
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
      --out-dir prepared \
      2>&1 | tee task.log
    echo "[$(date -u +%FT%TZ)] Glow cohort preparation completed"
  >>>
  output {
    File sample_map = "prepared/sample_map.tsv"
    File prepared_phenotypes = "prepared/phenotypes.tsv"
    File prepared_covariates = "prepared/covariates.tsv"
    File phenotype_batches = "prepared/phenotype_batches.tsv"
    File report = "prepared/preparation.json"
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

task ConvertGlowGenotypes {
  input {
    File pgen
    File pvar
    File psam
    File sample_map
    File? variant_ids
    Int source_index
    String role
    String docker_image
    Int cpu = 4
    Int memory_gb = 16
    Int extra_disk_gb = 100
  }
  Int disk_gb = ceil(size([pgen, pvar, psam, sample_map], "GB") + size(variant_ids, "GB")) + extra_disk_gb
  Int plink_memory_mb = memory_gb * 750
  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start Glow genotype conversion"
    python3 -m scripts.glow.convert \
      --pgen '~{sub(pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(psam, "'", "'\"'\"'")}' \
      --sample-map '~{sub(sample_map, "'", "'\"'\"'")}' \
      ~{if defined(variant_ids) then "--variant-ids '" + sub(select_first([variant_ids]), "'", "'\"'\"'") + "'" else ""} \
      --source-index ~{source_index} \
      --role '~{sub(role, "'", "'\"'\"'")}' \
      --threads ~{cpu} \
      --memory-mb ~{plink_memory_mb} \
      --out-dir converted \
      2>&1 | tee task.log
    echo "[$(date -u +%FT%TZ)] Glow genotype conversion completed"
  >>>
  output {
    File bgen = "converted/genotypes.bgen"
    File sample = "converted/genotypes.sample"
    File variants = "converted/variants.tsv"
    File report = "converted/conversion.json"
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

task RunGlowAnalysis {
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
    String project
    String region
    String subnetwork_uri
    String service_account
    String staging_bucket
    String staging_prefix
    String run_id
    String spark_image
    Int driver_cores = 4
    Int driver_memory_mb = 8192
    Int driver_overhead_mb = 4096
    Int executor_cores = 4
    Int executor_memory_mb = 8192
    Int executor_overhead_mb = 4096
    Int min_executors = 2
    Int max_executors = 4
    Int driver_disk_gb = 250
    Int executor_disk_gb = 250
    Int ttl_seconds = 21600
    Int poll_seconds = 30
    String docker_image
    Int cpu = 1
    Int memory_gb = 4
    Int extra_disk_gb = 100
  }
  Int disk_gb = ceil(size([model_bgen, model_sample, model_variants, sample_map, phenotypes, covariates, phenotype_batches], "GB") + size(association_bgens, "GB") + size(association_samples, "GB") + size(association_variants, "GB")) + extra_disk_gb
  command <<<
    set -euo pipefail
    echo "[$(date -u +%FT%TZ)] Start Glow remote controller"
    cat > association_bgens.txt <<'GLOW_BGEN_8de2c5e73b704d948777'
    ~{sep="\n" association_bgens}
    GLOW_BGEN_8de2c5e73b704d948777
    cat > association_samples.txt <<'GLOW_SAMPLE_07536efbe14c4cae982e'
    ~{sep="\n" association_samples}
    GLOW_SAMPLE_07536efbe14c4cae982e
    cat > association_variants.txt <<'GLOW_VARIANTS_c9a7b2de62fe427e9f56'
    ~{sep="\n" association_variants}
    GLOW_VARIANTS_c9a7b2de62fe427e9f56
    python3 -m scripts.glow.controller \
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
      --terra-image '~{sub(docker_image, "'", "'\"'\"'")}' \
      --project '~{sub(project, "'", "'\"'\"'")}' \
      --region '~{sub(region, "'", "'\"'\"'")}' \
      --subnetwork-uri '~{sub(subnetwork_uri, "'", "'\"'\"'")}' \
      --service-account '~{sub(service_account, "'", "'\"'\"'")}' \
      --staging-bucket '~{sub(staging_bucket, "'", "'\"'\"'")}' \
      --staging-prefix '~{sub(staging_prefix, "'", "'\"'\"'")}' \
      --run-id '~{sub(run_id, "'", "'\"'\"'")}' \
      --spark-image '~{sub(spark_image, "'", "'\"'\"'")}' \
      --driver-cores ~{driver_cores} \
      --driver-memory-mb ~{driver_memory_mb} \
      --driver-overhead-mb ~{driver_overhead_mb} \
      --executor-cores ~{executor_cores} \
      --executor-memory-mb ~{executor_memory_mb} \
      --executor-overhead-mb ~{executor_overhead_mb} \
      --min-executors ~{min_executors} \
      --max-executors ~{max_executors} \
      --driver-disk-gb ~{driver_disk_gb} \
      --executor-disk-gb ~{executor_disk_gb} \
      --ttl-seconds ~{ttl_seconds} \
      --poll-seconds ~{poll_seconds} \
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
    echo "[$(date -u +%FT%TZ)] Glow remote controller completed"
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

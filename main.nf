nextflow.enable.dsl=2

/*
 * APGAP influenza starter pipeline
 * Inputs:
 *   - Easy mode:   --reads_dir "/path/to/fastqs"
 *   - Advanced:    --reads     "/path/to/fastqs/*_{R1,R2}_001.fastq.gz"  (or similar)
 *
 * Workflow:
 *   Raw reads -> (optional primer trim) -> FastQC -> fastp -> FastQC(trimmed) -> MultiQC -> IRMA
 *
 * Primer trimming uses params:
 *   params.trim_primers (true/false)
 *   params.primers_5p_r1 (list of sequences)
 *   params.primers_5p_r2 (list of sequences)
 *
 * Put defaults in nextflow.config for shared use.
 */

process CUTADAPT_PRIMERS {
  tag "$sample_id"
  publishDir "${params.outdir}/trim/primers", mode: 'copy'
  cpus params.threads

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id),
          path("${sample_id}_R1.primertrim.fastq.gz"),
          path("${sample_id}_R2.primertrim.fastq.gz")

  container "quay.io/biocontainers/cutadapt:4.8--py310h4b81fae_0"

  script:
    def g_opts = (params.primers_5p_r1 ?: []).collect { p -> "-g ${p}" }.join(' ')
    def G_opts = (params.primers_5p_r2 ?: []).collect { p -> "-G ${p}" }.join(' ')
    def overlap = params.cutadapt_min_overlap ?: 10
    def errrate = params.cutadapt_error_rate  ?: 0.1

    """
    cutadapt \
      ${g_opts} ${G_opts} \
      --overlap ${overlap} \
      --error-rate ${errrate} \
      -j ${task.cpus} \
      -o ${sample_id}_R1.primertrim.fastq.gz \
      -p ${sample_id}_R2.primertrim.fastq.gz \
      ${r1} ${r2}
    """
}

process FASTQC_RAW {
  tag "$sample_id"
  publishDir "${params.outdir}/qc/fastqc/raw", mode: 'copy'
  cpus params.threads

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    path "*_fastqc.{zip,html}"

  container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

  """
  fastqc --threads ${task.cpus} ${r1} ${r2}
  """
}

process FASTP {
  tag "$sample_id"
  publishDir "${params.outdir}/trim/fastp", mode: 'copy'
  cpus params.threads

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id),
          path("${sample_id}_R1.trim.fastq.gz"),
          path("${sample_id}_R2.trim.fastq.gz"),
          path("${sample_id}.fastp.html"),
          path("${sample_id}.fastp.json")

  container "quay.io/biocontainers/fastp:0.23.2--h79da9fb_0"

  """
  fastp \
    -i ${r1} -I ${r2} \
    -o ${sample_id}_R1.trim.fastq.gz \
    -O ${sample_id}_R2.trim.fastq.gz \
    -h ${sample_id}.fastp.html \
    -j ${sample_id}.fastp.json \
    -w ${task.cpus}
  """
}

process FASTQC_TRIMMED {
  tag "$sample_id"
  publishDir "${params.outdir}/qc/fastqc/trimmed", mode: 'copy'
  cpus params.threads

  input:
    tuple val(sample_id), path(r1t), path(r2t)

  output:
    path "*_fastqc.{zip,html}"

  container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

  """
  fastqc --threads ${task.cpus} ${r1t} ${r2t}
  """
}

process MULTIQC {
  publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

  input:
    path qc_files

  output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

  container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"

  """
  multiqc -o . .
  """
}

process IRMA_RUN {
  tag "$sample_id"
  publishDir "${params.outdir}/assembly/irma", mode: 'copy'

  input:
    tuple val(sample_id), path(r1t), path(r2t)

  output:
    path "${sample_id}"

  container "ghcr.io/cdcgov/irma:v1.3.1"

  """
  IRMA ${params.irma_module} ${r1t} ${r2t} ${sample_id}
  """
}

workflow {

  /*
   * Build input pairs from either:
   *   --reads_dir  (easy; no glob knowledge needed)
   *   --reads      (advanced; explicit pattern)
   */
  def build_read_pairs = {
    if (params.reads_dir) {
      def d = params.reads_dir.toString().replaceAll(/\/$/, '')

      def candidates = [
        "${d}/*_{R1,R2}_001.fastq.gz",
        "${d}/*_{R1,R2}.fastq.gz",
        "${d}/*_{R1,R2}_001.fq.gz",
        "${d}/*_{R1,R2}.fq.gz"
      ]

      for (pat in candidates) {
        def files = file(pat)
        if (files && files.size() > 0) {
          log.info "Using reads_dir '${params.reads_dir}' with pattern: ${pat}"
          return Channel.fromFilePairs(pat, flat: true)
            .map { x -> tuple(x[0], x[1], x[2]) }
        }
      }

      error "No FASTQ pairs found in reads_dir='${params.reads_dir}'. Expected files like *_R1_001.fastq.gz and *_R2_001.fastq.gz (or *_R1.fastq.gz/_R2.fastq.gz)."
    }

    if (params.reads) {
      log.info "Using reads pattern: ${params.reads}"
      return Channel.fromFilePairs(params.reads, flat: true)
        .map { x -> tuple(x[0], x[1], x[2]) }
    }

    error "Provide either --reads_dir '/path/to/fastqs' (recommended) or --reads 'glob_pattern'."
  }

  read_pairs = build_read_pairs()

  // Optional primer trimming step (uses params.trim_primers and primer lists)
  reads_for_qc = read_pairs
  if (params.trim_primers) {
    reads_for_qc = CUTADAPT_PRIMERS(read_pairs)
  }

  // Raw QC (raw in the sense of "pre-fastp"; may already be primer-trimmed if enabled)
  raw_qc = FASTQC_RAW(reads_for_qc)

  // fastp trimming
  trimmed = FASTP(reads_for_qc)

  // trimmed_reads: (sid, r1t, r2t)
  trimmed_reads = trimmed.map { sid, r1t, r2t, html, json -> tuple(sid, r1t, r2t) }

  // fastp reports: html + json
  fastp_reports = trimmed
    .map { sid, r1t, r2t, html, json -> [html, json] }
    .flatten()

  // Trimmed QC
  trimmed_qc = FASTQC_TRIMMED(trimmed_reads)

  // MultiQC FIRST
  qc_files = raw_qc
    .mix(trimmed_qc)
    .mix(fastp_reports)
    .collect()

  mq = MULTIQC(qc_files)

  // ---- Gate IRMA on MultiQC completion, WITHOUT losing tuple structure ----
  done = mq.report.map { 1 }   // emits one "1" after MultiQC finishes

  trimmed_after_multiqc = trimmed_reads
    .combine(done)
    .map { it -> tuple(it[0], it[1], it[2]) }

  IRMA_RUN(trimmed_after_multiqc)
}

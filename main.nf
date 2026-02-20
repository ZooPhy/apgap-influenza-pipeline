nextflow.enable.dsl=2

/*
 * APGAP influenza starter pipeline
 * Raw reads -> FastQC -> fastp -> FastQC(trimmed) -> MultiQC -> IRMA
 */

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
  // fromFilePairs emits items like: [sample_id, r1, r2]
  read_pairs = Channel
    .fromFilePairs(params.reads, flat: true)
    .map { it -> tuple(it[0], it[1], it[2]) }

  // Raw QC
  raw_qc = FASTQC_RAW(read_pairs)

  // Trim
  trimmed = FASTP(read_pairs)

  // Split FASTP output:
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

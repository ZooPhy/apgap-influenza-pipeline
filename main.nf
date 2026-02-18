nextflow.enable.dsl=2

/*
  APGAP influenza starter QC pipeline
  Raw reads -> FastQC -> fastp trimming -> FastQC -> MultiQC
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
  publishDir "${params.outdir}/trim", mode: 'copy'
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
    tuple val(sample_id), path(r1), path(r2)

  output:
    path "*_fastqc.{zip,html}"

  container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

  """
  fastqc --threads ${task.cpus} ${r1} ${r2}
  """
}

process MULTIQC {
  publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

  input:
    path qc_files

  output:
    path "multiqc_report.html"
    path "multiqc_data"

  container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"

  """
  multiqc ${qc_files} -o .
  """
}

workflow {
  // fromFilePairs emits items that behave like: [ sample_id, r1, r2 ]
  read_pairs = Channel
    .fromFilePairs(params.reads, flat: true)
    .map { it -> tuple(it[0], it[1], it[2]) }

  // QC on raw reads
  raw_fastqc = FASTQC_RAW(read_pairs)

  // Trim + fastp reports
  trimmed = FASTP(read_pairs)

  // QC on trimmed reads (trimmed tuple includes html/json too)
  trimmed_fastqc = FASTQC_TRIMMED(
    trimmed.map { sid, r1t, r2t, html, json -> tuple(sid, r1t, r2t) }
  )

  // Feed FastQC outputs + fastp reports into MultiQC
  qc_files = raw_fastqc
    .mix(trimmed_fastqc)
    .mix(trimmed.map { sid, r1t, r2t, html, json -> html })
    .mix(trimmed.map { sid, r1t, r2t, html, json -> json })

  MULTIQC( qc_files.collect() )
}


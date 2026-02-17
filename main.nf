nextflow.enable.dsl=2

process FASTQC {
  tag "$sample_id"
  publishDir "${params.outdir}/qc/fastqc", mode: 'copy'
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
    path fastqc_reports

  output:
    path "multiqc_report.html"
    path "multiqc_data"

  container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"

  """
  multiqc ${fastqc_reports} -o .
  """
}

workflow {
  read_pairs = Channel.fromFilePairs(params.reads, flat: true)
  fastqc_out = FASTQC(read_pairs)
  MULTIQC( fastqc_out.collect() )
}

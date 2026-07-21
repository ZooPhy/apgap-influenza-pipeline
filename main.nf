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

  container params.cutadapt_container

  script:
    def g_opts = (params.primers_5p_r1 ?: []).collect { p -> "-g ${p}" }.join(' ')
    def G_opts = (params.primers_5p_r2 ?: []).collect { p -> "-G ${p}" }.join(' ')
    def overlap = params.cutadapt_min_overlap ?: 10
    def errrate = params.cutadapt_error_rate ?: 0.1

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

  container params.fastqc_container

  script:
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

  container params.fastp_container

  script:
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

  container params.fastqc_container

  script:
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

  container params.multiqc_container

  script:
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
    tuple val(sample_id), path("${sample_id}")

  container params.irma_container

  script:
  """
  IRMA ${params.irma_module} ${r1t} ${r2t} ${sample_id}
  """
}

process IRMA_COVERAGE {
  tag "$sample_id"
  publishDir "${params.outdir}/coverage/${sample_id}", mode: 'copy'

  input:
    tuple val(sample_id), path(sample_dir)
    path coverage_script

  output:
    tuple val(sample_id),
          path("${sample_id}_irma_coverage.tsv"),
          path("${sample_id}_pass_segments.txt"),
          path("${sample_id}.pass_segments.fa")

  container params.r_base_container

  script:
  """
  Rscript ${coverage_script} \
    ${sample_dir} \
    ${sample_id} \
    ${params.coverage_threshold}
  """
}

process BLASTN_IRMA_PASS_SEGMENTS {
  tag "$sample_id"
  publishDir "${params.outdir}/blast/${sample_id}", mode: 'copy'
  cpus params.threads

  container params.blast_container

  input:
    tuple val(sample_id),
          path(coverage_tsv),
          path(pass_segments_txt),
          path(pass_fasta)

  output:
    tuple val(sample_id),
          path("${sample_id}.irma_pass_segments.blast.tsv")

  script:
  """
  blastn \
    -query ${pass_fasta} \
    -db ${params.blast_db} \
    -out ${sample_id}.irma_pass_segments.blast.tsv \
    -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qlen slen evalue bitscore stitle" \
    -num_threads ${task.cpus}
  """
}

process VADR_IRMA_PASS_SEGMENTS {
  tag "$sample_id"
  publishDir "${params.outdir}/vadr/${sample_id}", mode: 'copy'
  cpus params.threads

  container params.vadr_container

  containerOptions = "--platform linux/amd64 -v ${projectDir}:${projectDir}"

  input:
    tuple val(sample_id),
          path(coverage_tsv),
          path(pass_segments_txt),
          path(pass_fasta)

  output:
    tuple val(sample_id),
          path("${sample_id}.vadr.*"),
          emit: vadr_results

  script:
  """
  v-annotate.pl \
    -f \
    --split \
    --cpu ${task.cpus} \
    -r \
    --atgonly \
    --xnocomp \
    --nomisc \
    --alt_fail extrant5,extrant3 \
    --mkey ${params.vadr_mkey} \
    ${pass_fasta} \
    ${sample_id}

  mv ${sample_id}/${sample_id}.vadr.* .
  """
}

process MAP_IRMA_TO_VADR_AA {
  tag "$sample_id"
  publishDir "${params.outdir}/variant_annotation/${sample_id}", mode: 'copy'

  input:
    tuple val(sample_id), path(vadr_files)
    path variant_annotation_script
    path genetic_code

  output:
    tuple val(sample_id), path("${sample_id}.irma_vadr_aa.tsv")

  container params.r_base_container

  script:
  """
  Rscript ${variant_annotation_script} \
    --sample_id ${sample_id} \
    --irma_dir ${projectDir}/${params.outdir}/assembly/irma/${sample_id} \
    --vadr_dir . \
    --genetic_code ${params.genetic_code} \
    --out ${sample_id}.irma_vadr_aa.tsv
  """
}

process IVAR_USING_IRMA_CONSENSUS {
  tag "$sample_id"
  publishDir "${params.outdir}/variants/ivar_irma_consensus/${sample_id}", mode: 'copy'
  cpus params.threads

  input:
    tuple val(sample_id), path(sample_dir), path(r1), path(r2)

  output:
    tuple val(sample_id),
          path("${sample_id}.irma.clean.fa"),
          path("${sample_id}.ivar.tsv"),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai")

  container params.ivar_container

  script:
  """
  cat ${sample_dir}/amended_consensus/${sample_id}_*.fa > ${sample_id}.irma.consensus.fa

  sed -E '
    s/^>.*_1.*/>PB2/
    s/^>.*_2.*/>PB1/
    s/^>.*_3.*/>PA/
    s/^>.*_4.*/>HA/
    s/^>.*_5.*/>NP/
    s/^>.*_6.*/>NA/
    s/^>.*_7.*/>MP/
    s/^>.*_8.*/>NS/
  ' ${sample_id}.irma.consensus.fa > ${sample_id}.irma.clean.fa

  bwa index ${sample_id}.irma.clean.fa
  samtools faidx ${sample_id}.irma.clean.fa

  bwa mem -t ${task.cpus} ${sample_id}.irma.clean.fa ${r1} ${r2} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

  samtools index ${sample_id}.sorted.bam

  samtools mpileup -A -d 0 -Q 0 -f ${sample_id}.irma.clean.fa ${sample_id}.sorted.bam | \
    ivar variants \
      -p ${sample_id}.ivar \
      -q ${params.ivar_min_qual} \
      -t ${params.ivar_min_freq}
  """
}

process IVAR_CONSENSUS_THEN_VARIANTS {
  tag "$sample_id"
  publishDir "${params.outdir}/variants/ivar_ref_consensus/${sample_id}", mode: 'copy'
  cpus params.threads

  container params.ivar_container

  input:
    tuple val(sample_id), path(r1), path(r2)
    path ref0

  output:
    tuple val(sample_id),
          path("${sample_id}.ivar.consensus.segmented.fa"),
          path("${sample_id}.ivar.variants.tsv"),
          path("${sample_id}.ivar_consensus.sorted.bam"),
          path("${sample_id}.ivar_consensus.sorted.bam.bai")

  script:
  """
  bwa index ${ref0}
  samtools faidx ${ref0}

  bwa mem -t ${task.cpus} ${ref0} ${r1} ${r2} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.ref.sorted.bam

  samtools index ${sample_id}.ref.sorted.bam

  mkdir -p segment_consensus

  for SEG in PB2 PB1 PA HA NP NA MP NS
  do
    samtools mpileup -A -d 0 -Q 0 -f ${ref0} -r \$SEG ${sample_id}.ref.sorted.bam | \
      ivar consensus \
        -p segment_consensus/${sample_id}.\${SEG} \
        -q ${params.ivar_min_qual} \
        -t ${params.ivar_consensus_freq} \
        -n N

    sed -E "s/^>.*/>\${SEG}/" \
      segment_consensus/${sample_id}.\${SEG}.fa \
      > segment_consensus/${sample_id}.\${SEG}.clean.fa
  done

  cat segment_consensus/${sample_id}.PB2.clean.fa \
      segment_consensus/${sample_id}.PB1.clean.fa \
      segment_consensus/${sample_id}.PA.clean.fa \
      segment_consensus/${sample_id}.HA.clean.fa \
      segment_consensus/${sample_id}.NP.clean.fa \
      segment_consensus/${sample_id}.NA.clean.fa \
      segment_consensus/${sample_id}.MP.clean.fa \
      segment_consensus/${sample_id}.NS.clean.fa \
      > ${sample_id}.ivar.consensus.segmented.fa

  bwa index ${sample_id}.ivar.consensus.segmented.fa
  samtools faidx ${sample_id}.ivar.consensus.segmented.fa

  bwa mem -t ${task.cpus} ${sample_id}.ivar.consensus.segmented.fa ${r1} ${r2} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.ivar_consensus.sorted.bam

  samtools index ${sample_id}.ivar_consensus.sorted.bam

  samtools mpileup -A -d 0 -Q 0 -f ${sample_id}.ivar.consensus.segmented.fa ${sample_id}.ivar_consensus.sorted.bam | \
    ivar variants \
      -p ${sample_id}.ivar.variants \
      -q ${params.ivar_min_qual} \
      -t ${params.ivar_min_freq}
  """
}

process IVAR_DEPTH {
  tag "$sample_id"
  publishDir "${params.outdir}/coverage/ivar/depth", mode: 'copy'
  cpus params.threads

  container params.ivar_container

  input:
    tuple val(sample_id), path(bam), path(bai)

  output:
    tuple val(sample_id), path("${sample_id}.ivar.depth.tsv")

  script:
  """
  samtools depth -a ${bam} > ${sample_id}.ivar.depth.tsv
  """
}

process IVAR_COVERAGE_SUMMARY {
  tag "$sample_id"
  publishDir "${params.outdir}/coverage/ivar", mode: 'copy'

  input:
    tuple val(sample_id), path(depth_tsv)
    path coverage_script

  output:
    tuple val(sample_id), path("${sample_id}_ivar_coverage.tsv")

  container params.r_base_container

  script:
  """
  Rscript ${coverage_script} \
    ${depth_tsv} \
    ${sample_id} \
    ${params.coverage_threshold}
  """
}

workflow {
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

  reads_for_qc = read_pairs
  if (params.trim_primers) {
    reads_for_qc = CUTADAPT_PRIMERS(read_pairs)
  }

  raw_qc = FASTQC_RAW(reads_for_qc)
  trimmed = FASTP(reads_for_qc)

  trimmed_reads = trimmed.map { sid, r1t, r2t, html, json -> tuple(sid, r1t, r2t) }
  fastp_reports = trimmed.map { sid, r1t, r2t, html, json -> [html, json] }.flatten()

  trimmed_qc = FASTQC_TRIMMED(trimmed_reads)

  qc_files = raw_qc
    .mix(trimmed_qc)
    .mix(fastp_reports)
    .collect()

  mq = MULTIQC(qc_files)
  done = mq.report.map { 1 }

  trimmed_after_multiqc = trimmed_reads
    .combine(done)
    .map { it -> tuple(it[0], it[1], it[2]) }

  irma_results = IRMA_RUN(trimmed_after_multiqc)

  irma_coverage = IRMA_COVERAGE(
    irma_results,
    file(params.coverage_script)
  )

  blast_results = BLASTN_IRMA_PASS_SEGMENTS(irma_coverage)
  vadr_results = VADR_IRMA_PASS_SEGMENTS(irma_coverage)

  variant_annotation_results = MAP_IRMA_TO_VADR_AA(
    vadr_results,
    file(params.variant_annotation_script),
    file(params.genetic_code)
  )

  /*
   * iVar sections disabled
   *
   * ivar_results = IVAR_CONSENSUS_THEN_VARIANTS(
   *   trimmed_after_multiqc,
   *   file(params.ivar_ref)
   * )
   *
   * ivar_bams = ivar_results.map { sid, cons_fa, variants_tsv, bam, bai ->
   *   tuple(sid, bam, bai)
   * }
   *
   * ivar_depth = IVAR_DEPTH(ivar_bams)
   *
   * IVAR_COVERAGE_SUMMARY(
   *   ivar_depth,
   *   file(params.ivar_coverage_script)
   * )
   */
}

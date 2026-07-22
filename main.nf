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
    def g_opts = (params.primers_5p_r1 ?: [])
      .collect { p -> "-g ${p}" }
      .join(' ')

    def G_opts = (params.primers_5p_r2 ?: [])
      .collect { p -> "-G ${p}" }
      .join(' ')

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
  	-i ${r1} \
  	-I ${r2} \
  	-o ${sample_id}_R1.trim.fastq.gz \
  	-O ${sample_id}_R2.trim.fastq.gz \
  	-h ${sample_id}.fastp.html \
  	-j ${sample_id}.fastp.json \
  	--detect_adapter_for_pe \
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
  cp -L ${coverage_script} irma_coverage_local.R

  Rscript irma_coverage_local.R \
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
    tuple val(sample_id), path(vadr_files), path(sample_dir)
    path variant_annotation_script
    path genetic_code

  output:
    tuple val(sample_id), path("${sample_id}.irma_vadr_aa.tsv")

  container params.r_base_container

  script:
  """
  cp -L ${variant_annotation_script} map_irma_to_vadr_aa_local.R
  cp -L ${genetic_code} genetic_code_local.csv

  Rscript map_irma_to_vadr_aa_local.R \
    --sample_id ${sample_id} \
    --irma_dir ${sample_dir} \
    --vadr_dir . \
    --genetic_code genetic_code_local.csv \
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


process SUBTYPE_EVIDENCE {
  tag "$sample_id"
  publishDir "${params.outdir}/subtype/${sample_id}", mode: 'copy'
  cpus 1
  memory '2 GB'
  stageInMode 'copy'

  input:
    tuple val(sample_id), path(sample_dir)

  output:
    tuple val(sample_id), path("${sample_id}.subtype_evidence.tsv")

  container params.samtools_container

  script:
  """
  printf 'sample_id\tsegment\tsubtype\tbam\treference_length\tmapped_reads\tmean_depth\tmedian_depth\tbreadth_at_threshold\n' \
    > ${sample_id}.subtype_evidence.tsv

  {
    find ${sample_dir} -type f -name 'A_HA_H*.bam'
    find ${sample_dir} -type f -name 'A_NA_N*.bam'
  } | sort -u | while IFS= read -r bam; do

    name=\$(basename "\$bam" .bam)
    segment=\$(printf '%s' "\$name" | cut -d_ -f2)
    subtype=\$(printf '%s' "\$name" | cut -d_ -f3)

    depth_file="\${name}.depth.tsv"
    sorted_file="\${name}.depths.sorted.txt"

    samtools depth -aa "\$bam" > "\$depth_file"
    cut -f3 "\$depth_file" | sort -n > "\$sorted_file"

    reference_length=\$(wc -l < "\$sorted_file" | tr -d ' ')
    mapped_reads=\$(samtools view -c -F 4 "\$bam")

    mean_depth=\$(awk '{s+=\$1} END {if (NR) printf "%.6f", s/NR; else print 0}' "\$sorted_file")
    median_depth=\$(awk '{a[NR]=\$1} END {if (NR==0) print 0; else if (NR%2) print a[(NR+1)/2]; else printf "%.6f", (a[NR/2]+a[NR/2+1])/2}' "\$sorted_file")
    breadth=\$(awk -v t=${params.subtype_min_median_depth} '{if (\$1>=t) n++} END {if (NR) printf "%.6f", n/NR; else print 0}' "\$sorted_file")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      '${sample_id}' "\$segment" "\$subtype" "\$name.bam" "\$reference_length" \
      "\$mapped_reads" "\$mean_depth" "\$median_depth" "\$breadth" \
      >> ${sample_id}.subtype_evidence.tsv
  done
  """
}

process SUMMARY_REPORT {
  tag "$sample_id"
  publishDir "${params.outdir}/summary", mode: 'copy'
  cpus 1
  memory '1 GB'
  stageInMode 'copy'

  input:
    tuple val(sample_id),
          path(coverage_tsv),
          path(pass_segments_txt),
          path(pass_fasta),
          path(blast_tsv),
          path(aa_tsv),
          path(subtype_tsv)
    path sample_summary_script

  output:
    tuple val(sample_id),
          path("${sample_id}.sample_summary.tsv"),
          path("${sample_id}.sample_summary.md")

  container params.r_base_container

  script:
  """
  cp -L ${coverage_tsv} coverage.tsv
  cp -L ${pass_segments_txt} pass_segments.txt
  cp -L ${blast_tsv} blast.tsv
  cp -L ${aa_tsv} aa.tsv
  cp -L ${subtype_tsv} subtype.tsv
  cp -L ${sample_summary_script} sample_summary_local.R

  Rscript sample_summary_local.R \
    --sample_id ${sample_id} \
    --coverage_tsv coverage.tsv \
    --pass_segments_txt pass_segments.txt \
    --blast_tsv blast.tsv \
    --aa_tsv aa.tsv \
    --host ${params.host} \
    --subtype_tsv subtype.tsv \
    --subtype_min_median_depth ${params.subtype_min_median_depth} \
    --subtype_min_breadth ${params.subtype_min_breadth} \
    --subtype_minor_fraction ${params.subtype_minor_fraction} \
    --out_tsv ${sample_id}.sample_summary.tsv \
    --out_md ${sample_id}.sample_summary.md
  """
}

workflow {
  def valid_hosts = ['human', 'bird', 'swine', 'environmental', 'other']

params.host = params.host
  .toString()
  .trim()
  .toLowerCase()

if (!valid_hosts.contains(params.host)) {
  error """
Invalid --host value: ${params.host}

Valid options:
  human
  bird
  swine
  environmental
  other
"""
}
  def build_read_pairs = {
    if (params.reads_dir) {
      def d = params.reads_dir.toString().replaceAll(/\/$/, '')

      def candidates = [
  		"${d}/*-R{1,2}.fastq.gz",
  		"${d}/*-R{1,2}_001.fastq.gz",
  		"${d}/*_{R1,R2}.fastq.gz",
  		"${d}/*_{R1,R2}_001.fastq.gz",
  		"${d}/*_{R1,R2}.fq.gz",
  		"${d}/*_{R1,R2}_001.fq.gz"
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

  // Raw FastQC always evaluates the original input reads.
  raw_qc = FASTQC_RAW(read_pairs)

  // Primer trimming is optional and occurs before fastp.
  reads_for_fastp = read_pairs
  if (params.trim_primers) {
    reads_for_fastp = CUTADAPT_PRIMERS(read_pairs)
  }

  trimmed = FASTP(reads_for_fastp)

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

  vadr_for_variant = vadr_results.map { sid, vadr_files ->
    tuple(sid, vadr_files)
  }

  irma_for_variant = irma_results.map { sid, sample_dir ->
    tuple(sid, sample_dir)
  }

  variant_annotation_inputs = vadr_for_variant
    .combine(irma_for_variant, by: 0)
    .map { sid, vadr_files, sample_dir ->
      tuple(sid, vadr_files, sample_dir)
    }

  variant_annotation_results = MAP_IRMA_TO_VADR_AA(
    variant_annotation_inputs,
    file(params.variant_annotation_script),
    file(params.genetic_code)
  )

  subtype_evidence = SUBTYPE_EVIDENCE(irma_results)

  coverage_for_summary = irma_coverage.map { sid, cov_tsv, pass_txt, pass_fa ->
    tuple(sid, cov_tsv, pass_txt, pass_fa)
  }

  blast_for_summary = blast_results.map { sid, blast_tsv ->
    tuple(sid, blast_tsv)
  }

  aa_for_summary = variant_annotation_results.map { sid, aa_tsv ->
    tuple(sid, aa_tsv)
  }

  subtype_for_summary = subtype_evidence.map { sid, subtype_tsv ->
    tuple(sid, subtype_tsv)
  }

  summary_inputs = coverage_for_summary
    .combine(blast_for_summary, by: 0)
    .combine(aa_for_summary, by: 0)
    .combine(subtype_for_summary, by: 0)
    .map { sid, cov_tsv, pass_txt, pass_fa, blast_tsv, aa_tsv, subtype_tsv ->
      tuple(sid, cov_tsv, pass_txt, pass_fa, blast_tsv, aa_tsv, subtype_tsv)
    }

  summary_report = SUMMARY_REPORT(
    summary_inputs,
    file("${projectDir}/scripts/sample_summary.R")
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

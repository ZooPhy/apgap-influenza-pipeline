# APGAP Influenza Pipeline

APGAP Influenza Pipeline is a Nextflow DSL2 workflow for influenza whole-genome analysis from paired-end Illumina sequencing data. The pipeline integrates **IRMA**, **BLAST**, and **VADR** to generate consensus genomes, identify influenza segments, annotate coding sequences, characterize amino acid variants, produce integrated per-sample summary reports, and generate a comprehensive run-level quality control report that summarizes sequencing performance across all samples.
The pipeline performs:

- Read quality control
- Automatic Illumina/Nextera adapter trimming with fastp
- Optional assay-specific 5' primer trimming with Cutadapt
- Genome assembly
- Coverage assessment
- Influenza segment identification
- Genome annotation
- Amino acid variant annotation
- Coverage-aware influenza A subtype determination
- Host-aware subtype interpretation
- Integrated summary reporting
- Run-level quality control reporting


---

# Workflow

```text
FASTQ Files
    │
    ▼
Raw Read Quality Control
(FastQC)
    │
    ▼
(Optional) Assay-Specific 5' Primer Trimming
(Cutadapt)
    │
    ▼
Adapter and Quality Trimming
(fastp)
    │
    ▼
Post-Trimming Quality Control
(FastQC → MultiQC)
    │
    ▼
IRMA Genome Assembly
    │
    ▼
Coverage Assessment
    ├── Coverage Report
    ├── PASS Segment List
    └── PASS Segment FASTA
    │
    ▼
BLAST
(Segment Identification)
    │
    ▼
VADR Annotation
(CDS and Gene Annotation)
    │
    ▼
Variant Annotation
(Map IRMA variants to VADR annotations)
    │
    ▼
Subtype Evidence
(IRMA HA/NA coverage and breadth)
    │
    ▼
Per-sample Summary Reports
    │
    ▼
Run Summary Report
(Multi-sample quality control and sequencing overview)
    │
    ▼
Final Results
    ├── Consensus Sequences
    ├── Coverage Statistics
    ├── Segment Identification
    ├── Gene/CDS Features
    ├── Amino Acid Variants
    ├── Mutation Effects
    ├── Subtype Evidence
    ├── Host-Subtype Review Flags
    └── Summary Reports
```

---

# Requirements

## Required Software

- Nextflow (DSL2)
- Java 17 or later

### Docker

- Docker Desktop or Docker Engine

### Apptainer

- Apptainer (or Singularity)

## BLAST Database

A local installation of BLAST+ is **not** required.

The `build_blast_db.sh` script automatically detects Docker, Apptainer, or Singularity and uses a containerized version of BLAST+ to build the influenza reference database. The same container runtime used by the pipeline is also used for database creation, ensuring consistent BLAST versions across platforms.

---

# First-Time Setup

Clone the repository:

```bash
git clone https://github.com/ZooPhy/apgap-influenza-pipeline.git
cd apgap-influenza-pipeline
```

Build the BLAST database:

```bash
./scripts/build_blast_db.sh
```

The setup script automatically:

1. Downloads `fluA_reference.fasta.zip` from the project's GitHub Releases page if it is not already present.
2. Extracts the influenza reference FASTA.
3. Builds the BLAST nucleotide database using a containerized version of `makeblastdb`.

The generated database will be stored in:

```text
FLU_DB/
└── fluA_db/
    ├── fluA_db.nhr
    ├── fluA_db.nin
    ├── fluA_db.nsq
    └── ...
```

The downloaded archive is retained in `FLU_DB/`, so future runs do not re-download it unless it is removed.

The BLAST database only needs to be built once and should **not** be committed to Git.

---

# Input

The pipeline accepts paired-end Illumina FASTQ files.

Example:

```text
data/
├── Sample1-R1.fastq.gz
├── Sample1-R2.fastq.gz
├── Sample2-R1.fastq.gz
└── Sample2-R2.fastq.gz
```

Supported naming conventions include both:

```text
*_R1.fastq.gz
*_R2.fastq.gz
```

and

```text
*-R1.fastq.gz
*-R2.fastq.gz
```

---

# Read Trimming

## Illumina and Nextera adapter trimming

Standard Illumina and Nextera adapter contamination is handled automatically by **fastp** during paired-end read processing. The pipeline does not require users to provide Illumina or Nextera adapter sequences for routine runs.

The fastp step performs adapter trimming and quality filtering before IRMA assembly. Adapter-trimmed reads and fastp reports are written under:

```text
results/trim/fastp/
```

## Optional assay-specific primer trimming

Cutadapt is reserved for known assay-specific primers, such as primers introduced during targeted amplicon generation. It is disabled by default.

Default configuration:

```groovy
trim_primers = false
primers_5p_r1 = []
primers_5p_r2 = []
```

Enable Cutadapt only when the exact 5' primer sequences used in the laboratory protocol are known:

```bash
nextflow run main.nf \
    -profile docker \
    --reads_dir data \
    --outdir results \
    --trim_primers true
```

The primer sequences should be supplied in `nextflow.config`, for example:

```groovy
primers_5p_r1 = ["ACTUAL_R1_PRIMER_SEQUENCE"]
primers_5p_r2 = ["ACTUAL_R2_PRIMER_SEQUENCE"]
```

Do not use generic influenza terminal primer sequences unless those exact primers were used during sample preparation. Do not use this Cutadapt step for ordinary Illumina or Nextera adapter removal; fastp already handles that stage.

Primer-trimmed reads are written under:

```text
results/trim/primers/
```

---

# Running the Pipeline

## Docker

```bash
nextflow run main.nf \
    -profile docker \
    --reads_dir data \
    --outdir results \
    --host human
```

## Apptainer

```bash
nextflow run main.nf \
    -profile apptainer \
    --reads_dir data \
    --outdir results \
    --host human
```

Resume a previous run:

```bash
nextflow run main.nf -resume
```

After successful completion:

- Review `results/run_summary/run_summary.md` for an overview of the entire sequencing run.
- Review `results/summary/` for detailed per-sample reports.

Cutadapt will only appear in the process list when `--trim_primers true` is enabled.
---

# Running on a SLURM Cluster

Example submission script:

```bash
#!/bin/bash
#SBATCH -p htc
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH -o jobout/pipeline_%j.out
#SBATCH -e jobout/pipeline_%j.err
#SBATCH --mail-type=ALL

module load openjdk
module load nextflow

export APPTAINER_CACHEDIR=$HOME/.apptainer/cache
mkdir -p "$APPTAINER_CACHEDIR"

export NXF_OPTS='-Xms1g -Xmx8g'

nextflow run main.nf \
    -profile apptainer \
    --reads_dir data \
    --outdir results \
    --host human \
    --threads ${SLURM_CPUS_PER_TASK} \
    -resume
```

Submit with:

```bash
sbatch run_pipeline.sh
```

---

# Output

```text
results/
├── qc/
├── trim/
├── assembly/
├── coverage/
├── blast/
├── vadr/
├── variant_annotation/
├── subtype/
├── summary/
└── run_summary/
```

## Output Directories

| Directory | Description |
|-----------|-------------|
| `qc/` | FastQC and MultiQC reports |
| `trim/` | fastp-trimmed reads and optional Cutadapt primer-trimmed reads |
| `assembly/` | IRMA assemblies and consensus genomes |
| `coverage/` | Coverage statistics and PASS segment FASTA files |
| `blast/` | Influenza segment identification |
| `vadr/` | Gene and CDS annotations |
| `variant_annotation/` | Amino acid variant annotations |
| `subtype/` | Per-sample HA and NA subtype evidence derived from IRMA subtype-specific BAM files |
| `summary/` | Integrated per-sample summary reports, including host-aware subtype interpretation |
| `run_summary/` | Run-level quality control report, figures, and aggregated summary tables |

---

# Coverage Assessment

Coverage is calculated independently for each influenza genome segment.

Segments with a median coverage greater than or equal to the specified threshold (default **20×**) are classified as **PASS**.

Generated files include:

- `sample_coverage.tsv`
- `sample_pass_segments.txt`
- `sample.pass_segments.fa`

Only PASS segments continue to BLAST, VADR, and amino acid annotation.


---

# Influenza A Subtype Determination

The pipeline determines a potential influenza A subtype from IRMA subtype-specific HA and NA BAM files, for example:

```text
A_HA_H3.bam
A_NA_N1.bam
```

For each subtype-specific BAM file, the pipeline calculates:

- Reference length
- Number of mapped reads
- Mean depth
- Median depth
- Breadth of coverage at the configured depth threshold

A subtype candidate is considered supported when it meets both configured thresholds:

```text
median depth >= subtype_min_median_depth
breadth >= subtype_min_breadth
```

Default settings are:

```groovy
subtype_min_median_depth = 20
subtype_min_breadth = 0.80
subtype_minor_fraction = 0.10
```

When one HA subtype and one NA subtype are supported, they are combined into a potential subtype, such as `H3N2` or `H3N1`.

If only HA or only NA is supported, the report identifies a partial subtype, for example:

```text
H3 (NA undetermined)
```

If more than one HA or NA subtype meets the absolute thresholds, the pipeline compares each secondary subtype with the dominant subtype. A secondary subtype is treated as meaningful when:

```text
secondary median depth / dominant median depth >= subtype_minor_fraction
```

When multiple meaningful HA or NA subtype signals remain, the sample is flagged as a potential mixed infection. This may reflect true coinfection, cross-sample contamination, barcode hopping, nonspecific assignment, or another sequencing artifact and requires review.

Subtype evidence files are written to:

```text
results/subtype/<sample>/<sample>.subtype_evidence.tsv
```

Each row contains:

```text
sample_id
segment
subtype
bam
reference_length
mapped_reads
mean_depth
median_depth
breadth_at_threshold
```

---

# Host-Aware Subtype Interpretation

The pipeline accepts a host category using `--host`.

The default is:

```text
human
```

Valid values are:

```text
human
bird
swine
environmental
other
```

Example:

```bash
nextflow run main.nf \
    -profile docker \
    --reads_dir data \
    --outdir results \
    --host human
```

The host category does not replace or suppress the analytical subtype call. Instead, it adds interpretation and review flags to the sample summary.

For human samples, `H1N1` and `H3N2` are treated as expected seasonal subtype combinations. Other combinations, such as `H3N1`, retain their analytical subtype call but are flagged for review because they may represent an unusual infection, reassortment, cross-sample contamination, barcode hopping, or nonspecific assignment.

For swine samples, `H1N1`, `H1N2`, and `H3N2` are treated as common subtype combinations. Other combinations are flagged for review.

For bird, environmental, and other samples, no narrow subtype whitelist is applied because subtype diversity is broader or the host context is less specific.

The summary also compares the median depth of the called HA and NA segments. A strong HA/NA depth imbalance is reported as a review flag but is not considered proof of contamination or barcode hopping.

---

# Variant Annotation

The variant annotation module combines:

- IRMA nucleotide variants
- BLAST segment identification
- VADR gene/CDS annotations

to report:

- Segment
- Gene
- CDS
- Codon number
- Codon position
- Reference nucleotide
- Alternate nucleotide
- Reference amino acid
- Alternate amino acid
- Amino acid substitution
- Synonymous or nonsynonymous classification

---

# Summary Reports

One of the final outputs of the pipeline is an integrated summary that combines results from multiple modules into a single report answering:

> **What did we find in this sample?**

Each sample summary includes:

- Host category
- Potential or partial influenza A subtype
- Host-subtype review warning, when applicable
- HA and NA median-depth comparison
- Missing or insufficient subtype evidence
- Low-level or mixed-subtype signals
- Number of genome segments meeting the coverage threshold
- PASS segment names
- Top BLAST match for each segment
- Number of amino acid substitutions
- Number of synonymous substitutions
- Number of nonsynonymous substitutions
- Stop-gained and stop-lost mutations (if present)

For each sample, the pipeline generates:

```text
results/summary/
├── sample1.sample_summary.tsv
├── sample1.sample_summary.md
├── sample2.sample_summary.tsv
├── sample2.sample_summary.md
└── ...
```

The Markdown report is intended for quick interpretation without inspecting each module individually. Subtype calls are reported as potential analytical results and should be interpreted together with host context, segment coverage, controls, and other samples from the sequencing run.

---

# Run Summary Report

After all samples have been processed, the pipeline automatically generates a run-level summary by aggregating all per-sample summary tables.

Outputs include:

```text
results/run_summary/
├── run_summary.md
├── run_summary.tsv
├── figures/
│   ├── run_summary.png
│   └── run_summary.pdf
└── tables/
    └── samples_requiring_review.tsv
```
The run summary provides a high-level overview of sequencing performance and analytical results across the entire sequencing run.

It includes:

- Overall sequencing statistics
- Counts of complete, near-complete, partial, and failed genomes
- Samples requiring immediate review
- Influenza subtype distribution
- Distribution of PASS genome segments
- Summary of review categories
- HA versus NA median depth comparison
- Amino acid variant burden
- Host category summary
- Influenza subtype summary
- Per-sample summary table

The report is intended to provide a rapid quality-control overview of an entire sequencing run while preserving detailed per-sample reports for downstream investigation.  
---

# Repository Structure

```text
.
├── conf/
├── data/
├── FLU_DB/
│   ├── fluA_reference.fasta.zip   (downloaded automatically if needed)
│   └── fluA_db/                   (generated automatically)
├── resources/
│   └── genetic_code_translation_table_1.csv
├── scripts/
│   ├── build_blast_db.sh
│   ├── irma_coverage.R
│   ├── map_irma_to_vadr_aa.R
│   ├── sample_summary.R
│   └── run_summary.R
├── main.nf
├── nextflow.config
└── README.md
```

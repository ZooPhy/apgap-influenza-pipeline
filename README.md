# apgap-influenza-pipeline

## Run locally

### Flags

- `--with-docker`
  Run nextflow with Docker instance
- `--reads_dir "<folder>"`
  Path to a folder containing paired-end FASTQs. This is the easiest way to run the pipeline.
  The pipeline will automatically try common Illumina naming patterns (e.g. `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`).
- `-c <config>`
  Load an additional Nextflow config file (merged with `nextflow.config`). Useful for shared primer profiles.
- `--trim_primers "<default false>"`
  A boolean flag to run CUTADAPT with primers specified in the param block of the config file. For example,
   ```
  primers_5p_r1 = ["AGCAAAAGCAGG"]   // Uni-12
  primers_5p_r2 = ["AGTAGAAACAAGG"]  // Uni-13
   ```
  
### Example runs  

```
Examples:
  nextflow run main.nf -with-docker --reads_dir "./my/reads/" --trim_primers true -c ./conf/primers.conf
  ```
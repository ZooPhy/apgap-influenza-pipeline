#!/usr/bin/env Rscript

# ============================================================
# Map IRMA nucleotide variants to VADR CDS annotations
# Output: codon number, codon position, AA change, effect
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    stop(paste("Missing required argument:", flag))
  }
  args[i + 1]
}

sample_id <- get_arg("--sample_id")
irma_dir <- get_arg("--irma_dir")
vadr_dir <- get_arg("--vadr_dir")
genetic_code_file <- get_arg("--genetic_code")
out_file <- get_arg("--out")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Read genetic code
# -----------------------------

genetic_code <- read.csv(
  genetic_code_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("Codon", "AminoAcid")

if (!all(required_cols %in% names(genetic_code))) {
  stop(
    "Genetic code file must contain columns: ",
    paste(required_cols, collapse = ", ")
  )
}

codon_table <- setNames(
  genetic_code$AminoAcid,
  genetic_code$Codon
)

translate_codon <- function(codon) {
  codon <- toupper(codon)

  if (nchar(codon) != 3 || grepl("[^ACGT]", codon)) {
    return("X")
  }

  aa <- codon_table[[codon]]

  if (is.null(aa)) {
    return("X")
  }

  aa
}

# -----------------------------
# 2. Read VADR FASTA
# -----------------------------

vadr_fa <- file.path(vadr_dir, paste0(sample_id, ".vadr.pass.fa"))
vadr_ftr <- file.path(vadr_dir, paste0(sample_id, ".vadr.ftr"))

if (!file.exists(vadr_fa)) {
  stop("Missing VADR pass FASTA: ", vadr_fa)
}

if (!file.exists(vadr_ftr)) {
  stop("Missing VADR feature table: ", vadr_ftr)
}

fasta_lines <- readLines(vadr_fa)
headers <- grep("^>", fasta_lines)
seqs <- list()

for (i in seq_along(headers)) {
  start <- headers[i]
  end <- if (i < length(headers)) headers[i + 1] - 1 else length(fasta_lines)

  seq_name <- sub("^>", "", fasta_lines[start])
  seq_name <- strsplit(seq_name, "\\s+")[[1]][1]

  seqs[[seq_name]] <- toupper(
    paste(fasta_lines[(start + 1):end], collapse = "")
  )
}

# -----------------------------
# 3. Read VADR .ftr CDS features
# -----------------------------

ftr_lines <- readLines(vadr_ftr)
ftr_lines <- ftr_lines[!grepl("^#", ftr_lines)]
ftr_lines <- ftr_lines[nchar(trimws(ftr_lines)) > 0]

features <- list()

for (line in ftr_lines) {
  fields <- strsplit(trimws(line), "\\s+")[[1]]

  if (length(fields) < 10) next
  if (fields[6] != "CDS") next

  seq_name <- fields[2]
  segment <- sub("^.*_", "", seq_name)
  protein <- fields[7]

  coord_tokens <- fields[grepl("\\.\\.", fields) & grepl(":", fields)]
  if (length(coord_tokens) == 0) next

  # First coordinate token = coordinates on sample sequence
  coords_text <- coord_tokens[1]
  coord_parts <- strsplit(coords_text, ",")[[1]]

  coords <- data.frame(
    start = integer(),
    end = integer(),
    strand = character(),
    stringsAsFactors = FALSE
  )

  for (coord in coord_parts) {
    m <- regmatches(
      coord,
      regexec("([0-9]+)\\.\\.([0-9]+):([+-])", coord)
    )[[1]]

    if (length(m) == 4) {
      coords <- rbind(
        coords,
        data.frame(
          start = as.integer(m[2]),
          end = as.integer(m[3]),
          strand = m[4],
          stringsAsFactors = FALSE
        )
      )
    }
  }

  features[[length(features) + 1]] <- list(
    seq_name = seq_name,
    segment = segment,
    protein = protein,
    coords_text = coords_text,
    coords = coords
  )
}

if (length(features) == 0) {
  stop("No CDS features parsed from: ", vadr_ftr)
}

# -----------------------------
# 4. Read all IRMA variant files
# -----------------------------

variant_files <- list.files(
  file.path(irma_dir, "tables"),
  pattern = "variants.txt$",
  full.names = TRUE
)

if (length(variant_files) == 0) {
  stop("No IRMA variant files found in: ", file.path(irma_dir, "tables"))
}

variants <- data.frame()

for (vf in variant_files) {
  x <- read.table(
    vf,
    header = TRUE,
    stringsAsFactors = FALSE,
    fill = TRUE,
    check.names = FALSE
  )

  if (nrow(x) == 0) next

  for (i in seq_len(nrow(x))) {
    irma_segment <- x$Reference_Name[i]

    parts <- strsplit(irma_segment, "_")[[1]]
    segs <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "M", "NS")
    hit <- intersect(parts, segs)

    if (length(hit) == 0) {
      segment <- irma_segment
    } else if (hit[1] == "M") {
      segment <- "MP"
    } else {
      segment <- hit[1]
    }

    variants <- rbind(
      variants,
      data.frame(
        irma_file = basename(vf),
        irma_segment = irma_segment,
        segment = segment,
        nt_pos = as.integer(x$Position[i]),
        depth = x$Total[i],
        nt_ref = toupper(x$Consensus_Allele[i]),
        nt_alt = toupper(x$Minority_Allele[i]),
        freq = x$Minority_Frequency[i],
        stringsAsFactors = FALSE
      )
    )
  }
}

if (nrow(variants) == 0) {
  warning("No variants found. Writing empty output.")
  write.table(
    variants,
    out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  quit(save = "no", status = 0)
}

# -----------------------------
# 5. Map IRMA variants to VADR CDS features
# -----------------------------

results <- data.frame()

for (i in seq_len(nrow(variants))) {
  v <- variants[i, ]

  for (feat in features) {
    if (feat$segment != v$segment) next

    coords <- feat$coords

    cds_pos <- NA_integer_
    offset <- 0

    for (j in seq_len(nrow(coords))) {
      start <- coords$start[j]
      end <- coords$end[j]

      if (v$nt_pos >= start && v$nt_pos <= end) {
        cds_pos <- offset + (v$nt_pos - start + 1)
        break
      }

      offset <- offset + (end - start + 1)
    }

    if (is.na(cds_pos)) next

    segment_seq <- seqs[[feat$seq_name]]
    if (is.null(segment_seq)) next

    cds_pieces <- c()

    for (j in seq_len(nrow(coords))) {
      cds_pieces <- c(
        cds_pieces,
        substr(segment_seq, coords$start[j], coords$end[j])
      )
    }

    cds_seq <- paste(cds_pieces, collapse = "")

    if (nchar(cds_seq) %% 3 != 0) {
      warning(
        paste(
          sample_id,
          feat$protein,
          "CDS length is not divisible by 3."
        )
      )
    }

    codon_number <- floor((cds_pos - 1) / 3) + 1
    codon_position <- ((cds_pos - 1) %% 3) + 1

    codon_start <- (codon_number - 1) * 3 + 1
    codon_end <- codon_start + 2

    if (codon_end > nchar(cds_seq)) next

    consensus_codon <- substr(cds_seq, codon_start, codon_end)

    consensus_codon_vec <- strsplit(consensus_codon, "")[[1]]

    # Verify that the IRMA consensus nucleotide agrees with the
    # nucleotide in the VADR consensus sequence.
    if (consensus_codon_vec[codon_position] != v$nt_ref) {

      warning(
        paste(
          sample_id,
          feat$protein,
          "position", v$nt_pos,
          "- IRMA consensus:", v$nt_ref,
          "VADR consensus:", consensus_codon_vec[codon_position]
        )
      )

      next
    }

    variant_codon_vec <- consensus_codon_vec

    variant_codon_vec[codon_position] <- v$nt_alt

    consensus_codon_final <- paste(consensus_codon_vec, collapse = "")
    variant_codon_final <- paste(variant_codon_vec, collapse = "")

    consensus_aa <- translate_codon(consensus_codon_final)
    variant_aa <- translate_codon(variant_codon_final)

    aa_change <- paste0(consensus_aa, codon_number, variant_aa)

    effect <- ifelse(
      consensus_aa == variant_aa,
      "synonymous",
      ifelse(
        variant_aa == "*",
        "stop_gained",
        ifelse(consensus_aa == "*", "stop_lost", "nonsynonymous")
      )
    )

    results <- rbind(
      results,
      data.frame(
        sample_id = sample_id,
        irma_file = v$irma_file,
        irma_segment = v$irma_segment,
        nt_pos = v$nt_pos,
        nt_ref = v$nt_ref,
        nt_alt = v$nt_alt,
        depth = v$depth,
        freq = v$freq,
        vadr_seq = feat$seq_name,
        protein = feat$protein,
        cds_coords = feat$coords_text,
        cds_pos = cds_pos,
        codon_number = codon_number,
        codon_position = codon_position,
        consensus_codon = consensus_codon_final,
        variant_codon = variant_codon_final,
        consensus_aa = consensus_aa,
        variant_aa = variant_aa,
        aa_change = aa_change,
        effect = effect,
        stringsAsFactors = FALSE
      )
    )
  }
}

# -----------------------------
# 6. Write output
# -----------------------------

write.table(
  results,
  out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript irma_coverage.R <sample_dir> <sample_id> <median_threshold>")
}

sample_dir <- args[1]
sample_id <- args[2]
median_threshold <- as.numeric(args[3])

if (is.na(median_threshold)) {
  stop("median_threshold must be numeric")
}

segment_from_name <- function(x) {
  x <- basename(x)
  x <- sub("-coverage\\.txt$", "", x)
  x <- sub("^A_", "", x)
  x <- toupper(x)

  # Canonicalize HA subtype names such as H1, H3, H5.
  if (grepl("^HA", x) || grepl("^H[0-9]+", x)) {
    return("HA")
  }

  # Canonicalize NA subtype names such as N1, N2, N5.
  if (grepl("^NA", x) || grepl("^N[0-9]+", x)) {
    return("NA")
  }

  for (s in c("PB2", "PB1", "PA", "NP", "MP", "NS")) {
    if (startsWith(x, s)) {
      return(s)
    }
  }

  strsplit(x, "_", fixed = TRUE)[[1]][1]
}

get_cov_stats <- function(f, sample_name, threshold) {
  df <- tryCatch(
    read.delim(f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(df) || nrow(df) == 0) return(NULL)

  if (!"Coverage Depth" %in% colnames(df)) {
    stop(paste("Coverage Depth column not found in", f))
  }

  depths <- suppressWarnings(as.numeric(df[["Coverage Depth"]]))
  depths <- depths[!is.na(depths)]

  if (length(depths) == 0) return(NULL)

  contig <- if ("Reference_Name" %in% colnames(df)) {
    ref_names <- unique(df[["Reference_Name"]])
    ref_names <- ref_names[!is.na(ref_names) & ref_names != ""]
    if (length(ref_names) > 0) ref_names[1] else sub("-coverage\\.txt$", "", basename(f))
  } else {
    sub("-coverage\\.txt$", "", basename(f))
  }

  med <- median(depths)
  mean_cov <- mean(depths)
  breadth <- mean(depths > 0)

  data.frame(
    sample = sample_name,
    segment = segment_from_name(contig),
    contig = contig,
    length = length(depths),
    median_coverage = round(med, 2),
    mean_coverage = round(mean_cov, 2),
    breadth_covered = round(breadth, 3),
    coverage_flag = ifelse(med >= threshold, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )
}

tables_dir <- file.path(sample_dir, "tables")

if (!dir.exists(tables_dir)) {
  stop(paste("tables directory not found in", sample_dir))
}

cov_files <- list.files(
  tables_dir,
  pattern = "-coverage\\.txt$",
  full.names = TRUE
)

if (length(cov_files) == 0) {
  stop(paste("No coverage files found in", tables_dir))
}

sample_res <- lapply(
  cov_files,
  get_cov_stats,
  sample_name = sample_id,
  threshold = median_threshold
)

sample_res <- Filter(Negate(is.null), sample_res)

if (length(sample_res) == 0) {
  stop(paste("No usable coverage files found in", tables_dir))
}

res <- do.call(rbind, sample_res)
valid_segments <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
unknown_segments <- setdiff(unique(res$segment), valid_segments)

if (length(unknown_segments) > 0) {
  warning(
    "Unrecognized segment names: ",
    paste(unknown_segments, collapse = ", ")
  )
}

seg_order <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
res$segment_order <- match(res$segment, seg_order)
res$segment_order[is.na(res$segment_order)] <- 999
res <- res[order(res$segment_order, res$contig), ]
res$segment_order <- NULL

coverage_tsv <- paste0(sample_id, "_irma_coverage.tsv")

write.table(
  res,
  file = coverage_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote IRMA coverage summary to: ", coverage_tsv)

pass_segments <- unique(res$segment[res$coverage_flag == "PASS"])

pass_segment_file <- paste0(sample_id, "_pass_segments.txt")

writeLines(pass_segments, pass_segment_file)

message("Wrote PASS segment list to: ", pass_segment_file)

segment_to_file_number <- c(
  "PB2" = "1",
  "PB1" = "2",
  "PA"  = "3",
  "HA"  = "4",
  "NP"  = "5",
  "NA"  = "6",
  "MP"  = "7",
  "NS"  = "8"
)

amended_dir <- file.path(sample_dir, "amended_consensus")
pass_fasta <- paste0(sample_id, ".pass_segments.fa")

file.create(pass_fasta)

if (!dir.exists(amended_dir)) {
  warning(paste("amended_consensus directory not found:", amended_dir))
} else {
  for (seg in pass_segments) {
    file_num <- segment_to_file_number[[seg]]

    if (is.null(file_num)) {
      warning(paste("No FASTA file mapping found for segment:", seg))
      next
    }

    fasta_file <- file.path(amended_dir, paste0(sample_id, "_", file_num, ".fa"))

    if (!file.exists(fasta_file)) {
      warning(paste("Consensus FASTA not found:", fasta_file))
      next
    }

    if (file.info(fasta_file)$size == 0) {
      warning(paste("Consensus FASTA is empty:", fasta_file))
      next
    }

    seq_lines <- readLines(fasta_file, warn = FALSE)

    seq_lines[grepl("^>", seq_lines)] <- paste0(">", sample_id, "_", seg)

    cat(seq_lines, sep = "\n", file = pass_fasta, append = TRUE)
    cat("\n", file = pass_fasta, append = TRUE)
  }
}

message("Wrote PASS consensus FASTA to: ", pass_fasta)

if (length(pass_segments) == 0) {
  warning("No segments passed coverage threshold. PASS FASTA will be empty.")
} else {
  message("PASS segments: ", paste(pass_segments, collapse = ", "))
}
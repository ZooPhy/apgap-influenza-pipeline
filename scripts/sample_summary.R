#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    stop(paste("Missing required argument:", flag))
  }
  args[i + 1]
}

sample_id <- get_arg("--sample_id")
coverage_tsv <- get_arg("--coverage_tsv")
pass_segments_txt <- get_arg("--pass_segments_txt")
blast_tsv <- get_arg("--blast_tsv")
aa_tsv <- get_arg("--aa_tsv")
subtype_tsv <- get_arg("--subtype_tsv")
out_tsv <- get_arg("--out_tsv")
out_md <- get_arg("--out_md")
host <- tolower(trimws(get_arg("--host")))

valid_hosts <- c(
  "human",
  "bird",
  "swine",
  "environmental",
  "other"
)

if (!host %in% valid_hosts) {
  stop(
    "Invalid --host value: ",
    host,
    ". Valid values: ",
    paste(valid_hosts, collapse = ", ")
  )
}


min_median_depth <- as.numeric(get_arg("--subtype_min_median_depth"))
min_breadth <- as.numeric(get_arg("--subtype_min_breadth"))
minor_fraction <- as.numeric(get_arg("--subtype_minor_fraction"))

safe_read_tsv <- function(path, header = TRUE) {
  if (!file.exists(path) || file.size(path) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  tryCatch(
    read.delim(
      path,
      header = header,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = "",
      fill = TRUE,
      na.strings = ""
    ),
    error = function(e) {
      data.frame(stringsAsFactors = FALSE)
    }
  )
}


get_called_depth <- function(evidence, segment, subtype) {
  if (
    nrow(evidence) == 0 ||
    !nzchar(subtype) ||
    !all(c("segment", "subtype", "median_depth") %in% names(evidence))
  ) {
    return(NA_real_)
  }

  rows <- evidence[
    !is.na(evidence$segment) &
      !is.na(evidence$subtype) &
      toupper(trimws(as.character(evidence$segment))) == segment &
      toupper(trimws(as.character(evidence$subtype))) == subtype,
    ,
    drop = FALSE
  ]

  if (nrow(rows) == 0) {
    return(NA_real_)
  }

  suppressWarnings(as.numeric(rows$median_depth[1]))
}

evaluate_host_subtype <- function(
  host,
  potential_subtype,
  subtype_status,
  ha_call,
  na_call,
  subtype_evidence
) {
  result <- list(
    warning = FALSE,
    level = "none",
    category = "",
    message = ""
  )

  if (subtype_status %in% c("mixed", "partial", "undetermined")) {
    return(result)
  }

  subtype <- toupper(trimws(potential_subtype))

  if (host == "human") {
    expected_human_subtypes <- c(
      "H1N1",
      "H3N2"
    )

    if (!subtype %in% expected_human_subtypes) {
      result$warning <- TRUE
      result$level <- "review"
      result$category <- "host_subtype_discordance"
      result$message <- paste0(
        subtype,
        " is not a typical seasonal human influenza A subtype. ",
        "Review subtype-specific depth, sample metadata, barcode assignments, ",
        "negative controls, and other samples in the sequencing run. ",
        "The result could reflect an unusual infection, reassortment, ",
        "cross-sample contamination, barcode hopping, or nonspecific assignment."
      )
    }
  }

  if (host == "swine") {
    common_swine_subtypes <- c(
      "H1N1",
      "H1N2",
      "H3N2"
    )

    if (!subtype %in% common_swine_subtypes) {
      result$warning <- TRUE
      result$level <- "review"
      result$category <- "host_subtype_discordance"
      result$message <- paste0(
        subtype,
        " is outside the configured common swine subtype set. ",
        "Review the subtype evidence and sample metadata."
      )
    }
  }

  if (host == "bird") {
    # Avian influenza A viruses are highly diverse, so do not use a narrow
    # whitelist. Mixed or unusual subtype combinations should be reviewed
    # through the existing evidence logic.
    return(result)
  }

  if (host %in% c("environmental", "other")) {
    # No host-specific subtype whitelist is applied.
    return(result)
  }

  result
}


clean_subtypes <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x) & toupper(x) != "NA"])
}

read_blast_tsv <- function(path) {
  x <- safe_read_tsv(path, header = FALSE)

  cols <- c(
    "qseqid", "sseqid", "pident", "qcovs", "length",
    "mismatch", "gapopen", "qlen", "slen",
    "evalue", "bitscore", "stitle"
  )

  if (nrow(x) > 0 && ncol(x) >= length(cols)) {
    names(x)[seq_along(cols)] <- cols
  }

  x
}

empty_segment_result <- function() {
  list(
    call = NA_character_,
    status = "undetermined",
    candidates = character(),
    low_level = character()
  )
}

classify_segment <- function(x, segment) {
  required <- c("segment", "subtype", "median_depth", "breadth_at_threshold")

  if (nrow(x) == 0 || !all(required %in% names(x))) {
    return(empty_segment_result())
  }

  x <- x[
    !is.na(x$segment) &
      toupper(trimws(as.character(x$segment))) == toupper(segment) &
      !is.na(x$subtype) &
      nzchar(trimws(as.character(x$subtype))) &
      toupper(trimws(as.character(x$subtype))) != "NA",
    ,
    drop = FALSE
  ]

  if (nrow(x) == 0) {
    return(empty_segment_result())
  }

  x$subtype <- trimws(as.character(x$subtype))
  x$median_depth <- suppressWarnings(as.numeric(x$median_depth))
  x$breadth <- suppressWarnings(as.numeric(x$breadth_at_threshold))

  x <- x[
    !is.na(x$median_depth) & is.finite(x$median_depth) &
      !is.na(x$breadth) & is.finite(x$breadth),
    ,
    drop = FALSE
  ]

  if (nrow(x) == 0) {
    return(empty_segment_result())
  }

  x <- x[order(x$subtype, -x$median_depth, -x$breadth), , drop = FALSE]
  x <- x[!duplicated(x$subtype), , drop = FALSE]
  x <- x[order(-x$median_depth, -x$breadth, x$subtype), , drop = FALSE]

  qualified <- x[
    x$median_depth >= min_median_depth & x$breadth >= min_breadth,
    ,
    drop = FALSE
  ]

  if (nrow(qualified) == 0) {
    return(list(
      call = NA_character_,
      status = "undetermined",
      candidates = character(),
      low_level = clean_subtypes(x$subtype)
    ))
  }

  dominant_depth <- qualified$median_depth[1]
  if (is.na(dominant_depth) || !is.finite(dominant_depth) || dominant_depth <= 0) {
    return(empty_segment_result())
  }

  qualified$relative_depth <- qualified$median_depth / dominant_depth

  meaningful <- qualified[
    !is.na(qualified$relative_depth) &
      qualified$relative_depth >= minor_fraction,
    ,
    drop = FALSE
  ]

  low_level <- qualified[
    !is.na(qualified$relative_depth) &
      qualified$relative_depth < minor_fraction,
    ,
    drop = FALSE
  ]

  candidates <- clean_subtypes(meaningful$subtype)
  low <- clean_subtypes(low_level$subtype)

  if (length(candidates) == 0) {
    return(list(
      call = NA_character_,
      status = "undetermined",
      candidates = character(),
      low_level = low
    ))
  }

  if (length(candidates) == 1) {
    return(list(
      call = candidates[1],
      status = "single",
      candidates = candidates,
      low_level = low
    ))
  }

  list(
    call = paste(candidates, collapse = "/"),
    status = "mixed",
    candidates = candidates,
    low_level = low
  )
}

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)

coverage <- safe_read_tsv(coverage_tsv)
blast <- read_blast_tsv(blast_tsv)
aa <- safe_read_tsv(aa_tsv)
subtype_evidence <- safe_read_tsv(subtype_tsv)
pass_segments <- character()
if (file.exists(pass_segments_txt) && file.size(pass_segments_txt) > 0) {
  pass_segments <- unique(trimws(readLines(pass_segments_txt, warn = FALSE)))
  pass_segments <- pass_segments[
    !is.na(pass_segments) & nzchar(pass_segments) & toupper(pass_segments) != "NA"
  ]
}


blast_segments <- character()
top_hits <- data.frame(stringsAsFactors = FALSE)

if (nrow(blast) > 0 && "qseqid" %in% names(blast)) {
  blast$qseqid <- as.character(blast$qseqid)
  blast_segments <- unique(blast$qseqid[!is.na(blast$qseqid) & nzchar(blast$qseqid)])

  if ("bitscore" %in% names(blast)) {
    blast$bitscore_num <- suppressWarnings(as.numeric(blast$bitscore))
    blast <- blast[order(blast$qseqid, -blast$bitscore_num), , drop = FALSE]
  }

  top_hits <- blast[
    !is.na(blast$qseqid) & nzchar(blast$qseqid) & !duplicated(blast$qseqid),
    ,
    drop = FALSE
  ]
}

ha <- classify_segment(subtype_evidence, "HA")
na <- classify_segment(subtype_evidence, "NA")

ha$candidates <- clean_subtypes(ha$candidates)
na$candidates <- clean_subtypes(na$candidates)
ha$low_level <- clean_subtypes(ha$low_level)
na$low_level <- clean_subtypes(na$low_level)

if (length(ha$candidates) == 1 && length(na$candidates) == 1) {
  potential_subtype <- paste0(ha$candidates[1], na$candidates[1])
  subtype_status <- "single"
  ha_call <- ha$candidates[1]
  na_call <- na$candidates[1]
  missing_subtype_segments <- ""
} else if (length(ha$candidates) > 1 || length(na$candidates) > 1) {
  potential_subtype <- "Potential mixed infection"
  subtype_status <- "mixed"
  ha_call <- paste(ha$candidates, collapse = "/")
  na_call <- paste(na$candidates, collapse = "/")
  missing_subtype_segments <- ""
} else {
  ha_call <- if (length(ha$candidates) == 1) ha$candidates[1] else ""
  na_call <- if (length(na$candidates) == 1) na$candidates[1] else ""

  missing <- character()
  if (!nzchar(ha_call)) missing <- c(missing, "HA")
  if (!nzchar(na_call)) missing <- c(missing, "NA")
  missing_subtype_segments <- paste(missing, collapse = ",")

  if (nzchar(ha_call) && !nzchar(na_call)) {
    potential_subtype <- paste0("Partial subtype: ", ha_call, " (NA undetermined)")
    subtype_status <- "partial"
  } else if (!nzchar(ha_call) && nzchar(na_call)) {
    potential_subtype <- paste0("Partial subtype: ", na_call, " (HA undetermined)")
    subtype_status <- "partial"
  } else {
    potential_subtype <- "Undetermined"
    subtype_status <- "undetermined"
  }
}

effect_count <- function(effect_name) {
  if (nrow(aa) == 0 || !"effect" %in% names(aa)) return(0L)
  sum(aa$effect == effect_name, na.rm = TRUE)
}

coverage_segments <- if (nrow(coverage) > 0) nrow(coverage) else length(pass_segments)
pass_count <- length(pass_segments)
blast_hits <- nrow(blast)
aa_changes <- nrow(aa)
syn <- effect_count("synonymous")
nonsyn <- effect_count("nonsynonymous")
stop_gained <- effect_count("stop_gained")
stop_lost <- effect_count("stop_lost")

ha_called_depth <- get_called_depth(
  subtype_evidence,
  "HA",
  ha_call
)

na_called_depth <- get_called_depth(
  subtype_evidence,
  "NA",
  na_call
)

ha_na_depth_ratio <- if (
  !is.na(ha_called_depth) &&
  !is.na(na_called_depth) &&
  max(ha_called_depth, na_called_depth) > 0
) {
  min(ha_called_depth, na_called_depth) /
    max(ha_called_depth, na_called_depth)
} else {
  NA_real_
}

host_subtype_depth_imbalance <- (
  subtype_status == "single" &&
  !is.na(ha_na_depth_ratio) &&
  ha_na_depth_ratio < 0.01
)

host_assessment <- evaluate_host_subtype(
  host = host,
  potential_subtype = potential_subtype,
  subtype_status = subtype_status,
  ha_call = ha_call,
  na_call = na_call,
  subtype_evidence = subtype_evidence
)

summary_tsv <- data.frame(
  sample_id = sample_id,
  host = host,
  potential_subtype = potential_subtype,
  subtype_status = subtype_status,
  ha_call = ha_call,
  na_call = na_call,
  missing_subtype_segments = missing_subtype_segments,
  ha_low_level = paste(ha$low_level, collapse = ","),
  na_low_level = paste(na$low_level, collapse = ","),
  subtype_min_median_depth = min_median_depth,
  subtype_min_breadth = min_breadth,
  subtype_minor_fraction = minor_fraction,
  host_subtype_warning = host_assessment$warning,
  host_subtype_warning_level = host_assessment$level,
  host_subtype_warning_category = host_assessment$category,
  host_subtype_warning_message = host_assessment$message,
  ha_called_median_depth = ha_called_depth,
  na_called_median_depth = na_called_depth,
  ha_na_depth_ratio = ha_na_depth_ratio,
  host_subtype_depth_imbalance = host_subtype_depth_imbalance,
  coverage_segments = coverage_segments,
  pass_segments = pass_count,
  blast_hits = blast_hits,
  blast_segments = paste(blast_segments, collapse = ","),
  aa_changes = aa_changes,
  synonymous_changes = syn,
  nonsynonymous_changes = nonsyn,
  stop_gained = stop_gained,
  stop_lost = stop_lost,
  stringsAsFactors = FALSE
)

write.table(summary_tsv, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

subtype_md_line <- if (subtype_status == "partial") {
  if (nzchar(ha_call) && !nzchar(na_call)) {
    paste0("**Partial influenza A subtype: ", ha_call, " (NA undetermined)**")
  } else if (!nzchar(ha_call) && nzchar(na_call)) {
    paste0("**Partial influenza A subtype: ", na_call, " (HA undetermined)**")
  } else {
    "**Partial influenza A subtype: Undetermined**"
  }
} else {
  paste0("**Potential influenza A subtype: ", potential_subtype, "**")
}

md <- c(
  paste0("# Sample Summary: ", sample_id),
  "",
  paste0("**Host category: ", host, "**"),
  "",
  subtype_md_line,
  "",
  paste0(
    "Subtype evidence thresholds: median depth >= ", min_median_depth,
    "x, breadth >= ", round(min_breadth * 100),
    "%, secondary-to-dominant depth fraction >= ", minor_fraction, "."
  ),
  "",
  paste0(
    "This sample had **", pass_count, " PASS segments** out of **",
    coverage_segments, " covered segments**."
  )
)

if (isTRUE(host_assessment$warning)) {
  md <- c(
    md,
    "",
    "## Host-subtype review warning",
    paste0("- Host category: **", host, "**"),
    paste0("- Subtype call: **", potential_subtype, "**"),
    paste0("- ", host_assessment$message)
  )
}

if (
  subtype_status == "single" &&
  !is.na(ha_called_depth) &&
  !is.na(na_called_depth)
) {
  md <- c(
    md,
    "",
    "## HA/NA depth comparison",
    paste0("- ", ha_call, " median depth: **", ha_called_depth, "x**"),
    paste0("- ", na_call, " median depth: **", na_called_depth, "x**"),
    paste0(
      "- Lower-to-higher HA/NA depth ratio: **",
      round(ha_na_depth_ratio, 4),
      "**"
    )
  )

  if (host_subtype_depth_imbalance) {
    md <- c(
      md,
      "- The HA and NA depth values are strongly imbalanced. Review for low-level carryover, barcode hopping, contamination, or uneven segment amplification."
    )
  }
}

if (subtype_status == "partial") {
  md <- c(
    md,
    "",
    "## Incomplete subtype evidence",
    paste0("- Supported HA call: **", if (nzchar(ha_call)) ha_call else "none", "**"),
    paste0("- Supported NA call: **", if (nzchar(na_call)) na_call else "none", "**"),
    paste0("- Missing or insufficient segment evidence: **", missing_subtype_segments, "**")
  )
}

if (subtype_status == "mixed") {
  md <- c(
    md,
    "",
    "## Subtype warning",
    paste0("- Qualifying HA candidates: ", if (length(ha$candidates)) paste(ha$candidates, collapse = ", ") else "none"),
    paste0("- Qualifying NA candidates: ", if (length(na$candidates)) paste(na$candidates, collapse = ", ") else "none"),
    "- More than one HA and/or NA subtype met the configured evidence thresholds.",
    "- This may represent coinfection, cross-sample contamination, barcode hopping, or nonspecific assignment and requires review."
  )
}

if (length(ha$low_level) > 0 || length(na$low_level) > 0) {
  low <- character()
  if (length(ha$low_level) > 0) low <- c(low, paste0("HA: ", paste(ha$low_level, collapse = ", ")))
  if (length(na$low_level) > 0) low <- c(low, paste0("NA: ", paste(na$low_level, collapse = ", ")))
  md <- c(md, "", "## Low-level subtype signals", paste0("- ", low))
}

if (length(pass_segments) > 0) {
  md <- c(md, "", "## PASS segments", paste0("- ", pass_segments))
}

if (
  nrow(subtype_evidence) > 0 &&
  all(c("subtype", "segment", "median_depth", "breadth_at_threshold") %in% names(subtype_evidence))
) {
  display_evidence <- subtype_evidence[
    !is.na(subtype_evidence$subtype) &
      nzchar(trimws(as.character(subtype_evidence$subtype))) &
      toupper(trimws(as.character(subtype_evidence$subtype))) != "NA" &
      !is.na(subtype_evidence$segment) &
      nzchar(trimws(as.character(subtype_evidence$segment))),
    ,
    drop = FALSE
  ]

  if (nrow(display_evidence) > 0) {
    md <- c(md, "", "## IRMA subtype evidence")

    for (i in seq_len(nrow(display_evidence))) {
      breadth_value <- suppressWarnings(as.numeric(display_evidence$breadth_at_threshold[i]))
      breadth_percent <- if (is.na(breadth_value)) "NA" else paste0(round(breadth_value * 100, 1), "%")

      md <- c(
        md,
        paste0(
          "- **", display_evidence$subtype[i], "** (", display_evidence$segment[i],
          "): median depth ", display_evidence$median_depth[i], "x; breadth at ",
          min_median_depth, "x = ", breadth_percent
        )
      )
    }
  }
}

if (nrow(top_hits) > 0) {
  md <- c(md, "", "## Top BLAST hits")

  for (i in seq_len(nrow(top_hits))) {
    seg <- if ("qseqid" %in% names(top_hits)) top_hits$qseqid[i] else "NA"
    hit <- if ("sseqid" %in% names(top_hits)) top_hits$sseqid[i] else "NA"
    bits <- if ("bitscore" %in% names(top_hits)) top_hits$bitscore[i] else "NA"
    md <- c(md, paste0("- **", seg, "** -> **", hit, "** (bitscore: ", bits, ")"))
  }
} else {
  md <- c(md, "", "No BLAST hits were reported.")
}

md <- c(
  md,
  "",
  "## Variant summary",
  paste0("- Amino acid changes: **", aa_changes, "**"),
  paste0("- Synonymous: **", syn, "**"),
  paste0("- Nonsynonymous: **", nonsyn, "**"),
  paste0("- Stop gained: **", stop_gained, "**"),
  paste0("- Stop lost: **", stop_lost, "**")
)

writeLines(md, out_md)

if (!file.exists(out_tsv)) stop("Summary TSV was not created: ", out_tsv)
if (!file.exists(out_md)) stop("Summary Markdown was not created: ", out_md)

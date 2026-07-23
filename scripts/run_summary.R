#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

###############################################################################
# Argument parsing
###############################################################################


summary_files <- list.files(
  pattern = "_summary\\.tsv$",
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No summary TSV files supplied.")
}

###############################################################################
# Output directories
###############################################################################

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Read all sample summaries
###############################################################################

message("Reading ", length(summary_files), " sample summaries...")

summary <-
  map_dfr(
    summary_files,
    ~ read_tsv(
      .x,
      show_col_types = FALSE,
      progress = FALSE
    )
  )
  
###############################################################################
# Basic sanity checks
###############################################################################

required_columns <- c(
  "sample_id",
  "host",
  "potential_subtype",
  "subtype_status",
  "pass_segments",
  "coverage_segments",
  "host_subtype_warning",
  "host_subtype_depth_imbalance",
  "ha_called_median_depth",
  "na_called_median_depth",
  "aa_changes"
)

missing_cols <- setdiff(required_columns, names(summary))

if (length(missing_cols) > 0) {

  stop(
    "Missing required columns:\n",
    paste(missing_cols, collapse = ", ")
  )

}

###############################################################################
# Derived fields
###############################################################################
  
  summary <-
  summary |>
  mutate(
    genome_status = case_when(
      pass_segments == 8 ~ "Complete",
      pass_segments >= 6 ~ "Near-complete",
      pass_segments >= 1 ~ "Partial",
      TRUE ~ "Failed"
    ),
    review_required =
      host_subtype_warning |
      host_subtype_depth_imbalance |
      subtype_status != "single" |
      pass_segments < 8
  )
  
  summary$genome_status <-
  factor(
    summary$genome_status,
    levels = c(
      "Complete",
      "Near-complete",
      "Partial",
      "Failed"
    )
  )

###############################################################################
# Review reasons
###############################################################################

review <-
  summary |>
  rowwise() |>
  mutate(

    review_reason =
      paste(

        c(

          if (host_subtype_warning)
            "Host/subtype discordance",

          if (host_subtype_depth_imbalance)
            "HA/NA depth imbalance",

          if (subtype_status == "mixed")
            "Mixed subtype",

          if (subtype_status == "partial")
            "Partial subtype",

          if (subtype_status == "undetermined")
            "Undetermined subtype",

          if (pass_segments < 8)
            paste(pass_segments, "PASS segments")

        ),

        collapse = "; "

      )

  ) |>
  ungroup() |>
  filter(review_required)

###############################################################################
# Write summary tables
###############################################################################

write_tsv(
  summary,
  "run_summary.tsv"
)

write_tsv(
  review,
  "tables/samples_requiring_review.tsv"
)

message("Wrote run_summary.tsv")

###############################################################################
# Common plotting theme
###############################################################################

theme_pipeline <- function() {

  theme_bw(base_size = 12) +

    theme(

      panel.grid.minor = element_blank(),

      plot.title =
        element_text(face = "bold"),

      legend.position = "right",

      strip.background =
        element_rect(fill = "grey95")

    )

}

###############################################################################
# Panel A - Influenza subtype distribution
###############################################################################

subtype_summary <-
  summary |>
  count(
    potential_subtype,
    sort = TRUE
  )

p_subtype <-
  ggplot(
    subtype_summary,
    aes(
      x = reorder(
        potential_subtype,
        n
      ),
      y = n
    )
  ) +
  geom_col(fill = "#4472C4") +
  coord_flip() +
  labs(
    title = "Influenza subtype",
    x = NULL,
    y = "Samples"
  ) +
  theme_pipeline()

###############################################################################
# Panel B - PASS segments
###############################################################################

p_pass <-
  ggplot(summary, aes(factor(pass_segments))) +
  geom_bar(fill = "#70AD47") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.25
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "PASS segments",
    x = "PASS segments",
    y = "Samples"
  ) +
  theme_pipeline()

###############################################################################
# Panel C - Review categories
###############################################################################

warning_df <-
  tibble(
    category = c(
      "Host/subtype",
      "HA/NA imbalance",
      "Mixed subtype",
      "Partial subtype",
      "<8 PASS"
    ),
    count = c(
      sum(summary$host_subtype_warning,
          na.rm = TRUE),
      sum(summary$host_subtype_depth_imbalance,
          na.rm = TRUE),
      sum(summary$subtype_status == "mixed",
          na.rm = TRUE),
      sum(summary$subtype_status == "partial",
          na.rm = TRUE),
      sum(summary$pass_segments < 8,
          na.rm = TRUE)
    )
  )

p_warning <-
  ggplot(
    warning_df,
    aes(
      reorder(category, count),
      count
    )
  ) +
  geom_col(fill = "#C0504D") +
  coord_flip() +
  labs(
    title = "Review categories",
    x = NULL,
    y = "Samples"
  ) +
  theme_pipeline()

###############################################################################
# Panel D - HA versus NA depth
###############################################################################

scatter_df <-
  summary |>
  filter(
    !is.na(ha_called_median_depth),
    !is.na(na_called_median_depth),
    ha_called_median_depth > 0,
    na_called_median_depth > 0
  )

p_depth <-
  ggplot(
    scatter_df,
    aes(
      ha_called_median_depth,
      na_called_median_depth,
      colour = genome_status,
      shape = review_required
    )
  ) +
  geom_point(
  size = 3,
  alpha = 0.8
) +
geom_abline(
  slope = 1,
  intercept = 0,
  linetype = 2,
  linewidth = 0.4,
  colour = "grey80"
) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  labs(
    title = "HA vs NA depth",
    x = "HA median depth",
    y = "NA median depth"
  ) +
  theme_pipeline()

###############################################################################
# Panel E - Amino acid changes
###############################################################################

p_variants <-
  ggplot(
    summary,
    aes(
      genome_status,
      aa_changes,
      fill = genome_status
    )
  ) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Variant burden",
    x = NULL,
    y = "Amino acid changes"
  ) +
  theme_pipeline() +
  theme(
    legend.position = "none"
  )

###############################################################################
# Combine figure
###############################################################################

combined <-
  (p_subtype | p_pass) /
  (p_warning | p_depth) /
  p_variants +
  plot_annotation(
    title = "APGAP Influenza Pipeline Run Summary",
    tag_levels = "A"
  )

ggsave(
  filename = "figures/run_summary.png",
  plot = combined,
  width = 12,
  height = 12,
  dpi = 300
)

ggsave(
  filename = "figures/run_summary.pdf",
  plot = combined,
  width = 12,
  height = 12
)

message("Wrote figures/run_summary.png")
message("Wrote figures/run_summary.pdf")

###############################################################################
# Markdown report
###############################################################################

overview <- c(
  "# APGAP Influenza Pipeline Run Summary",
  "",
  "## Overview",
  "",
  "| Metric | Value |",
  "|--------|------:|",
  sprintf("| Samples processed | %d |", nrow(summary)),
  sprintf("| Complete genomes (8 PASS) | %d |",
          sum(summary$pass_segments == 8, na.rm = TRUE)),
  sprintf("| Near-complete genomes (6-7 PASS) | %d |",
          sum(summary$pass_segments >= 6 &
                summary$pass_segments < 8,
              na.rm = TRUE)),
  sprintf("| Partial genomes (1-5 PASS) | %d |",
          sum(summary$pass_segments > 0 &
                summary$pass_segments < 6,
              na.rm = TRUE)),
  sprintf("| Failed genomes (0 PASS) | %d |",
          sum(summary$pass_segments == 0,
              na.rm = TRUE)),
  ""
)

figure_section <- c(
  "## Run summary",
  "",
  "![](figures/run_summary.png)",
  "",
  "**Figure 1.** Summary of sequencing run quality and subtype results.",
  "",
  "**A.** Influenza subtype distribution.",
  "**B.** Distribution of PASS genomic segments.",
  "**C.** Review categories.",
  "**D.** HA versus NA median sequencing depth.",
  "**E.** Amino acid variant burden by genome completeness.",
  ""
)

###############################################################################
# Prioritize review list
###############################################################################

review <- review |>
  mutate(
    severity =
  		1000 * host_subtype_warning +
  		100 * host_subtype_depth_imbalance +
   		50 * (subtype_status == "mixed") +
   		25 * (subtype_status == "partial") +
   		10 * (pass_segments < 8)
  	) |>
  arrange(desc(severity), sample_id)

###############################################################################
# Samples requiring immediate review
###############################################################################

review_md <- c(
  "",
  "## Samples requiring immediate review",
  "",
  "| Sample | Why flagged |",
  "|--------|-------------|"
)

if (nrow(review) == 0) {

  review_md <- c(
    review_md,
    "| None | No samples require review. |"
  )

} else {

  review_display <- head(review, 10)

  for (i in seq_len(nrow(review_display))) {

    review_md <- c(
      review_md,
      sprintf(
        "| %s | %s |",
        review_display$sample_id[i],
        review_display$review_reason[i]
      )
    )

  }

  if (nrow(review) > 10) {

    review_md <- c(
      review_md,
      "",
      sprintf(
        "*Only the first 10 of %d samples requiring review are shown. See `tables/samples_requiring_review.tsv` for the complete list.*",
        nrow(review)
      )
    )

  }

}

review_md <- c(review_md, "")


###############################################################################
# Subtype summary
###############################################################################

subtype_md <- c(
  "",
  "## Influenza subtype calls",
  "",
  "| Subtype | Samples |",
  "|----------|-------:|"
)

for (i in seq_len(nrow(subtype_summary))) {

  subtype_md <- c(
    subtype_md,
    sprintf(
      "| %s | %d |",
      subtype_summary$potential_subtype[i],
      subtype_summary$n[i]
    )
  )

}

subtype_md <- c(subtype_md, "")


###############################################################################
# Host summary
###############################################################################

host_summary <-
  summary |>
  count(host, sort = TRUE)

host_md <- c(
  "",
  "## Host categories",
  "",
  "| Host | Samples |",
  "|------|-------:|"
)

for (i in seq_len(nrow(host_summary))) {

  host_md <- c(
    host_md,
    sprintf(
      "| %s | %d |",
      host_summary$host[i],
      host_summary$n[i]
    )
  )

}

host_md <- c(host_md, "")

###############################################################################
# Per-sample summary
###############################################################################

sample_md <- c(
  "",
  "## Sample summary",
  "",
  "| Sample | Host | Subtype  | PASS | Genome  | AA changes | Status |",
  "|--------|------|----------|-----:|---------|-----------:|--------|"
)

summary_table <-
  summary |>
  mutate(
    sample_num = readr::parse_number(sample_id),
    sample_prefix = stringr::str_remove(sample_id, "[0-9]+$")
  ) |>
  arrange(sample_prefix, sample_num) |>
  select(-sample_num, -sample_prefix)

for (i in seq_len(nrow(summary_table))) {

  sample_md <- c(
    sample_md,
    sprintf(
      "| %s | %s | %s | %d | %s | %d | %s |",
      summary_table$sample_id[i],
      summary_table$host[i],
      summary_table$potential_subtype[i],
      summary_table$pass_segments[i],
      summary_table$genome_status[i],
      summary_table$aa_changes[i],
      ifelse(summary_table$review_required[i], "Review", "OK")
    )
  )

}

sample_md <- c(sample_md, "")

###############################################################################
# Run QC summary
###############################################################################

qc_md <- c(
  "",
  "## Run QC summary",
  "",
  "| Metric | Count |",
  "|--------|------:|",
  sprintf("| Samples requiring review | %d |",
          sum(summary$review_required, na.rm = TRUE)),
  sprintf("| Host/subtype warnings | %d |",
          sum(summary$host_subtype_warning, na.rm = TRUE)),
  sprintf("| HA/NA depth imbalances | %d |",
          sum(summary$host_subtype_depth_imbalance, na.rm = TRUE)),
  sprintf("| Mixed subtype calls | %d |",
          sum(summary$subtype_status == "mixed", na.rm = TRUE)),
  sprintf("| Partial subtype calls | %d |",
          sum(summary$subtype_status == "partial", na.rm = TRUE)),
  sprintf("| Undetermined subtype calls | %d |",
          sum(summary$subtype_status == "undetermined", na.rm = TRUE))
)

###############################################################################
# Write Markdown report
###############################################################################

writeLines(
  c(
    overview,
    qc_md,
    review_md,
    figure_section,
    subtype_md,
    host_md,
    sample_md
  ),
  con = "run_summary.md"
)

message("Processed ", nrow(summary), " samples.")
message("Wrote run_summary.md")
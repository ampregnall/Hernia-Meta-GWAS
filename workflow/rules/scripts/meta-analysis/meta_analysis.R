# Fixed-effects inverse-variance weighted meta-analysis with Cochran's Q.
# Adapted from https://github.com/ampregnall/nf-meta-gwas/blob/main/bin/meta_analysis.R
# for Snakemake: inputs come from the rule, and chromosomes are processed in a loop
# (1-22 and 23 = X) instead of one job per chromosome.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(arrow)
library(dplyr)
library(glue)

parquet_files <- unlist(snakemake@input)
output <- snakemake@output[[1]]
trait_type <- snakemake@params[["trait_type"]]
chromosomes <- snakemake@params[["chromosomes"]]

# -----------------------------------------------------------------------------
# Column schema
# -----------------------------------------------------------------------------
CHR_COL <- "CHR"
POS_COL <- "POS"
RSID_COL <- "rsID"
EA_COL <- "EA"
NEA_COL <- "NEA"
EAF_COL <- "EAF"
BETA_COL <- "BETA"
SE_COL <- "SE"
N_COL <- "N"
N_CASE_COL <- "N_CASE"
N_CONTROL_COL <- "N_CONTROL"

by_cols <- c(CHR_COL, POS_COL, RSID_COL, EA_COL, NEA_COL)
eaf_n_col <- paste0("N_", EAF_COL) # temporary column for EAF weighting

if (trait_type == "binary") {
  n_cols_to_sum <- c(EAF_COL, N_CASE_COL, N_CONTROL_COL, N_COL)
} else {
  n_cols_to_sum <- c(EAF_COL, N_COL)
}

# -----------------------------------------------------------------------------
# Open dataset
# -----------------------------------------------------------------------------
message(glue("Meta-analyzing {length(parquet_files)} files:"))
message(paste(parquet_files, collapse = "\n"))

parquet_schema <- arrow::open_dataset(parquet_files[1])$schema

# Override any dictionary-encoded fields to utf8 so all files share one schema
unified_fields <- lapply(parquet_schema$fields, function(field) {
  if (field$type$id == arrow::Type$DICTIONARY) {
    arrow::field(field$name, arrow::utf8())
  } else {
    field
  }
})

ds <- arrow::open_dataset(parquet_files, schema = arrow::schema(unified_fields))

# -----------------------------------------------------------------------------
# IVW meta-analysis, one chromosome at a time
# -----------------------------------------------------------------------------
meta_chromosome <- function(chr) {
  message(glue("Processing chromosome {chr} ({trait_type})"))

  res <- ds |>
    filter(.data[[CHR_COL]] == chr) |>
    filter(
      !is.na(.data[[BETA_COL]]),
      !is.na(.data[[SE_COL]]),
      .data[[SE_COL]] > 0
    ) |>
    mutate(
      W = 1 / (.data[[SE_COL]])^2,
      B = .data[[BETA_COL]] * W,
      WB2 = (B^2) / W,
      # Weight EAF by N before grouping so we can compute a N-weighted mean later
      across(
        any_of(EAF_COL),
        ~ .x * .data[[N_COL]]
      )
    ) |>
    group_by(across(all_of(by_cols))) |>
    summarise(
      N_CONTRIBUTIONS = n(),
      # Accumulate total N used for EAF weighting
      across(
        any_of(EAF_COL),
        ~ sum(if_else(is.na(.x), 0, .data[[N_COL]]), na.rm = TRUE),
        .names = "N_{.col}"
      ),
      across(any_of(POS_COL), ~ min(.x, na.rm = TRUE)),
      WB2 = sum(WB2, na.rm = TRUE),
      across(any_of(c("W", "B")), ~ sum(.x, na.rm = TRUE)),
      across(any_of(n_cols_to_sum), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    collect() |>
    mutate(
      B_w = B,
      B = B / W,
      SE = 1 / sqrt(W),
      # Recover N-weighted mean EAF
      across(
        any_of(EAF_COL),
        ~ .x / .data[[eaf_n_col]]
      )
    ) |>
    select(-any_of(eaf_n_col)) |>
    mutate(
      Z = B / SE,
      P = 2 * stats::pnorm(-abs(Z)),
      MLOG10P = -log10(2) -
        stats::pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE) / log(10),
      # Cochran's Q heterogeneity
      Q = WB2 - (B_w^2) / W,
      Q_DF = pmax(N_CONTRIBUTIONS - 1L, 0L),
      Q_PVAL = stats::pchisq(Q, df = Q_DF, lower.tail = FALSE),
      Q_PVAL = if_else(Q_DF == 0L, 1, Q_PVAL),
      I2 = if_else(Q > 0, pmax((Q - Q_DF) / Q, 0), 0)
    )

  if (nrow(res) == 0) {
    message(glue("No variants remaining after filtering for chromosome {chr}"))
  } else {
    message(glue("Chromosome {chr}: {nrow(res)} variants"))
  }

  res
}

meta_results <- lapply(chromosomes, meta_chromosome) |>
  bind_rows() |>
  arrange(CHR, POS)

# -----------------------------------------------------------------------------
# Final column selection and ordering
# -----------------------------------------------------------------------------
n_cols <- if (trait_type == "binary") {
  c("N_CASE", "N_CONTROL", "N")
} else {
  "N"
}

meta_results <- meta_results |>
  select(
    CHR,
    POS = all_of(POS_COL),
    rsID = all_of(RSID_COL),
    EA = all_of(EA_COL),
    NEA = all_of(NEA_COL),
    BETA = B,
    SE,
    P,
    Z,
    EAF = all_of(EAF_COL),
    N_CONTRIBUTIONS,
    all_of(n_cols),
    Q,
    Q_DF,
    Q_PVAL,
    I2,
    MLOG10P
  ) |>
  mutate(
    SNPID = as.character(glue("{CHR}:{POS}:{NEA}:{EA}")),
    MAF = pmin(EAF, 1 - EAF)
  ) |>
  select(SNPID, everything())

message(glue("Total: {nrow(meta_results)} variants"))

arrow::write_parquet(meta_results, output)

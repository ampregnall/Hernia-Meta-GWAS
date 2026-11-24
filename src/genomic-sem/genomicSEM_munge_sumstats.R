library(data.table)
library(GenomicSEM)
library(tidyverse)

# List all METAL formatted files
sumstats <- c(
  "data/processed/meta-formatted/Ventral-EUR_meta_formatted.txt",
  "data/processed/meta-formatted/Inguinal-EUR_meta_formatted.txt",
  "data/processed/meta-formatted/Umbilical-EUR_meta_formatted.txt",
  "data/processed/meta-formatted/Femoral-EUR_meta_formatted.txt",
  "data/processed/meta-formatted/Diaphragmatic-EUR_meta_formatted.txt"
)

# Process files.
# Note that for binary traits that reflect a meta-analysis across multiple cohorts
# this should reflect the sum of effective sample sizes across contributing cohorts
# or the sample size column in the summary statistics should reflect the SNP-specific,
# sum of effective sample sizes
# Reference: https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Genome%E2%80%90wide-Models#step-1-munge-the-summary-statistics

purrr::walk(sumstats, function(x) {
  message(glue::glue("Loading file: {basename(x)}"))
  prefix <- stringr::str_split_i(basename(x), pattern = "_", i = 1)
  df <- vroom::vroom(x, progress = TRUE) %>%
    dplyr::mutate(MAF = dplyr::if_else(EAF < 0.5, EAF, 1 - EAF)) %>%
    dplyr::mutate(N_EFF = 4 / ((1 / N_CASE) + (1 / N_CONTROL))) %>% # Calculate sum of effective sample size as above
    dplyr::select(
      SNP = rsID,
      A1 = EA,
      A2 = NEA,
      effect = BETA,
      SE,
      P,
      N = N_EFF,
      Z,
      MAF
    )
  vroom::vroom_write(
    df,
    glue::glue("data/processed/genomic-sem/{prefix}_GenomicSEM_input.txt"),
    delim = "\t"
  )
})

# List all adipose trait files
sumstats_adipose <- list(
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_asat_bgen_stats.gz",
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_gfat_bgen_stats.gz",
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_vat_bgen_stats.gz",
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_asatadjbmi3_bgen_stats.gz",
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_vatadjbmi3_bgen_stats.gz",
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_gfatadjbmi3_bgen_stats.gz",
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_asatGfatRatio_bgen_stats.gz",
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_vatAsatRatio_bgen_stats.gz",
  "data/meta/mendelian-randomization/bmi-mvmr-instruments/0321_vatGfatRatio_bgen_stats.gz"
)

n_adipose <- c(
  38965,
  38965,
  38965,
  37641,
  37641,
  37641,
  37641,
  37641,
  37641
)

purrr::walk2(sumstats_adipose, n_adipose, function(x, n) {
  message(glue::glue("Loading file: {basename(x)}"))
  prefix <- stringr::str_split_i(basename(x), pattern = "_", i = 2)
  df <- vroom::vroom(x, progress = TRUE) |>
    mutate(
      N = n,
      Z = BETA / SE,
      MAF = dplyr::if_else(A1FREQ < 0.5, A1FREQ, 1 - A1FREQ)
    ) |>
    select(
      SNP,
      A1 = ALLELE1,
      A2 = ALLELE0,
      effect = BETA,
      SE,
      P = P_BOLT_LMM,
      N,
      Z,
      MAF
    )
  vroom::vroom_write(
    df,
    glue::glue("data/processed/genomic-sem/{prefix}_GenomicSEM_input.txt"),
    delim = "\t"
  )
})


# List files from previous step
sumstats <- list.files(
  "data/processed/genomic-sem/input",
  pattern = "*_input.txt",
  full.names = TRUE
)

# Define hyperparameters
hm3 <- "/project/voltron/Resources/eur_w_ld_chr/w_hm3.snplist"

trait.names <- str_split_i(basename(sumstats), "_", 1)

N <- rep(NA, length(sumstats))
info.filter <- 0.9 # Imputation quality
maf.filter <- 0.01

# Munge sumstats
GenomicSEM::munge(
  files = sumstats,
  hm3 = hm3,
  trait.names = trait.names,
  N = N,
  info.filter = info.filter,
  maf.filter = maf.filter
)

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(data.table)
library(tidyverse)

# Load custom scripts
source(snakemake@params[["creds"]])
source(snakemake@params[["locus"]])

lead <- fread(snakemake@input[["lead"]])

if (nrow(lead) == 0) {
  credset <- data.frame(
    SNPID = character(),
    rsID = character(),
    LOCUS = character(),
    CHR = numeric(),
    POS = numeric(),
    EA = character(),
    NEA = character(),
    BETA = numeric(),
    SE = numeric(),
    NEAREST_GENE = character(),
    BF = numeric(),
    BF_PIP = numeric()
  )
} else {
  sumstats <- fread(snakemake@input[["sumstats"]])

  sumstats <- locus_extract(
    sumstats,
    sumstats_chr_col = CHR,
    sumstats_pos_col = POS,
    locus_df = lead,
    locus_chr_col = CHR,
    locus_pos_col = POS,
    locus_gene_col = GENE,
    locus_size = snakemake@params[["locus_size"]]
  )

  credset <- calc_credset(
    sumstats,
    locus_marker_col = locus_marker,
    effect_col = BETA,
    se_col = SE,
    samplesize_col = N,
    cred_interval = snakemake@params[["credible_interval"]]
  )

  credset <- credset %>%
    dplyr::select(
      SNPID,
      rsID,
      LOCUS = locus_marker,
      CHR = chromosome,
      POS = position,
      EA,
      NEA,
      BETA,
      SE,
      NEAREST_GENE = gene,
      BF = bf,
      BF_PIP = posterior_prob
    )
}

fwrite(credset, snakemake@output[["credset"]], sep = "\t")

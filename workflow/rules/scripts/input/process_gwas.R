# Clean, lift over, and harmonize one set of summary statistics. Writes a table of
# variant counts at each step.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Retrieve Snakemake objects
input <- snakemake@input[["sumstats"]]
dbsnp_tarball <- snakemake@input[["dbsnp_tarball"]]
output <- snakemake@output[["sumstats"]]
counts_out <- snakemake@output[["counts"]]
sumstats_id <- snakemake@wildcards[["id"]]
utils <- snakemake@params[["utils"]]
chainfile <- snakemake@params[["chainfile"]]
study <- snakemake@params[["study"]]
n_cases <- snakemake@params[["n_cases"]]
n_controls <- snakemake@params[["n_controls"]]

# Load required libraries
library(dplyr)
library(data.table)
library(MungeSumstats)
library(BSgenome.Hsapiens.NCBI.GRCh38)

# Load custom functions
source(utils)

count_variants <- function(df, step) {
  chrom <- as.character(df[[intersect(c("CHROM", "CHR"), names(df))[1]]])
  data.table(
    id = sumstats_id,
    step = step,
    n_variants = nrow(df),
    n_chrX = sum(chrom %in% c("X", "23"))
  )
}

# Load summary statistics
df <- fread(input)
counts <- list(data.table(
  id = sumstats_id,
  step = "input",
  n_variants = nrow(df),
  n_chrX = NA_integer_
))

# Study-specific cleaning (column names, liftover, missing values, chromosomes)
if (study == "Pan-UKBB") {
  df <- liftover_37_to_38(df, chr = chr, pos = pos, chainfile)
  df <- process_ukbb(df, n_cases, n_controls)
} else if (study == "FinnGen-R12") {
  df <- process_finngen(df, n_cases, n_controls)
} else if (study == "AOU") {
  df <- process_aou(df, n_cases, n_controls)
} else if (study == "MVP") {
  df <- process_mvp(df, n_cases, n_controls)
} else if (study == "BJP") {
  df <- liftover_37_to_38(df, chr = CHR, pos = POS, chainfile)
  df <- process_bjp(df, n_cases, n_controls)
} else {
  stop("No cleaning function for study: ", study)
}
counts[[length(counts) + 1]] <- count_variants(df, "study_cleaning")

# Harmonize to GRCh38. X is kept (MungeSumstats drops X, Y, and MT by default).
# rsIDs come from the local dbSNP SNPlocs tarball; `dbSNP` is ignored when a
# tarball is supplied.
df <- MungeSumstats::format_sumstats(
  df,
  ref_genome = "GRCh38",
  dbSNP_tarball = dbsnp_tarball,
  rmv_chr = c("Y", "MT"),
  bi_allelic_filter = FALSE,
  # Flip FRQ with BETA when alleles are swapped to match the reference, so EAF stays
  # the frequency of the effect allele (non-biallelic SNPs are not flipped)
  allele_flip_frq = TRUE,
  return_data = TRUE,
  return_format = "data.table"
)
counts[[length(counts) + 1]] <- count_variants(df, "mungesumstats")

# Variants removed at each step (the minor allele count filter is applied after
# LDSC correction and its count is added there)
counts <- data.table::rbindlist(counts)
counts[, n_removed := shift(n_variants) - n_variants]

fwrite(df, output, sep = "\t")
fwrite(counts, counts_out, sep = "\t")

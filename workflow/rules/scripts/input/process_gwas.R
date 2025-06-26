# Retrieve Snakemake objects
input <- snakemake@input[[1]]
output <- snakemake@output[[1]]
utils <- snakemake@params[["utils"]]
chainfile <- snakemake@params[["chainfile"]]
study <- snakemake@params[["study"]]
genome <- snakemake@params[["genome"]]
n_cases <- snakemake@params[["n_cases"]]
n_controls <- snakemake@params[["n_controls"]]

# Load required libraries
library(readr)
library(dplyr)
library(data.table)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(BSgenome.Hsapiens.NCBI.GRCh38)

# Load custom functions
source(utils)

# Load summary statistics
df <- fread(input)

# # Perform liftover on relevant files
# if (genome == "GrCh37") {
#   df <- liftover_37_to_38(df, chainfile)
# }

# Process studies
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
}

df.formatted <- MungeSumstats::format_sumstats(
  df, 
  ref_genome = "GrCh38", 
  dbSNP = 144,
  bi_allelic_filter = FALSE, 
  allele_flip_frq = FALSE,
  return_data = TRUE,
  return_format = "data.table")

fwrite(df.formatted, output, sep = "\t")
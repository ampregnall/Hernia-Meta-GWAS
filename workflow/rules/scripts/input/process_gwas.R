# Retrieve Snakemake objects
input <- snakemake@input
output <- snakemake@output
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

getwd()

# Load custom functions
source(utils)

# Load summary statistics
df <- fread(input)

# Perform liftover on relevant files
if (genome == "GrCh37") {
  df <- liftover_37_to_38(df, chainfile)
}

# Process studies
if (study == "Pan-UKBB") {
  df <- process_ukbb(df, n_cases, n_controls)
} else if (study == "FinnGen-R12") {
  df <- process_finngen(df, n_cases, n_controls)
} else if (study == "AOU") {
  df <- process_aou(df, n_cases, n_controls)
} else if (study == "MVP") {
  df <- process_mvp(df, n_cases, n_controls)
}

# Save data
fwrite(df, output, sep = "\t")

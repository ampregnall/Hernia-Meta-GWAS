# Load required libraries
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(dtplyr)
source("workflow/rules/scripts/input/process_gwas_utils.R")

# Command line input
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
outfile <- basename(input_file)
study <- word(basename(input_file), 1, sep = "_")
genome <- word(basename(input_file), 4, sep = "_")
genome <- str_remove(genome, ".txt")

# Load data
df <- fread(input_file)

if (study == "Pan-UKBB") {
  df <- process_ukbb(df)
} else if (study == "FinnGen-R12") {
  df <- process_finngen(df)
} else if (study == "AOU") {
  df <- process_aou(df)
} else if (study == "MVP") {
  df <- process_mvp(df)
}

if (genome == "GrCh37") {
  df <- liftover_37_to_38(df)
  outfile <- paste0(str_remove(outfile, "GrCh37.txt"), "GrCh38.txt")
}

# Save data
output_file <- paste0("data/processed/gwas-summary-statistics/", outfile)
fwrite(df, output_file, sep = "\t")
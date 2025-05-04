# Load required libraries
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(dtplyr)
source("/workflow/rules/scripts/process_gwas_utils.R")

# Command line input
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
study <- stringr::word(basename(input_file), 1, sep = "_")
print(input_file)

# Load data
df <- fread(input_file) %>% lazy_dt()

if (study == "Pan-UKBB") {
  df <- process_ukbb(df)
} else if (study == "FinnGen-R12") {
  df <- process_finngen(df)
} else if (study == "AOU") {
  df <- process_aou(df)
} else if (study == "MVP") {
  df <- process_mvp(df)
}

# Save data
output_file <- paste0("data/processed/gwas-summary-statistics/", study)
fwrite(df, output_file, sep = "\t")
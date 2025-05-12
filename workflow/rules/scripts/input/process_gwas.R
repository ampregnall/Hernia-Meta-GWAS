# Load required libraries
library(readr)
library(dplyr)
library(data.table)
source("workflow/rules/scripts/input/process_gwas_utils.R")
chainfile <- "data/meta/liftover/hg19ToHg38.over.chain"

# Command line input
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

# Extract study parameters
study <- stringr::word(basename(input_file), 1, sep = "_")
phenotype <- stringr::word(basename(input_file), 2, sep = "_")
ancestry <- stringr::word(basename(input_file), 3, sep = "_")
genome <- stringr::word(basename(input_file), 4, sep = "_")
genome <- stringr::str_remove(genome, ".txt")
outfile <- basename(input_file)

# Load metadata
meta <- readxl::read_xlsx("data/meta/studies/hernia_gwas_summary_stats.xlsx", sheet = "data") %>% 
  filter(Study == study & Phenotype == phenotype & Population == ancestry)

# Get number of cases and control
n_cases = meta$Cases
n_controls = meta$Controls

# Load summary statistics
df <- fread(input_file)

if (genome == "GrCh37") {
  df <- liftover_37_to_38(df, chainfile)
  outfile <- paste0(stringr::str_remove(outfile, "GrCh37.txt"), "GrCh38.txt")
}

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
output_file <- paste0("data/processed/gwas-summary-statistics/", outfile)
fwrite(df, output_file, sep = "\t")
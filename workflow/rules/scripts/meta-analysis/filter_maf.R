library(data.table)
library(dplyr)

# Command line input
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
outfile <- basename(input_file)

# Load summary statistics
df <- fread(input_file)

# Process files
df <- df %>% 
  dplyr::mutate(MAF = dplyr::if_else(EFF_ALL_FREQ < 0.5, EFF_ALL_FREQ, 1 - EFF_ALL_FREQ)) %>%
  dplyr::filter(MAF > 50 / (2* N_CASES))

# Save data
output_file <- paste0("data/processed/gwas-inspector-sumstats/", outfile)
fwrite(df, output_file, sep = "\t")
input <- snakemake@input[[1]]
output <- snakemake@output[[1]]

# Load libraries
library(data.table)
library(dplyr)

# Load summary statistics
df <- fread(input)

# Process files
df <- df %>%
  dplyr::mutate(MAF = dplyr::if_else(EAF < 0.5, EAF, 1 - EAF)) %>%
  dplyr::filter(MAF > 50 / (2 * N_CASE))

# Save data
fwrite(df, output, sep = "\t")

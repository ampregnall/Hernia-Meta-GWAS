input <- snakemake@input[[1]]
output <- snakemake@output[[1]]
utils <- snakemake@params[["utils"]]

# Load libraries
library(data.table)
library(dplyr)

# Load summary statistics
df <- fread(input)

# Process files
df <- df %>%
  dplyr::mutate(MAF = dplyr::if_else(EFF_ALL_FREQ < 0.5, EFF_ALL_FREQ, 1 - EFF_ALL_FREQ)) %>%
  dplyr::filter(MAF > 50 / (2 * N_CASES))

# Save data
fwrite(df, output, sep = "\t")

library(data.table)
library(tidyverse)

# Load custom scripts
source(snakemake@params[["locus"]])

sumstats <- fread(snakemake@input[["sumstats"]])
lead <- fread(snakemake@input[["lead"]])
rs <- fread(snakemake@input[["rsparsepro"]]) 
rs <- rs %>% dplyr::select(rsID = RSID, RSP_PIP = PIP, RSP_CS = cs)

# Extract all variants within 500 kb of lead variants
sumstats <- locus_extract(sumstats,
                          sumstats_chr_col = CHR, 
                          sumstats_pos_col = POS, 
                          locus_df = lead,
                          locus_chr_col = CHR,
                          locus_pos_col = POS,
                          locus_gene_col = GENE,
                          locus_size = 1e6)

sumstats <- sumstats %>% dplyr::filter(rsID %in% rs$rsID) %>%
  dplyr::select(SNPID, rsID, LOCUS = locus_marker, CHR = chromosome, POS = position,
                EA, NEA, BETA, SE, NEAREST_GENE = gene) %>% 
  dplyr::left_join(rs, by = "rsID")

fwrite(sumstats, snakemake@output[["finemapping"]], sep = "\t")
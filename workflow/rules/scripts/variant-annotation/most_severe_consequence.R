library(data.table)
library(tidyverse)
source(snakemake@params[["otarget"]])

df <- fread(snakemake@input[["credset"]])

snps <- df %>% dplyr::pull(SNPID) %>%
  stringr::str_replace_all(., ":", "_")

snps_annotated <- annotateVariantConsequences(snps)

snps_annotated <- snps_annotated %>% 
  dplyr::mutate(SNPID = stringr::str_replace_all(SNPID, "_", ":"))

fwrite(snps_annotated, snakemake@output[["out"]])
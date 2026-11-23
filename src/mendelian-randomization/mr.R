library(TwoSampleMR)
library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
phenotype <- args[1]

# Load instruments
mrinstrument <- as.data.frame(fread("data/meta/mendelian-randomization/MR-Instruments.txt", sep = "\t"))

# Load outcome data
outcome <- TwoSampleMR::read_outcome_data(
  glue::glue("data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt"), 
  sep = "\t",
  snps = mrinstrument$rsID, 
  snp_col = "rsID", 
  beta_col = "BETA", 
  se_col = "SE",
  eaf_col = "EAF", 
  pval_col = "P", 
  samplesize_col = "N", 
  ncase_col = "N_CASE", 
  ncontrol_col = "N_CONTROL", 
  chr_col = "CHR", 
  pos_col = "POS", 
  effect_allele_col = "EA", 
  other_allele_col = "NEA")

mrinstrument <- mrinstrument |> 
  dplyr::filter(SNP %in% outcome$SNP & !is.na(SNP)) |> 
  dplyr::group_split(Phenotype) 

mrinstrument <- purrr::map(mrinstrument, function(df) {
  
  df <- dplyr::rename(df, rsid = SNP)
  
  df_clump <- try(
    ieugwasr::ld_clump(
      dat       = df,
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile     = "/project/voltron/Resources/1000-Genomes-Phase3/MRC-IEU/EUR"
    )
  )
  
  # On error or no rows, return empty tibble with correct columns
  if (inherits(df_clump, "try-error") || is.null(df_clump) || nrow(df_clump) == 0L) {
    df_clump <- df
  } 
  return(df_clump)
  
}) |> dplyr::bind_rows()

mrinstrument <- dplyr::select(mrinstrument, -c(id))

# Format instrument as TwoSampleME exposure
exposure <- TwoSampleMR::format_data(
  mrinstrument, 
  type = "exposure", 
  phenotype_col = "Phenotype", 
  snp_col = "rsid",
  beta_col = "beta", 
  se_col = "se", 
  eaf_col = "eaf", 
  effect_allele_col = "effect_allele", 
  other_allele_col = "other_allele",
  pval_col = "pval", 
  samplesize_col = "samplesize"
)

dat <- TwoSampleMR::harmonise_data(exposure, outcome)

res <- mr(dat)
res$fdr <- p.adjust(res$pval, method = "fdr")
res_het <- mr_heterogeneity(dat)
res_pleio <- mr_pleiotropy_test(dat)

fwrite(res, glue::glue("data/processed/mendelian-randomization/{phenotype}_mr.txt"), sep = "\t")
fwrite(res_het, glue::glue("data/processed/mendelian-randomization/{phenotype}_mr_heterogeneity.txt"), sep = "\t")
fwrite(res_pleio, glue::glue("data/processed/mendelian-randomization/{phenotype}_mr_pleiotropy.txt"), sep = "\t")
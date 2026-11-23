library(TwoSampleMR)
library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
phenotype <- args[1]

# Load instruments
mrinstrument <- TwoSampleMR::read_exposure_data(
  filename = "data/meta/mendelian-randomization/DrugBank_MR_Instruments.txt", 
  phenotype_col = "Phenotype", 
  snp_col = "rsID", 
  beta_col = "BETA", 
  se_col = "SE",
  effect_allele_col = "EA", 
  other_allele_col = "NEA", 
  pval_col = "P", 
  pos_col = "POS", 
  chr_col = "CHR",
  sep = "\t"
)

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

dat <- TwoSampleMR::harmonise_data(mrinstrument, outcome)

res <- mr(dat)
fwrite(res, glue::glue("data/processed/mendelian-randomization/{phenotype}_DrugTargets_mr.txt"), sep = "\t")

res_het <- mr_heterogeneity(dat)
fwrite(res_het, glue::glue("data/processed/mendelian-randomization/{phenotype}_DrugTargets_mr_heterogeneity.txt"), sep = "\t")

res_pleio <- mr_pleiotropy_test(dat)
fwrite(res_pleio, glue::glue("data/processed/mendelian-randomization/{phenotype}_DrugTargets_mr_pleiotropy.txt"), sep = "\t")
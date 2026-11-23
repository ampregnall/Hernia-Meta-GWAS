library(data.table)
library(GenomicSEM)

### PROCESS BMI SUMMARY STATISTICS WITH GENOMIC-SEM PACKAGE ###
df <- fread("data/meta/lava/BMI_sumstats.txt.gz") # From GCST009004. A1 is effect allele

df <- df |>
  dplyr::mutate(MAF = dplyr::if_else(FRQ < 0.5, FRQ, 1 - FRQ)) |>
  dplyr::mutate(N = 806834, Z = BETA / SE) |>
  dplyr::select(SNP, A1, A2, effect = BETA, SE, P, N, Z, MAF)

vroom::vroom_write(
  df,
  glue::glue("data/processed/genomic-sem/BMI_GenomicSEM_input.txt"),
  delim = "\t"
)

# List files from previous step
sumstats <- "data/processed/genomic-sem/BMI_GenomicSEM_input.txt"

# Define hyperparameters
hm3 <- "/project/voltron/Resources/eur_w_ld_chr/w_hm3.snplist"
trait.names <- c("BMI")
N <- c(NA)
info.filter <- 0.9 # Imputation quality
maf.filter <- 0.01

# Munge sumstats
GenomicSEM::munge(
  files = sumstats,
  hm3 = hm3,
  trait.names = trait.names,
  N = N,
  info.filter = info.filter,
  maf.filter = maf.filter
)

### USE GENOMIC-SEM FOR COVARIANCE ###
sumstats <- list.files(
  "data/processed/genomic-sem/munge",
  pattern = "\\.gz$",
  full.names = TRUE
)


sample.prev <- c(
  NA, # BMI
  0.052510, # Diaphragmatic
  0.003094, # Femoral
  0.061529, # Inguinal
  0.019831, # Umbilical
  0.014628
) # Ventral

population.prev <- c(
  NA, # BMI
  0.00018, # Diaphragmatic
  0.008, # Femoral
  0.04, # Inguinal
  0.03, # Umbilical
  0.017
) # Ventral


trait.names <- c(
  "BMI",
  "Diaphragmatic",
  "Femoral",
  "Inguinal",
  "Umbilical",
  "Ventral"
)

ld <- "/project/voltron/Resources/eur_w_ld_chr/"
wld <- "/project/voltron/Resources/eur_w_ld_chr/"

ldsc_output <- GenomicSEM::ldsc(
  traits = sumstats,
  sample.prev = sample.prev,
  population.prev = population.prev,
  ld = ld,
  wld = wld,
  trait.names = trait.names
)

cov <- ldsc_output$I
rownames(cov) <- trait.names
colnames(cov) <- trait.names
cov <- round(cov2cor(cov), 5)
lava_phenotype_order <- c(
  "Ventral",
  "Inguinal",
  "Umbilical",
  "Femoral",
  "Diaphragmatic",
  "BMI"
)
cov <- cov[lava_phenotype_order, lava_phenotype_order]

write.table(cov, "data/processed/lava/LAVA_sample_overlap.txt", quote = F)

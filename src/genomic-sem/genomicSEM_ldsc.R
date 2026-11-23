library(data.table)
library(GenomicSEM)
library(psych)

# -----------------------
# Inputs
# -----------------------
sumstats <- list.files("data/processed/genomic-sem/munge", pattern = "\\.gz$", full.names = TRUE)

trait.names <- c(
  "ASAT",
  "ASATAdjBMI",
  "ASATGFATRatio",
  "BMI",
  "Diaphragmatic",
  "Femoral", 
  "GFAT",
  "GFATAdjBMI",
  "Inguinal",
  "Umbilical",
  "VAT",
  "VATAdjBMI",
  "VATASATRatio",
  "VATGFATRatio",
  "Ventral")

sample.prev <- c(
  NA,
  NA,
  NA,
  NA,
  0.50, # 50% for N_EFF
  0.50,
  NA,
  NA,
  0.50,
  0.50,
  NA,
  NA,
  NA,
  NA,
  0.50
)
 
population.prev <- c(
  NA,
  NA,
  NA,
  NA,
  0.00018, # Diaphragmatic 
  0.008,   # Femoral
  NA,
  NA,
  0.04,    # Inguinal
  0.03,    # Umbilical
  NA,
  NA,
  NA,
  NA,
  0.017
  )   # Ventral



ld  <- "/project/voltron/Resources/eur_w_ld_chr/"
wld <- "/project/voltron/Resources/eur_w_ld_chr/"

out_ldsc <- "data/processed/genomic-sem/ldsc/GenomicSEM_hernia_adipose_traits.RData"

# -----------------------
# 1) LDSC
# -----------------------
ldsc_output <- GenomicSEM::ldsc(
  traits = sumstats,
  sample.prev = sample.prev,
  population.prev = population.prev,
  ld = ld,
  wld = wld,
  trait.names = trait.names
)

save(ldsc_output, file = out_ldsc)

# Genetic covariance matrix for hernia subtypes
hernia_traits <- c("Ventral", "Inguinal", "Umbilical", "Femoral", "Diaphragmatic")

S <- ldsc_output$S
colnames(S) <- rownames(S) <- trait.names
S_hernia <- S[hernia_traits, hernia_traits]

# -----------------------
# 2) EFA using factanal (1-2 factors; varimax + promax)
# -----------------------
run_factanal <- function(S, nfactors, rotation, n.obs) {
  # factanal expects covmat as either a matrix or list(cov=..., n.obs=...)
  # and standardizes the covariance matrix internally
  fa <- stats::factanal(
    covmat   = list(cov = S, n.obs = n.obs),
    factors  = nfactors,
    rotation = rotation
  )
  fa
}

n_obs_efa <- 100000L

efa_solutions <- list()
for (k in 1:2) {
  for (rot in c("varimax", "promax")) {
    key <- paste0("k", k, "_", rot)
    efa_solutions[[key]] <- run_factanal(S_hernia, nfactors = k, rotation = rot, n.obs = n_obs_efa)
    cat("\n====================\nEFA:", key, "\n====================\n")
    print(efa_solutions[[key]]$loadings)
  }
}

efa_loadings <- purrr::imap(efa_solutions, function(fa, key) {
  L <- unclass(fa$loadings)
  data.frame(model = key, trait = rownames(L), L, row.names = NULL)
}) |> dplyr::bind_rows()

fwrite(efa_loadings, "data/processed/genomic-sem/EFA_loadings.txt", sep = "\t")

# -----------------------
# 3) CFA in Genomic-SEM
# -----------------------
# usermodel() only uses the traits named in the model, so the adipose traits in
# ldsc_output are ignored except where they are explicitly included

# Common factor model
common_factor <- '
  CF =~ Ventral + Inguinal + Umbilical + Femoral + Diaphragmatic
'

# Two factor model. F2 has two indicators, so its loadings are constrained equal
# for identification. Residual variances for ventral and femoral hernia are
# constrained > 0.001 to avoid negative (Heywood) estimates.
two_factor_constrained <- '
  F1 =~ Ventral + Umbilical + Diaphragmatic
  F2 =~ a*Inguinal + a*Femoral
  Ventral ~~ c*Ventral
  Femoral ~~ d*Femoral
  c > .001
  d > .001
'

# Two factor model with F2 regressed on BMI, VAT, and ASAT
two_factor_constrained_all <- '
  F1 =~ Ventral + Umbilical + Diaphragmatic
  F2 =~ a*Inguinal + a*Femoral
  Ventral ~~ z*Ventral
  Femoral ~~ x*Femoral
  F2 ~ BMI + VAT + ASAT
  z > .001
  x > .001
'

models <- list(
  common_factor = common_factor,
  two_factor_constrained = two_factor_constrained,
  two_factor_constrained_all = two_factor_constrained_all
)

fits <- purrr::map(models, function(model) {
  GenomicSEM::usermodel(
    covstruc = ldsc_output,
    model = model,
    estimation = "DWLS",
    std.lv = TRUE
  )
})

cfa_modelfit <- purrr::imap(fits, ~ data.frame(model = .y, .x$modelfit)) |> dplyr::bind_rows()
cfa_results <- purrr::imap(fits, ~ data.frame(model = .y, .x$results)) |> dplyr::bind_rows()

print(cfa_modelfit)

fwrite(cfa_modelfit, "data/processed/genomic-sem/CFA_modelfit.txt", sep = "\t")
fwrite(cfa_results, "data/processed/genomic-sem/CFA_results.txt", sep = "\t")

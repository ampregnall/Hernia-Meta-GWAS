# ---- RETICULATE SETUP ----
# Prefer a Python explicitly given by the job environment; if not, fall back to the
# Python on PATH (works great when the Snakemake rule uses `conda:`).
py_exec <- Sys.getenv("RETICULATE_PYTHON", unset = "")

if (py_exec == "") {
  # Use Python from the current shell's PATH (e.g., Snakemake conda env)
  py_exec <- unname(Sys.which("python"))
  if (is.na(py_exec) || py_exec == "") {
    stop("No Python found. Set RETICULATE_PYTHON or ensure python is on PATH.")
  }
}

# Sanity check before binding reticulate
if (!file.exists(py_exec)) {
  stop(sprintf("Configured Python does not exist: %s", py_exec))
}

# Lock reticulate to that interpreter
Sys.setenv(RETICULATE_PYTHON = py_exec)
library(reticulate)
reticulate::use_python(py_exec, required = TRUE)
cfg <- reticulate::py_config()
message("Using Python: ", cfg$python, " (", cfg$version, ")")

library(bigsnpr)
library(tidyverse)
library(reticulate)
library(logger)

# Extract snakemake objects
sumstats_path <- snakemake@input[["sumstats"]]
lead_path <- snakemake@input[["lead"]]
locus_utils <- snakemake@params[["utils"]]
rsparse <- snakemake@params[["rsparsepro"]]
stats_output <- snakemake@output[["stats"]]
rsparsepro_output <- snakemake@output[["rsparsepro"]]
bed <- snakemake@input[["bed"]]

# Load custom scripts
source(locus_utils)
reticulate::source_python(rsparse)

# Load bed file
logger::log_info("Loading bed file")
obj <- bigsnpr::snp_attach(bed)
G <- obj$genotypes

logger::log_info("Loading summary statistics")
sumstats <- data.table::fread(sumstats_path)

logger::log_info("Loading lead variants")
lead <- data.table::fread(lead_path)

windowed <- locus_extract(sumstats_df = sumstats,
                          sumstats_chr_col = CHR,
                          sumstats_pos_col = POS,
                          locus_df = lead,
                          locus_chr_col = CHR,
                          locus_pos_col = POS,
                          locus_gene_col = GENE,
                          locus_size = 1e6) # CHECK LOCUS SIZE
  
# Get marker names
markers <- windowed %>% dplyr::pull(locus_marker) %>% unique()
  
snps <- purrr::map(markers, function(marker) {
    
    # Get Z scores for RSparsePro
    snps <- windowed %>% 
      dplyr::filter(locus_marker == marker) %>%
      dplyr::select(rsID, Z) %>%
      dplyr::distinct(rsID, .keep_all = TRUE) %>%
      dplyr::arrange(rsID) %>%
      dplyr::rename(RSID = rsID)
    
    return(snps)

})

names(snps) <- markers

# Calculate SNP Correlation Matrix ----------------------------------------

logger::log_info("Calculating SNP correlation matrices")
snp_correlations <- purrr::map(snps, function(x) {
  n_snps_locus <- nrow(x)
  snp_sub <- which(obj$map$marker.ID %in% x$RSID)
  corr <- snp_cor(G, ind.col = snp_sub)
  x <- dplyr::filter(x, RSID %in% obj$map$marker.ID)
  n_snps_ref <- nrow(x)
  delta <- n_snps_locus - n_snps_ref
  logdf <- data.frame("n_snps_locus" = n_snps_locus,
                      "n_snps_ref" = n_snps_ref,
                      "delta" = delta)
  return(list("snps" = x, "snp_correlations" = as.matrix(corr), "logdf" = logdf))
})

stats <- snp_correlations %>% 
  purrr::map("logdf") %>%
  dplyr::bind_rows(.id = "LOCUS") 

fwrite(stats, stats_output, sep = "\t")

# Run RSparsePro ----------------------------------------------------------

logger::log_info("Running RSparsePro...")
rsparsepro_results <- purrr::map(test, function(x) {
  calculate_credible_sets(zfile = x[[1]], ld = x[[2]])})

results <- rsparsepro_results %>% dplyr::bind_rows(.id = "LOCUS")

fwrite(results, rsparsepro_output, sep = "\t")
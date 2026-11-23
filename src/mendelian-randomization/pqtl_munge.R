#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MungeSumstats)
  library(glue)
  library(tools)
})

# ---- Parse command line argument ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript run_mungesumstats_pqtl.R <input_file>")
}

input_path <- args[1]

if (!file.exists(input_path)) {
  stop(glue("Input file does not exist: {input_path}"))
}

# ---- Prepare output path ----
out_dir <- "data/meta/pqtl"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

file_name <- basename(input_path)

# Remove .txt.gz or .gz safely
file_stub <- sub("\\.txt\\.gz$", "", file_name)
file_stub <- sub("\\.gz$", "", file_stub)

output_path <- file.path(out_dir, paste0(file_stub, "_GRCh38.tsv.gz"))

message(glue("Input:  {input_path}"))
message(glue("Output: {output_path}"))

# ---- Run MungeSumstats ----
format_sumstats(
  path = input_path,
  ref_genome = NULL,
  convert_ref_genome = "GRCh38",
  save_path = output_path
)

message("Finished successfully.")

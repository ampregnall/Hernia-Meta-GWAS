# Install hyprcoloc and dependencies if needed
if (!requireNamespace("hyprcoloc", quietly = TRUE)) {
  # remotes::install_version("RcppEigen", version = "0.3.3.9.3")
  remotes::install_github("jrs95/hyprcoloc", build_vignettes = TRUE, dependencies = FALSE)
}

# Load packages
library(hyprcoloc)
library(tidygenomics)
library(data.table)
library(tidyverse)

# Load custom scripts
locus <- snakemake@params[["locus"]]
source(locus)

### LOAD SNAKEMAKE OBJECTS
eqtl <- snakemake@params[["eqtl"]]
sumstats <- snakemake@input[["sumstats"]]
lead <- snakemake@input[["lead"]]
phenotype <- stringr::str_split_i(basename(sumstats), "-", 1) # Extract phenotype from filename
min_overlap <- snakemake@params[["min_overlap"]]

### Load eQTL data
eqtl <- data.table::fread(eqtl, sep = "\t")

### Join. Select relevant columns.
eqtl <- eqtl %>%
  dplyr::mutate(trait = glue::glue("{gene_id} eQTL")) %>%
  dplyr::select(trait, rsID = rsid, CHR = chromosome, POS = position, 
                BETA = beta, SE = se, EA = alt, NEA = ref, P = pvalue) %>%
  dplyr::filter(CHR != "X") %>%
  dplyr::mutate(CHR = as.integer(CHR))
  
### Load sumstats
sumstats <- data.table::fread(sumstats) %>%
    dplyr::mutate(trait = phenotype) %>%
    dplyr::select(trait, rsID, CHR, POS, BETA, SE, EA, NEA, P) 

### Load lead variants
lead_vars <- data.table::fread(lead) %>%
  dplyr::mutate(START = POS - 500000, END = POS + 500000) %>%
  tidygenomics::genome_cluster(by=c("CHR", "START", "END"), max_distance = 0) %>% 
  dplyr::slice_min(order_by = P, by=c("cluster_id"))

# Get all data +/- 500kb of lead variants
sumstats_windowed <- sumstats %>%
  locus_extract_dtable(sumstats_chr_col = CHR, 
                sumstats_pos_col = POS, 
                locus_df = lead_vars,
                locus_chr_col = CHR,
                locus_pos_col = POS,
                locus_gene_col = GENE,
                locus_size = 1e6)

# Get all data +/- 500kb of lead variants
eqtl_windowed <- eqtl %>%
  locus_extract_dtable(sumstats_chr_col = CHR, 
                       sumstats_pos_col = POS, 
                       locus_df = lead_vars,
                       locus_chr_col = CHR,
                       locus_pos_col = POS,
                       locus_gene_col = GENE,
                       locus_size = 1e6)

valid_rsid <- sumstats_windowed %>%
  dplyr::distinct(rsID) %>%
  dplyr::pull(rsID)

snps_per_gene <- eqtl_windowed %>%
  dplyr::filter(rsID %in% valid_rsid) %>%
  dplyr::distinct(trait, rsID) %>%
  dplyr::group_by(trait) %>%
  dplyr::summarize(snp_set = list(rsID), .groups="drop") %>%
  dplyr::mutate(n_overlap_gwas = lengths(lapply(snp_set, function(s) intersect(s, valid_rsid)))) %>%
  dplyr::filter(n_overlap_gwas > min_overlap)

# Run Coloc
hyprcoloc_res <- eqtl_windowed %>%
  dplyr::filter(trait %in% snps_per_gene$trait) %>%
  dplyr::group_nest(trait) %>%
  mutate(hyprcoloc_res = map2(data, trait, function(df, trait_id) {
    
    message("Running HyPrColoc for trait: ", trait_id)
    locus_id <- dplyr::distinct(df, locus_marker) %>% dplyr::pull()
    
    df$trait <- trait_id
    
    df <- df %>%
      dplyr::bind_rows(sumstats_windowed %>% dplyr::filter(locus_marker == locus_id)) %>%
      dplyr::group_by(rsID, trait) %>%
      dplyr::filter(n() == 1) %>%  # keep only rows with unique rsID-trait combos
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(rsID) & rsID != "" & rsID %in% valid_rsid)
    
    .betas <- df %>%
      select(rsID, BETA, trait) %>%
      dplyr::distinct(rsID, BETA, trait) %>%
      pivot_wider(names_from = trait, values_from = BETA) %>%
      column_to_rownames(var = "rsID") %>%
      tidyr::drop_na() %>%
      as.matrix()
    
    .ses <- df %>%
      dplyr::select(rsID, SE, trait) %>%
      dplyr::distinct(rsID, SE, trait) %>%
      pivot_wider(names_from = trait, values_from = SE) %>%
      column_to_rownames(var = "rsID") %>%
      tidyr::drop_na() %>%
      as.matrix()
    
    .complete_idx <- complete.cases(.betas)
    
    .betas <- .betas[.complete_idx, , drop = FALSE]
    .ses <- .ses[.complete_idx, , drop = FALSE]
    
    .ids <- rownames(.betas)
    .trait_names <- colnames(.betas)
    
    .result <- hyprcoloc::hyprcoloc(.betas, .ses, trait.names = .trait_names, snp.id = .ids)
    return(.result[[1]])
    
  }
  ))

# Extract results
hyprcoloc_res_all <- hyprcoloc_res %>%
  dplyr::select(-data) %>%
  tidyr::unnest(hyprcoloc_res) %>%  
  dplyr::filter(traits != "None")

# Save results
fwrite(x = hyprcoloc_res_all, 
       file = snakemake@output[["eqtl"]],
       sep = "\t")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Install hyprcoloc and dependencies if needed
if (!requireNamespace("hyprcoloc", quietly = TRUE)) {
  # remotes::install_version("RcppEigen", version = "0.3.3.9.3")
  remotes::install_github(
    "jrs95/hyprcoloc",
    build_vignettes = TRUE,
    dependencies = FALSE
  )
}

# Load packages
library(hyprcoloc)
library(tidygenomics)
library(data.table)
library(tidyverse)

# Load custom scripts
creds <- snakemake@params[["creds"]]
locus <- snakemake@params[["locus"]]
source(creds)
source(locus)

### LOAD SNAKEMAKE OBJECTS
phenotypes <- snakemake@params[["phenotypes"]]
sumstats_files <- snakemake@input[["sumstats"]]
lead_files <- snakemake@input[["lead"]]

### LOAD SUMSTATS
sumstats <- purrr::map2(phenotypes, sumstats_files, function(phenotype, file) {
  data.table::fread(file) %>%
    dplyr::mutate(trait = phenotype) %>%
    dplyr::select(trait, rsID, CHR, POS, BETA, SE, EA, NEA, P)
}) %>%
  dplyr::bind_rows()

### LOAD LEAD VARIANTS
lead_vars <- purrr::map2(phenotypes, lead_files, function(phenotype, file) {
  vroom::vroom(file) %>%
    dplyr::mutate(phenotype = phenotype)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(START = POS - 500000, END = POS + 500000) %>%
  tidygenomics::genome_cluster(
    by = c("CHR", "START", "END"),
    max_distance = 0
  ) %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::summarise(
    CHR = first(CHR),
    START = min(START),
    END = max(END),
    ID = glue::glue("{CHR}:{START}:{END}"),
    NTRAITS = n_distinct(phenotype),
    TRAITS = list(unique(phenotype))
  ) %>%
  dplyr::filter(NTRAITS > 1)

sumstats_windowed <- purrr::map(1:nrow(lead_vars), function(i) {
  lead_vars <- lead_vars[i, ]
  traits_to_include <- unlist(lead_vars$TRAITS)

  df <- sumstats %>%
    dplyr::filter(trait %in% traits_to_include) %>%
    dplyr::filter(CHR == lead_vars$CHR) %>%
    dplyr::filter(dplyr::between(POS, lead_vars$START, lead_vars$END)) %>%
    dplyr::mutate(locus_marker = lead_vars$ID)

  return(df)
}) %>%
  dplyr::bind_rows()


# Run Coloc
hyprcoloc_res <- sumstats_windowed %>%
  dplyr::group_nest(locus_marker) %>%
  mutate(
    hyprcoloc_res = map2(data, locus_marker, function(df, locus_id) {
      message("Running HyPrColoc for locus: ", locus_id)

      df <- df %>%
        dplyr::group_by(rsID, trait) %>%
        dplyr::filter(n() == 1) %>% # keep only rows with unique rsID-trait combos
        dplyr::ungroup()

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

      .ids <- rownames(.betas)

      .trait_names <- colnames(.betas)

      .result <- hyprcoloc::hyprcoloc(
        .betas,
        .ses,
        trait.names = .trait_names,
        snp.id = .ids,
        snpscores = TRUE
      )
      .credset <- cred.sets(.result)
      return(list(coloc = .result[[1]], credset = .credset))
    })
  )


# Extract results
hyprcoloc_res_all <- hyprcoloc_res %>%
  dplyr::select(-data) %>%
  tidyr::unnest_wider(hyprcoloc_res) %>%
  tidyr::unnest(coloc) %>%
  dplyr::filter(traits != "None") %>%
  dplyr::select(-credset)

# Extract credible sets
credset <- hyprcoloc_res %>%
  dplyr::select(-data) %>%
  tidyr::unnest_wider(hyprcoloc_res) %>%
  tidyr::unnest(coloc) %>%
  dplyr::filter(traits != "None") %>%
  tidyr::unnest_longer(credset, indices_include = TRUE) %>%
  mutate(
    credset = map(credset, ~ tibble(rsID = names(.x), PIP = as.numeric(.x)))
  ) %>%
  unnest(credset) %>%
  dplyr::select(locus_marker, traits, rsID, PIP, credset_id) %>%
  dplyr::inner_join(
    sumstats_windowed %>%
      dplyr::select(
        rsID,
        CHR,
        POS,
        EA,
        NEA
      ) %>%
      dplyr::distinct(rsID, .keep_all = TRUE),
    by = "rsID"
  ) %>%
  dplyr::rename(
    CS = credset_id,
    TRAITS = traits,
    LOCUS_MARKER = locus_marker
  ) %>%
  dplyr::select(LOCUS_MARKER, TRAITS, rsID, CHR, POS, EA, NEA, PIP, CS)

# Save results
fwrite(
  x = hyprcoloc_res_all,
  file = snakemake@output[["colocalization"]],
  sep = "\t"
)

fwrite(x = credset, file = snakemake@output[["credset"]], sep = "\t")

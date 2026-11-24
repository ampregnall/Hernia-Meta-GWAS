library(dplyr)

# Identify lead variants from the AFR and AMR meta-analyses and cross-reference them
# against the lead variants from the multi-population (ALL) and EUR meta-analyses.
# NB: variants are matched on exact SNPID, not on overlapping loci.
#
# SHARED column:
#   UNIQUE    - not a lead variant in the ALL or EUR meta-analyses
#   ALL       - lead variant in the ALL but not the EUR meta-analysis
#   EUR       - lead variant in the EUR but not the ALL meta-analysis
#   EUR & ALL - lead variant in both the ALL and EUR meta-analyses

# Femoral hernia summary statistics are only available for EUR populations
phenotypes <- c("Ventral", "Inguinal", "Umbilical", "Diaphragmatic")
lead_dir <- "data/processed/lead-variants/meta"
out_file <- "data/processed/lead-variants/population-specific/population_specific_lead_variants.txt"

read_lead_variants <- function(phenotype, population) {
  readr::read_delim(
    glue::glue("{lead_dir}/{phenotype}-{population}_lead_variants.txt"),
    col_types = "cnnccnnnnnnnnnc"
  ) %>%
    dplyr::mutate(PHENOTYPE = phenotype, POP = population, .before = 1)
}

process_lead_variants <- function(phenotype) {
  all <- read_lead_variants(phenotype, "ALL")
  eur <- read_lead_variants(phenotype, "EUR")
  noneur <- dplyr::bind_rows(
    read_lead_variants(phenotype, "AMR"),
    read_lead_variants(phenotype, "AFR")
  )

  noneur %>%
    dplyr::mutate(
      SHARED = dplyr::case_when(
        SNPID %in% eur$SNPID & SNPID %in% all$SNPID ~ "EUR & ALL",
        SNPID %in% eur$SNPID ~ "EUR",
        SNPID %in% all$SNPID ~ "ALL",
        TRUE ~ "UNIQUE"
      ),
      .after = POP
    )
}

population_specific <- purrr::map(phenotypes, process_lead_variants) %>%
  dplyr::bind_rows()

fs::dir_create(dirname(out_file))
data.table::fwrite(population_specific, out_file, sep = "\t")

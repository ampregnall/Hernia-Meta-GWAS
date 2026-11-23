library(data.table)
library(tidyverse)

create_gene_prioritization_bf_df <- function(phenotype) {

  df <- fread(glue::glue("data/processed/finemapping/{phenotype}-ALL_bayesian_finemapping.txt"))
  gene_map <- fread("data/meta/pops/utils/gene_annot.txt")

  eqtl <- fread(glue::glue("data/processed/variant-annotation/eqtl/{phenotype}_eqtl_colocalization.txt")) %>%
    dplyr::mutate(ENSGID = stringr::str_split_i(trait, " ", 1)) %>%
    dplyr::left_join(., gene_map, by = "ENSGID") %>%
    dplyr::filter(!is.na(NAME)) %>%
    dplyr::select(rsID = candidate_snp,
                  EQTL_GENE = NAME,
                  EQTL_PIP = posterior_prob,
                  EQTL_PIP_SNP = posterior_explained_by_snp)

  if (nrow(eqtl) > 0) {
    df <- df %>%
      dplyr::left_join(eqtl, by = "rsID")
  } else {
    # Add empty columns to maintain consistent structure
    df <- df %>%
      dplyr::mutate(
        EQTL_GENE = NA_character_,
        EQTL_PIP = NA_real_,
        EQTL_PIP_SNP = NA_real_
      )
  }

  # Pops
  pops <- fread(glue::glue("data/processed/pops/{phenotype}.preds")) %>%
    dplyr::left_join(gene_map, by = "ENSGID")

  df <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(UPPER = POS + 5e5, LOWER = max(0, POS - 5e5))

  by <- dplyr::join_by(CHR, overlaps(x$LOWER, x$UPPER, y$START, y$END))

  df <- dplyr::full_join(df, pops, by) %>%
    dplyr::group_by(rsID) %>%
    dplyr::arrange(desc(PoPS_Score)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(rsID)) %>%
    dplyr::select(SNPID, rsID, LOCUS, CHR, POS, EA, NEA, BETA, SE, NEAREST_GENE, BF,
                  BF_PIP, EQTL_GENE, EQTL_PIP, EQTL_PIP_SNP,
                  POPS_SCORE = PoPS_Score, POPS_GENE = NAME)

  ms <- fread(glue::glue("data/processed/vep/output/{phenotype}-ALL.txt"), skip = "#Uploaded_variation") %>%
    dplyr::select(SNPID = `#Uploaded_variation`, MS_CONSEQUENCE = Consequence, ENSGID = Gene) %>%
    dplyr::left_join(gene_map, by = "ENSGID") %>%
    dplyr::select(!ENSGID)

  df <- df %>%
    dplyr::left_join(ms, by="SNPID")

  return(df)

}

ventral <- create_gene_prioritization_bf_df(phenotype = "Ventral")
inguinal <- create_gene_prioritization_bf_df(phenotype = "Inguinal")
umbilical <- create_gene_prioritization_bf_df(phenotype = "Umbilical")
femoral <- create_gene_prioritization_bf_df(phenotype = "Femoral")
diaphragmatic <- create_gene_prioritization_bf_df(phenotype = "Diaphragmatic")

fwrite(ventral, "data/processed/variant-annotation/credsets-annotated/Ventral_bayesian_annotated.txt", sep = "\t")
fwrite(inguinal, "data/processed/variant-annotation/credsets-annotated/Inguinal_bayesian_annotated.txt", sep = "\t")
fwrite(umbilical, "data/processed/variant-annotation/credsets-annotated/Umbilical_bayesian_annotated.txt", sep = "\t")
fwrite(femoral, "data/processed/variant-annotation/credsets-annotated/Femoral_bayesian_annotated.txt", sep = "\t")
fwrite(diaphragmatic, "data/processed/variant-annotation/credsets-annotated/Diaphragmatic_bayesian_annotated.txt", sep = "\t")

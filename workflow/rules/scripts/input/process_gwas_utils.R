# Load libraries
library(magrittr)

# Autosomes plus X. Y and MT are not analyzed.
VALID_CHROMS <- c(as.character(1:22), "X")

# Standardize chromosome labels across studies to 1-22 and "X". Studies code X as
# "X", "chrX", or 23; numeric conversion would turn X into NA and drop it.
standardise_chrom <- function(chrom) {
  chrom <- toupper(sub("^chr", "", as.character(chrom), ignore.case = TRUE))
  chrom[chrom == "23"] <- "X"
  chrom
}

liftover_37_to_38 <- function(df, chr, pos, chainfile) {
  df <- df %>%
    dplyr::mutate("{{chr}}" := paste0("chr", {{ chr }})) %>%
    dplyr::mutate(start = {{ pos }}, end = start) %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      starts.in.df.are.0based = FALSE
    ) %>%
    rtracklayer::liftOver(., rtracklayer::import.chain(chainfile)) %>%
    tidyr::as_tibble()

  return(df)
}

process_bjp <- function(df, n_case, n_control) {
  df <- df %>%
    dplyr::select(!POS) %>% # Remove original POS column
    dplyr::mutate(
      CHROM = standardise_chrom(seqnames)
    ) %>%
    dplyr::mutate(N_CASE = n_case, N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::rename(
      POS = start,
      P = p.value,
      EAF = AF.Cases,
      EFFECT_ALLELE = Allele2,
      OTHER_ALLELE = Allele1
    ) %>%
    dplyr::select(
      CHROM,
      POS,
      EFFECT_ALLELE,
      OTHER_ALLELE,
      EAF,
      BETA,
      SE,
      P,
      N_CASE,
      N_CONTROL,
      N
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(CHROM = standardise_chrom(CHROM)) %>%
    dplyr::filter(CHROM %in% VALID_CHROMS) %>%
    dplyr::mutate(
      MARKER = glue::glue("{CHROM}:{POS}:{OTHER_ALLELE}:{EFFECT_ALLELE}")
    ) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()

  return(df)
}

process_ukbb <- function(df, n_case, n_control) {
  df <- df %>%
    dplyr::mutate(
      CHROM = standardise_chrom(seqnames)
    ) %>%
    dplyr::mutate(N_CASE = n_case, N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::rename(
      POS = start,
      BETA = dplyr::starts_with("beta"),
      SE = dplyr::starts_with("se_"),
      neglog10pval = dplyr::starts_with("neglog10"),
      EAF = dplyr::starts_with("af_cases"),
      EFFECT_ALLELE = alt,
      OTHER_ALLELE = ref
    ) %>%
    dplyr::mutate(P = 10^(-1 * neglog10pval)) %>%
    dplyr::select(
      CHROM,
      POS,
      EFFECT_ALLELE,
      OTHER_ALLELE,
      EAF,
      BETA,
      SE,
      P,
      N_CASE,
      N_CONTROL,
      N
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(CHROM = standardise_chrom(CHROM)) %>%
    dplyr::filter(CHROM %in% VALID_CHROMS) %>%
    dplyr::mutate(
      MARKER = glue::glue("{CHROM}:{POS}:{OTHER_ALLELE}:{EFFECT_ALLELE}")
    ) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()

  return(df)
}

process_mvp <- function(df, n_case, n_control) {
  df <- df %>%
    dplyr::mutate(N_CASE = n_case, N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::rename(
      CHROM = CHR,
      POS = POS,
      EFFECT_ALLELE = ALT,
      OTHER_ALLELE = REF,
      EAF = AAF
    ) %>%
    dplyr::select(
      CHROM,
      POS,
      EFFECT_ALLELE,
      OTHER_ALLELE,
      EAF,
      BETA,
      SE,
      P,
      N_CASE,
      N_CONTROL,
      N
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(CHROM = standardise_chrom(CHROM)) %>%
    dplyr::filter(CHROM %in% VALID_CHROMS) %>%
    dplyr::mutate(
      MARKER = glue::glue("{CHROM}:{POS}:{OTHER_ALLELE}:{EFFECT_ALLELE}")
    ) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()

  return(df)
}

process_finngen <- function(df, n_case, n_control) {
  df <- df %>%
    dplyr::mutate(N_CASE = n_case, N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::select(
      `#chrom`,
      pos,
      alt,
      ref,
      af_alt,
      beta,
      sebeta,
      pval,
      N_CASE,
      N_CONTROL,
      N
    ) %>%
    tidyr::drop_na() %>%
    dplyr::rename(
      CHROM = `#chrom`,
      POS = pos,
      EFFECT_ALLELE = alt,
      OTHER_ALLELE = ref,
      EAF = af_alt,
      BETA = beta,
      SE = sebeta,
      P = pval
    ) %>%
    dplyr::mutate(CHROM = standardise_chrom(CHROM)) %>%
    dplyr::filter(CHROM %in% VALID_CHROMS) %>%
    dplyr::mutate(
      MARKER = glue::glue("{CHROM}:{POS}:{OTHER_ALLELE}:{EFFECT_ALLELE}")
    ) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()

  return(df)
}

process_aou <- function(df, n_case, n_control) {
  df <- df %>%
    dplyr::mutate(N_CASE = n_case, N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::mutate(CHROM = standardise_chrom(CHR)) %>%
    dplyr::rename(
      EFFECT_ALLELE = Allele2,
      OTHER_ALLELE = Allele1,
      EAF = AF_Allele2,
      P = Pvalue
    ) %>%
    dplyr::select(
      CHROM,
      POS,
      EFFECT_ALLELE,
      OTHER_ALLELE,
      EAF,
      BETA,
      SE,
      P,
      N_CASE,
      N_CONTROL,
      N
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(CHROM = standardise_chrom(CHROM)) %>%
    dplyr::filter(CHROM %in% VALID_CHROMS) %>%
    dplyr::mutate(
      MARKER = glue::glue("{CHROM}:{POS}:{OTHER_ALLELE}:{EFFECT_ALLELE}")
    ) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()

  return(df)
}

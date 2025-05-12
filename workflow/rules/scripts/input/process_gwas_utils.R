# Load libraries
library(magrittr)

liftover_37_to_38 <- function(df, chainfile) {
  
  df <- df %>%
    dplyr::mutate(chr = glue::glue("chr{chr}")) %>%
    dplyr::mutate(start = pos, end = start) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = FALSE) %>%
    rtracklayer::liftOver(., rtracklayer::import.chain(chainfile)) %>%
    tidyr::as_tibble()
  
  return(df)
  
}

process_ukbb <- function(df, n_case, n_control) {
  
  df <- df %>%
    dplyr::mutate(CHROM = as.numeric(stringr::str_replace(seqnames, "chr", ""))) %>%
    dplyr::mutate(N_CASE = n_case, N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::rename(POS = pos, BETA = dplyr::starts_with("beta"), 
                  SE = dplyr::starts_with("se_"), 
                  neglog10pval = dplyr::starts_with("neglog10"),
                  af_cases = dplyr::starts_with("af_cases"),
                  af_controls = dplyr::starts_with("af_controls"),
                  EFFECT_ALLELE = alt, OTHER_ALLELE = ref) %>%
    dplyr::mutate(EAF = (N_CASE*af_cases + N_CONTROL*af_controls)/(N)) %>%
    dplyr::mutate(P = 10^(-1*neglog10pval)) %>%
    dplyr::select(CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, EAF, BETA, SE, P, N_CASE, N_CONTROL, N) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(MARKER = glue::glue("{CHROM}_{POS}")) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()
  
    return(df)
  
}

process_mvp <- function(df, n_case, n_control) {
  
  df <- df %>% 
    dplyr::mutate(BETA = log(or), SE = abs(BETA / qnorm(pval / 2))) %>%
    dplyr::mutate(N_CASE = n_case,N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::rename(CHROM = chrom, POS = pos, EFFECT_ALLELE = ea, 
                  OTHER_ALLELE = ref, EAF = af, P = pval) %>%
    dplyr::select(CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, 
                  EAF, BETA, SE, P, N_CASE, N_CONTROL, N) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(MARKER = glue::glue("{CHROM}_{POS}")) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()
  
  return(df)
      
}

process_finngen <- function(df, n_case, n_control) {
  df <- df %>%
    dplyr::mutate(N_CASE = n_case, N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::select(`#chrom`, pos, alt, ref, af_alt, beta, sebeta, pval, N_CASE, N_CONTROL, N) %>%
    tidyr::drop_na() %>%
    dplyr::rename(CHROM = `#chrom`, POS = pos, EFFECT_ALLELE = alt, OTHER_ALLELE = ref, 
                  EAF = af_alt, BETA = beta, SE = sebeta, P = pval) %>%
    dplyr::mutate(MARKER = glue::glue("{CHROM}_{POS}")) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()
  
  return(df)
  
}

process_aou <- function(df, n_case, n_control) {
  
  df <- df %>%
    dplyr::mutate(N_CASE = n_case, N_CONTROL = n_control) %>%
    dplyr::mutate(N = N_CASE + N_CONTROL) %>%
    dplyr::mutate(CHROM = as.numeric(stringr::str_replace(CHR, "chr", ""))) %>%
    dplyr::rename(EFFECT_ALLELE = Allele2, OTHER_ALLELE = Allele1, EAF = AF_Allele2, P = Pvalue) %>%
    dplyr::select(CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, EAF, BETA, SE, P, N_CASE, N_CONTROL, N) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(MARKER = glue::glue("{CHROM}_{POS}")) %>%
    dplyr::select(MARKER, everything()) %>%
    tidyr::as_tibble()
  
  return(df)

}

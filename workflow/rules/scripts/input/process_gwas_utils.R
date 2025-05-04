# Load libraries
library(data.table)
library(dtplyr)
library(magrittr)
library(tidyr)
library(dplyr)

process_ukbb <- function(df) {
  df <- df %>%
    rename(beta = 7, se = 8, neglog10pval = 9) %>%
    mutate(pvalue = 10^(-1*neglog10pval)) %>%
    select(chr, pos, ref, alt, beta, se, pvalue) %>%
    drop_na() %>%
    as_tibble()
  
  return(df)
  
}

process_mvp <- function(df) {
  df <- df %>% 
    separate_wider_delim(cols = ci, delim = ",", names = c("ci_low", "ci_high")) %>%
    mutate(beta = log(or), se = (as.numeric(ci_high) - or) / 1.96) %>%
    select(chrom, pos, ref, alt, beta, se, pval) %>%
    drop_na() %>%
    rename(chr = chrom, pvalue = pval) %>%
    as_tibble()
  
  return(df)
      
}

process_finngen <- function(df) {
  df <- df %>%
    select(`#chrom`, pos, ref, alt, beta, sebeta, pval) %>%
    drop_na() %>%
    rename(chr = `#chrom`, se = sebeta, pvalue = pval) %>%
    as_tibble()
  
  return(df)
  
}

process_aou <- function(df) {
  
  df <- df %>%
    select(CHR, POS, Allele1, Allele2, BETA, SE, Pvalue) %>%
    drop_na() %>%
    rename(chr = CHR, pos = POS, ref = Allele1, alt = Allele2, beta = BETA, se = SE, pvalue = Pvalue) %>%
    as_tibble()
  
  return(df)

}
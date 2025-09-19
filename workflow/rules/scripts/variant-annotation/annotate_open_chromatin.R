library(GenomicRanges)
library(data.table)
library(rtracklayer)
library(dplyr)

### Load chromatin bed files
h3k4me1_p <- rtracklayer::import("data/meta/open-chromatin/Penis_Foreskin_Fibroblast_Primary_Cells.H3K4me1.GrCh38.bed")
h3k4me3_p <- rtracklayer::import("data/meta/open-chromatin/Penis_Foreskin_Fibroblast_Primary_Cells.H3K4me3.GrCh38.bed")
h3k27me3_p <- rtracklayer::import("data/meta/open-chromatin/Penis_Foreskin_Fibroblast_Primary_Cells.H3K27me3.GrCh38.bed")

h3k4me1_b <- rtracklayer::import("data/meta/open-chromatin/Breast_Fibroblast_Primary_Cells.H3K4me1.GrCh38.bed")
h3k4me3_b <- rtracklayer::import("data/meta/open-chromatin/Breast_Fibroblast_Primary_Cells.H3K27me3.GrCh38.bed")
h3k27me3_b <- rtracklayer::import("data/meta/open-chromatin/Breast_Fibroblast_Primary_Cells.H3K27me3.GrCh38.bed")

### Load finemapping results 
snps <- fread(snakemake@input[["snps_in"]])

annotate_open_chromatin <- function(df) {
  snps_gr <- GRanges(
    seqnames = ifelse(grepl("^chr", df$CHR), 
                      df$CHR, 
                      paste0("chr", df$CHR)),
    ranges = IRanges(
      start = df$POS,  # GRanges is 1-based, so no need to subtract 1
      end = df$POS     # SNPs are single positions
    ),
    rsID = df$rsID
  )
  
  h3k4me1_p_overlaps <- findOverlaps(snps_gr, h3k4me1_p)
  h3k4me3_p_overlaps <- findOverlaps(snps_gr, h3k4me3_p)
  h3k27me3_p_overlaps <- findOverlaps(snps_gr, h3k27me3_p)
  
  snps_in_h3k4me1_p <- snps_gr[queryHits(h3k4me1_p_overlaps)] %>% tibble::as_tibble()
  snps_in_h3k4me3_p <- snps_gr[queryHits(h3k4me3_p_overlaps)] %>% tibble::as_tibble()
  snps_in_h3k27me3_p <- snps_gr[queryHits(h3k27me3_p_overlaps)] %>% tibble::as_tibble()
  
  h3k4me1_b_overlaps <- findOverlaps(snps_gr, h3k4me1_b)
  h3k4me3_b_overlaps <- findOverlaps(snps_gr, h3k4me3_b)
  h3k27me3_b_overlaps <- findOverlaps(snps_gr, h3k27me3_b)
  
  snps_in_h3k4me1_b <- snps_gr[queryHits(h3k4me1_b_overlaps)] %>% tibble::as_tibble()
  snps_in_h3k4me3_b <- snps_gr[queryHits(h3k4me3_b_overlaps)] %>% tibble::as_tibble()
  snps_in_h3k27me3_b <- snps_gr[queryHits(h3k27me3_b_overlaps)] %>% tibble::as_tibble()
  
  df <- df %>%
    dplyr::mutate(inH3K4ME1_penis = case_when(rsID %in% snps_in_h3k4me1_p$rsID ~ 1, TRUE ~ 0)) %>%
    dplyr::mutate(inH3K4ME3_penis = case_when(rsID %in% snps_in_h3k4me3_p$rsID ~ 1, TRUE ~ 0)) %>%
    dplyr::mutate(inH3K27ME3_penis = case_when(rsID %in% snps_in_h3k27me3_p$rsID ~ 1, TRUE ~ 0))  %>%
    dplyr::mutate(inH3K4ME1_breast = case_when(rsID %in% snps_in_h3k4me1_b$rsID ~ 1, TRUE ~ 0)) %>%
    dplyr::mutate(inH3K4ME3_breast = case_when(rsID %in% snps_in_h3k4me3_b$rsID ~ 1, TRUE ~ 0)) %>%
    dplyr::mutate(inH3K27ME3_breast = case_when(rsID %in% snps_in_h3k27me3_b$rsID ~ 1, TRUE ~ 0)) %>%
    dplyr::select(rsID, inH3K4ME1_penis, inH3K4ME3_penis, inH3K27ME3_penis, 
                  inH3K4ME1_breast, inH3K4ME3_breast, inH3K27ME3_breast)
  
  return(df)
    
                    
}

snps <- annotate_open_chromatin(snps)

fwrite(snps, snakemake@output[["snps_out"]], sep = "\t")
suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(glue)
  library(purrr)
  library(hyprcoloc)
})

# ---- MR hits (fast + minimal passes) ----
mr_files <- c(
  "data/processed/mendelian-randomization/archive/Inguinal_mr.txt",
  "data/processed/mendelian-randomization/archive/Ventral_mr.txt",
  "data/processed/mendelian-randomization/archive/Umbilical_mr.txt",
  "data/processed/mendelian-randomization/archive/Diaphragmatic_mr.txt",
  "data/processed/mendelian-randomization/archive/Femoral_mr.txt"
)

mr_hits <- vroom::vroom(mr_files, id = "Phenotype", show_col_types = FALSE) |>
  mutate(
    Phenotype = str_split_i(basename(Phenotype), "_", 1)
  ) |>
  filter(
    str_detect(exposure, "Fenland"),
    (nsnp > 1 & method == "Inverse variance weighted") |
      (nsnp == 1 & method == "MR Egger")
  ) |>
  group_by(Phenotype) |>
  mutate(fdr = p.adjust(pval, method = "fdr")) |>
  ungroup() |>
  filter(fdr < 0.05) |>
  mutate(Gene = str_split_i(exposure, ":", 4)) |>
  distinct(Phenotype, Gene)   # only need unique pairs for coloc

# ---- pQTL file map ----
pqtl_files <- list(
  "VPS26A" = "data/meta/pqtl/VP26A_19293_6_GRCh38.tsv.gz",
  "H6PD" = "data/meta/pqtl/G6PE_7161_25_GRCh38.tsv.gz",
  "EFEMP1" = "data/meta/pqtl/FBLN3_8480_29_GRCh38.tsv.gz",
  "ABO" = "data/meta/pqtl/BGAT_9253_52_GRCh38.tsv.gz",
  "LGALS9" = "data/meta/pqtl/LEG9_9197_4_GRCh38.tsv.gz",
  "LILRB3" = "data/meta/pqtl/LIRB3_11334_7_GRCh38.tsv.gz",
  "GPX1" = "data/meta/pqtl/Glutathione_peroxidase_15591_28_GRCh38.tsv.gz",
  "LGALS3" = "data/meta/pqtl/Galectin_3_3066_12_GRCh38.tsv.gz",
  "ITIH1" = "data/meta/pqtl/ITI_heavy_chain_H1_7955_195_GRCh38.tsv.gz",
  "MAPK3" = "data/meta/pqtl/ERK_1_2855_49_GRCh38.tsv.gz",
  "YEATS4" = "data/meta/pqtl/YETS4_17474_106_GRCh38.tsv.gz",
  "CFHR1" = "data/meta/pqtl/FHR1_15468_14_GRCh38.tsv.gz",
  "NAE1" = "data/meta/pqtl/ULA1_18177_49_GRCh38.tsv.gz",
  "ERAP1" = "data/meta/pqtl/ARTS1_4964_67_GRCh38.tsv.gz",
  "BTN3A3" = "data/meta/pqtl/BT3A3_17692_2_GRCh38.tsv.gz",
  "B3GLCT" = "data/meta/pqtl/B3GLT_5727_35_GRCh38.tsv.gz",
  "CSGALNACT2" = "data/meta/pqtl/CGAT2_10772_21_GRCh38.tsv.gz",
  "TXNDC5" = "data/meta/pqtl/TXND5_11212_7_GRCh38.tsv.gz",
  "TAPBP" = "data/meta/pqtl/TPSN_12378_71_GRCh38.tsv.gz",
  "DNAJB4" = "data/meta/pqtl/DNJB4_18884_22_GRCh38.tsv.gz",
  "CBL" = "data/meta/pqtl/CBL_12016_60_GRCh38.tsv.gz",
  "SDF2" = "data/meta/pqtl/SDF2_15333_11_GRCh38.tsv.gz",
  "PRTN3" = "data/meta/pqtl/Proteinase_3_13720_95_GRCh38.tsv.gz",
  "MFAP4" = "data/meta/pqtl/MFAP4_5636_10_GRCh38.tsv.gz",
  "MRC2" = "data/meta/pqtl/MRC2_3041_55_GRCh38.tsv.gz",
  "CTRB2" = "data/meta/pqtl/CTRB2_5648_28_GRCh38.tsv.gz",
  "VEGFC" = "data/meta/pqtl/VEGF_C_3132_1_GRCh38.tsv.gz",
  "CD209" = "data/meta/pqtl/DC_SIGN_3029_52_GRCh38.tsv.gz",
  "SHBG" = "data/meta/pqtl/SHBG_4929_55_GRCh38.tsv.gz",
  "LIMA1" = "data/meta/pqtl/LIMA1_11543_84_GRCh38.tsv.gz",
  "TGFBI" = "data/meta/pqtl/BGH3_3283_21_GRCh38.tsv.gz",
  "PDIA5" = "data/meta/pqtl/PDIA5_5593_11_GRCh38.tsv.gz",
  "CNTFR" = "data/meta/pqtl/CNTFR_alpha_14101_2_GRCh38.tsv.gz",
  "WWOX" = "data/meta/pqtl/WWOX_17724_3_GRCh38.tsv.gz",
  "APOC3" = "data/meta/pqtl/Apo_C_III_6461_54_GRCh38.tsv.gz",
  "IDUA" = "data/meta/pqtl/IDUA_3169_70_GRCh38.tsv.gz",
  "ISLR2" = "data/meta/pqtl/ISLR2_13124_20_GRCh38.tsv.gz",
  "ADAMTS5" = "data/meta/pqtl/ADAMTS_5_3168_8_GRCh38.tsv.gz",
  "LZIC" = "data/meta/pqtl/LZIC_19338_3_GRCh38.tsv.gz",
  "APOL1" = "data/meta/pqtl/Apo_L1_11510_31_GRCh38.tsv.gz",
  "APOA5" = "data/meta/pqtl/Apo_A_V_15363_32_GRCh38.tsv.gz",
  "OSMR" = "data/meta/pqtl/OSMR_10892_8_GRCh38.tsv.gz",
  "FBLN1" = "data/meta/pqtl/fibulin_1_6470_19_GRCh38.tsv.gz",
  "LY75" = "data/meta/pqtl/LY75_16620_26_GRCh38.tsv.gz",
  "TRIL" = "data/meta/pqtl/TRIL_6527_1_GRCh38.tsv.gz",
  "HHIP" = "data/meta/pqtl/HHIP_10833_64_GRCh38.tsv.gz",
  "ENPP2" = "data/meta/pqtl/ENPP2_16892_23_GRCh38.tsv.gz",
  "TNFRSF6B" = "data/meta/pqtl/DcR3_5070_76_GRCh38.tsv.gz",
  "TNXB" = "data/meta/pqtl/Tenascin_X_5698_60_GRCh38.tsv.gz"
)

# ---- Gene coordinates (join once) ----
gene_location <- data.table::fread("data/meta/pops/utils/gene_annot.txt",
                                   select = c("ENSGID","NAME","CHR","START","END","TSS"))

gene_location_manual <- data.frame(
  ENSGID = c("ENSG00000187676", "ENSG00000175164", "ENSG00000255690"),
  NAME   = c("B3GLCT", "ABO", "TRIL"),
  CHR    = c(13, 9, 7),
  START  = c(131199975, 133233278, 28992974),
  END    = c(131332276, 133276024, 28997946),
  TSS    = c(NA, NA, NA)
)

genes <- bind_rows(gene_location, gene_location_manual) |>
  distinct(NAME, .keep_all = TRUE)

pairs <- mr_hits |>
  inner_join(genes, by = c("Gene" = "NAME")) |>
  mutate(pqtl_path = unname(pqtl_files[Gene])) |>
  filter(!is.na(pqtl_path)) |>
  distinct(Phenotype, Gene, CHR, START, END, pqtl_path)

# ---- Fast matrix builder using data.table::dcast ----
.make_mats <- function(dt) {
  # dt columns: rsID, trait, BETA, SE
  dt <- dt[!is.na(rsID) & rsID != ""]
  dt <- unique(dt, by = c("rsID","trait"))
  
  betas <- dcast(dt, rsID ~ trait, value.var = "BETA")
  ses   <- dcast(dt, rsID ~ trait, value.var = "SE")
  
  rb <- betas$rsID
  rs <- ses$rsID
  betas$rsID <- NULL
  ses$rsID <- NULL
  
  betas <- as.matrix(betas)
  ses   <- as.matrix(ses)
  
  rownames(betas) <- rb
  rownames(ses)   <- rs
  
  # Keep only shared SNPs
  common <- intersect(rownames(betas), rownames(ses))
  if (length(common) == 0) return(NULL)
  
  betas <- betas[common, , drop = FALSE]
  ses   <- ses[common,   , drop = FALSE]
  
  # Remove rows with ANY NA in either matrix
  complete_rows <- complete.cases(betas) & complete.cases(ses)
  
  if (!any(complete_rows)) return(NULL)
  
  betas <- betas[complete_rows, , drop = FALSE]
  ses   <- ses[complete_rows,   , drop = FALSE]
  
  list(
    betas = betas,
    ses   = ses
  )
}

# ---- Coloc runner (NO break/next; returns NULL on skip) ----
pqtl_colocalization <- function(Phenotype, Gene, CHR, START, END, pqtl_path) {
  message("Running HyPrColoc for: ", Phenotype, " / ", Gene)
  
  sum_path <- glue("data/processed/metal-formatted/{Phenotype}-ALL_metal_formatted.txt")
  if (!file.exists(sum_path)) return(NULL)
  if (!file.exists(pqtl_path)) return(NULL)
  
  pad <- 250000L
  lo  <- START - pad
  hi  <- END + pad
  
  # Read minimal columns only
  sumstats <- data.table::fread(sum_path, select = c("CHR","POS","rsID","BETA","SE"))
  pqtl     <- data.table::fread(pqtl_path, select = c("CHR","BP","SNP","BETA","SE"))
  
  # Filter to locus + de-dup
  sum_locus <- sumstats[CHR == CHR & POS >= lo & POS <= hi]
  if (nrow(sum_locus) == 0) return(NULL)
  sum_locus[, trait := Phenotype]
  sum_locus <- unique(sum_locus, by = c("rsID","trait"))
  sum_locus <- sum_locus[, .(rsID, trait, BETA, SE)]
  
  pqtl_locus <- pqtl[CHR == CHR & BP >= lo & BP <= hi]
  if (nrow(pqtl_locus) == 0) return(NULL)
  pqtl_locus[, `:=`(rsID = SNP, trait = Gene)]
  pqtl_locus <- unique(pqtl_locus, by = c("rsID","trait"))
  pqtl_locus <- pqtl_locus[, .(rsID, trait, BETA, SE)]
  
  dt <- rbindlist(list(sum_locus, pqtl_locus), use.names = TRUE, fill = TRUE)
  mats <- .make_mats(dt)
  if (is.null(mats)) return(NULL)
  
  # Need >= 2 traits and some SNPs
  #if (ncol(mats$betas) < 2 || nrow(mats$betas) < 10) return(NULL)
  
  out <- hyprcoloc::hyprcoloc(
    mats$betas, mats$ses,
    trait.names = colnames(mats$betas),
    snp.id      = rownames(mats$betas)
  )
  message("HyPrColoc run succesfully for: ", Phenotype, " / ", Gene)
  
  # Attach identifiers so you can trace results
  out$results$Phenotype <- Phenotype
  out$results$Gene <- Gene
  out
}

# ---- Run across pairs safely ----
pqtl_colocalization_safe <- purrr::possibly(pqtl_colocalization, otherwise = NULL)

res <- purrr::pmap(pairs, pqtl_colocalization_safe)
results <- res |>
  purrr::compact() |>     # remove NULLs
  purrr::flatten() |>     # remove the extra list layer
  purrr::list_rbind() 


vroom::vroom_write(results, "data/processed/mendelian-randomization/MR_PQTL_colocalization.txt", delim = "\t")
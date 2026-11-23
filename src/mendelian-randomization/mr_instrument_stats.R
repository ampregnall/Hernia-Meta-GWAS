files <- list.files("data/meta/mendelian-randomization", full.names = TRUE, pattern = ".txt")

df <- vroom::vroom(files)

df_split <- dplyr::group_split(df, Phenotype)

instrument_stats <- purrr::map(df_split, function(data) {
  
  data$rsid <- data$SNP
  pheno <- data[1,]$Phenotype
  samplesize <- mean(data$samplesize)
  snps_pre = nrow(data)
  # data <- dplyr::rename(data, pvalue = pval)
  
  df_clump <- try(
    ieugwasr::ld_clump(
      dat       = data,
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile     = "/project/voltron/Resources/1000-Genomes-Phase3/MRC-IEU/EUR"
    ),
    silent = TRUE
  )
  
  # On error or no rows, return empty tibble with correct columns
  if (inherits(df_clump, "try-error") || is.null(df_clump) || nrow(df_clump) == 0L) {
    df_clump <- data
  } 
  
  snps_post <- nrow(df_clump)
  df_clump$R <- calculate_r2(df_clump$beta, df_clump$se, df_clump$eaf, samplesize)
  total_R <- sum(df_clump$R)
  f <- calculate_f_from_r2(total_R, samplesize, k = snps_post)
  
  return(
    data.frame(
      Phenotype = pheno,
      samplesize = samplesize,
      instrument_snps_preclump = snps_pre,
      instrument_snps_postclump = snps_post,
      f_stat = f
    )
  )
  
}) |> dplyr::bind_rows()


utils <- snakemake@params[["utils"]]
dir <- snakemake@params[["data_dir"]]
metal_dir <- snakemake@params[["metal"]]
mega_dir <- snakemake@params[["mrmega"]]

source(utils)
studies <- list.files(dir, full.names = TRUE)

phenotypes <- snakemake@config$create_meta_analysis_configs$phenotypes
populations <- snakemake@config$create_meta_analysis_configs$populations

for (phenotype in phenotypes) {
  files <- grep(phenotype, studies, value = TRUE)

  for (population in populations) {
    if (population == "ALL") {
      metal_config(
        config_name = paste(phenotype, "ALL", sep = "-"),
        output_dir = metal_dir,
        study_files = files,
        MARKERLABEL = "SNPID",
        CHROMOSOMELABEL = "CHR",
        POSITIONLABEL = "POS",
        EFFECT_ALLELE = "EA",
        OTHER_ALLELE = "NEA",
        EFFECTLABEL = "BETA",
        STDERR = "SE",
        FREQLABEL = "EAF",
        NCASE = "N_CASE",
        NCONTROL = "N_CONTROL",
        SAMPLESIZE = "N"
      )
      
      writeLines(files, paste0(mega_dir, "/", phenotype, "-ALL.in"))
    } else {
      study_files <- grep(population, files, value = TRUE)
      
      if (length(study_files == 0)) {
        next
      }

      metal_config(
        config_name = paste(phenotype, population, sep = "-"),
        output_dir = metal_dir,
        study_files = study_files,
        MARKERLABEL = "SNPID",
        CHROMOSOMELABEL = "CHR",
        POSITIONLABEL = "POS",
        EFFECT_ALLELE = "EA",
        OTHER_ALLELE = "NEA",
        EFFECTLABEL = "BETA",
        STDERR = "SE",
        FREQLABEL = "EAF",
        NCASE = "N_CASE",
        NCONTROL = "N_CONTROL",
        SAMPLESIZE = "N"
      )
      
      writeLines(study_files, paste0(mega_dir, "/", phenotype, "-", population, ".in"))
    }
  }
}

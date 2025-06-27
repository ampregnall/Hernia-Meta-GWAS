utils <- snakemake@params[["utils"]]
dir <- snakemake@params[["data_dir"]]
metal_dir <- snakemake@params[["metal"]]
mega_dir <- snakemake@params[["mrmega"]]

source(utils)
studies <- list.files(dir, full.names = TRUE)

phenotypes <- snakemake@config$create_metal_config$phenotypes
populations <- snakemake@config$create_metal_config$populations

for (phenotype in phenotypes) {
  files <- grep(phenotype, studies, value = TRUE)

  for (population in populations) {
    if (population == "ALL") {
      metal_config(
        config_name = paste(phenotype, "ALL", sep = "-"),
        output_dir = metal_dir,
        study_files = files,
        MARKERLABEL = "MARKER",
        CHROMOSOMELABEL = "CHR",
        POSITIONLABEL = "POSITION",
        EFFECT_ALLELE = "EFFECT_ALL",
        OTHER_ALLELE = "OTHER_ALL",
        EFFECTLABEL = "EFFECT",
        STDERR = "STDERR",
        FREQLABEL = "EFF_ALL_FREQ",
        NCASE = "N_CASE",
        NCONTROL = "N_CONTROL",
        SAMPLESIZE = "N_TOTAL"
      )
      
      writeLines(files, paste0(mega_dir, "/", phenotype, "-ALL.in"))
    } else {
      study_files <- grep(population, files, value = TRUE)

      metal_config(
        config_name = paste(phenotype, population, sep = "-"),
        output_dir = metal_dir,
        study_files = study_files,
        MARKERLABEL = "MARKER",
        CHROMOSOMELABEL = "CHR",
        POSITIONLABEL = "POSITION",
        EFFECT_ALLELE = "EFFECT_ALL",
        OTHER_ALLELE = "OTHER_ALL",
        EFFECTLABEL = "EFFECT",
        STDERR = "STDERR",
        FREQLABEL = "EFF_ALL_FREQ",
        NCASE = "N_CASE",
        NCONTROL = "N_CONTROL",
        SAMPLESIZE = "N_TOTAL"
      )
      
    }
  }
}

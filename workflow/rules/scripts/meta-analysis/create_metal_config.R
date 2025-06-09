source("workflow/rules/scripts/metal/metal_config.R")
studies <- list.files("data/processed/gwas-inspector-sumstats", full.names = TRUE)

phenotypes <- c("Abdominal", "Femoral", "Diaphragmatic", "Ventral", "Umbilical", "Inguinal")
ancestries <- c("EUR", "AMR", "EAS", "CSA", "AFR", "MID")

for (phenotype in phenotypes) {
  files <- grep(phenotype, studies, value = TRUE)
  
  metal_config(
    config_name = paste(phenotype, "ALL", sep="-"),
    output_dir = "data/processed/metal",
    study_files = files,
    MARKERLABEL = "varID",
    CHROMOSOMELABEL = "CHR",
    POSITIONLABEL = "POSITION",
    EFFECT_ALLELE = "EFFECT_ALL",
    OTHER_ALLELE = "OTHER_ALL",
    EFFECTLABEL = "EFFECT",
    STDERR = "STDERR",
    FREQLABEL = "EFF_ALL_FREQ",
    NCASE = "N_CASE",
    NCONTROL = "N_CONTROL",
    SAMPLESIZE = "N_TOTAL")
  
    for (ancestry in ancestries) {
    study_files <- grep(ancestry, files, value = TRUE)
    
    if (length(study_files) > 1) {
      metal_config(
        config_name = paste(phenotype, ancestry, sep="-"),
        output_dir = "data/processed/metal",
        study_files = study_files,
        MARKERLABEL = "varID",
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

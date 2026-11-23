utils <- snakemake@params[["utils"]]
dir <- snakemake@params[["data_dir"]]
metal_dir <- snakemake@params[["metal"]]

source(utils)
studies <- list.files(dir, full.names = TRUE)

phenotypes <- snakemake@config$create_meta_analysis_configs$phenotypes
populations <- snakemake@config$create_meta_analysis_configs$populations

write_metal_config <- function(config_name, study_files) {
  metal_config(
    config_name = config_name,
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
}

for (phenotype in phenotypes) {
  files <- grep(phenotype, studies, value = TRUE)

  for (population in populations) {
    if (population == "ALL") {
      write_metal_config(paste(phenotype, "ALL", sep = "-"), files)
    } else {
      study_files <- grep(population, files, value = TRUE)

      # Population-specific meta-analyses require >= 2 sumstats. Single-study EAS
      # sumstats are also passed through METAL so they share the same output format.
      # Must match the valid_pairs / metal_pairs logic in workflow/Snakefile.
      if (length(study_files) >= 2 || (population == "EAS" && length(study_files) == 1)) {
        write_metal_config(paste(phenotype, population, sep = "-"), study_files)
      }
    }
  }
}

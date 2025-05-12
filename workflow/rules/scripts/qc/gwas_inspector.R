# if (!requireNamespace("GWASinspector", quietly = TRUE)) {
#   install.packages("GWASinspector")
# }
    
library(GWASinspector)

# Command line input
args <- commandArgs(trailingOnly = TRUE)
config <- args[1]

### START DEVELOPMENT ##########################################################
# config <- "results/gwas-inspector/config/abdominal_config.ini"
# config <- "results/gwas-inspector/config/femoral_config.ini"
# config <- "results/gwas-inspector/config/ventral_config.ini"
# config <- "results/gwas-inspector/config/umbilical_config.ini"
# config <- "results/gwas-inspector/config/inguinal_config.ini"
# config <- "results/gwas-inspector/config/diaphragmatic_config.ini"
### END DEVELOPEMENT

job <- setup_inspector(config)
run_inspector(job)

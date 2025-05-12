#!/bin/bash
#BSUB -J gwas_inspector[1-6]              # Job name and array index
#BSUB -o logs/gwas_inspector_%I.out  # Standard output
#BSUB -e logs/gwas_inspector_%I.err  # Standard error
#BSUB -n 1                           # Number of cores
#BSUB -R "rusage[mem=32GB]"          # Memory requirement (adjust as needed)
#BSUB -W 24:00                        # Runtime limit (hh:mm)
#BSUB -N                             # Send email at job end

# 1. Create the array of input files
#cd /project/damrauer_wachtel/Hernia-PRS

# 2. Select the file for this job index
myArray=($(ls -1 results/gwas-inspector/config/*.ini))
config_file="${myArray[${LSB_JOBINDEX} - 1]}"

# 3. Run your R script with the selected input
export R_LIBS_USER=$HOME/R/rocker-rstudio/bioconductor-tidyverse_3.17
export SINGULARITY_BIND="/lsf:/lsf, /project/:/project/, /appl/:/appl/, /lsf/:/lsf/, /scratch/:/scratch, /static:/static"
module load singularity
singularity exec /project/voltron/rstudio/containers/bioconductor-tidyverse_3.17.sif Rscript workflow/rules/scripts/qc/gwas_inspector.R "${config_file}"
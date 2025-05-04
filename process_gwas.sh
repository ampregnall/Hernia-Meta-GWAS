#!/bin/bash
#BSUB -J process_gwas[1-55]        # Job name and array index
#BSUB -o logs/process_gwas_%I.out  # Standard output
#BSUB -q damrauer_normal           # Set queue name
#BSUB -e logs/process_gwas_%I.err  # Standard error
#BSUB -n 1                         # Number of cores
#BSUB -R "rusage[mem=20GB]"        # Memory requirement (adjust as needed)
#BSUB -W 1:00                      # Runtime limit (hh:mm)
#BSUB -B                           # Send email at job start
#BSUB -N                           # Send email at job end

# 1. Create the array of input files
myArray=($(ls -1 /data/raw/gwas-summary-statistics/*.txt))

# 2. Select the file for this job index
input_file="${myArray[${LSB_JOBINDEX} - 1]}"

# 3. Run your R script with the selected input
module load R/4.3
cd /project/damrauer_wachtel/Hernia-PRS
Rscript /workflow/rules/scripts/input/process_gwas.R "${input_file}"
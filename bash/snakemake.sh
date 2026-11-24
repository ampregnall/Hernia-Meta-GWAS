#!/bin/bash
#BSUB -J snakemake                   # Job name and array index
#BSUB -o logs/snakemake.out          # Standard output
#BSUB -e logs/snakemake.err          # Standard error
#BSUB -n 1                           # Number of cores
#BSUB -R "rusage[mem=4GB]"           # Memory requirement (adjust as needed)
#BSUB -W 24:00                       # Runtime limit (hh:mm)
#BSUB -N                             # Send email at job end

# Submit from the repository root: bsub < bash/snakemake.sh
# Snakemake and the LSF executor plugin are installed by pixi (pixi.toml).
cd "${LS_SUBCWD:-$(pwd)}"

pixi run workflow

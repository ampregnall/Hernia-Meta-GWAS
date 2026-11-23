#!/bin/bash
#BSUB -J mr_mega                     # Job name
#BSUB -o logs/mr_mega.out            # Standard output
#BSUB -e logs/mr_mega.err            # Standard error
#BSUB -n 1                           # Number of cores
#BSUB -R "rusage[mem=32GB]"          # Memory requirement (adjust as needed)
#BSUB -W 24:00                       # Runtime limit (hh:mm)
#BSUB -N                             # Send email at job end

# Multi-population meta-regression with MR-MEGA. Run after the Snakemake workflow,
# which produces the MAF-filtered summary statistics used as input.
# Requires a conda environment named "gwaslab" created from workflow/rules/envs/gwaslab.yaml:
#   conda env create -f workflow/rules/envs/gwaslab.yaml

cd /project/damrauer_wachtel/Hernia-Meta-GWAS

module load MR-MEGA/0.2
module load miniconda
source activate gwaslab

PHENOTYPES=(Abdominal Ventral Inguinal Diaphragmatic Umbilical Femoral)
INPUT_DIR=data/processed/sumstats-maf-filtered
OUT_DIR=data/processed/mr-mega
LEAD_DIR=data/processed/lead-variants/mr-mega
FIG_DIR=results/mqq-plots/mr-mega

mkdir -p ${OUT_DIR}/config ${LEAD_DIR} ${FIG_DIR}

for PHENOTYPE in "${PHENOTYPES[@]}"; do
    META=${PHENOTYPE}-ALL
    ls -1 ${INPUT_DIR}/*_${PHENOTYPE}_*_LDSC.txt > ${OUT_DIR}/config/${META}.in

    # MR-MEGA requires more cohorts than PCs + 2; use the maximum number of PCs (K - 3)
    NUM_DATASETS=$(wc -l < ${OUT_DIR}/config/${META}.in)
    NUM_PCS=$((NUM_DATASETS - 3))

    if [ "$NUM_PCS" -lt 1 ]; then
        echo "Skipping ${META}: ${NUM_DATASETS} cohorts is too few for MR-MEGA"
        continue
    fi

    echo "Running MR-MEGA for ${META} with ${NUM_DATASETS} cohorts and ${NUM_PCS} PCs"

    MR-MEGA --filelist ${OUT_DIR}/config/${META}.in --out ${OUT_DIR}/${META} \
        --qt \
        --name_pos POS \
        --name_chr CHR \
        --name_n N_CASE \
        --name_beta BETA \
        --name_eaf EAF \
        --name_nea NEA \
        --name_ea EA \
        --name_marker SNPID \
        --name_se SE \
        --pc $NUM_PCS \
        --no_std_names

    python src/mr-mega/lead_vars_mqq_mrmega.py \
        --input ${OUT_DIR}/${META}.result \
        --lead-vars ${LEAD_DIR}/${META}_lead_variants.txt \
        --fig ${FIG_DIR}/${META}_mqq_plot.png
done

# Hernia Meta-GWAS

Code for **"Decoding the genetic architecture of hernia through genome-wide association and multi-trait analyses"** (Pregnall *et al.*).

We performed multi-population genome-wide association meta-analyses of five hernia subtypes (ventral, inguinal, umbilical, femoral, and diaphragmatic) using summary statistics from the All of Us Research Program, FinnGen, BioBank Japan, the Pan-UK Biobank, and the Million Veteran Program. Downstream analyses include statistical finemapping, multi-trait and QTL colocalization, variant-to-gene prioritization, heritability and genetic correlation estimation, local genetic correlation (LAVA), genomic structural equation modeling (Genomic-SEM), and Mendelian randomization.

## Repository structure

```
.
├── config/
│   ├── config.yaml          # Workflow parameters and helper script paths
│   └── studies.txt          # Input summary statistics: study, phenotype, population, genome build, cases, controls
├── profile/default/         # Snakemake profile for the LSF cluster
├── workflow/                # Snakemake workflow: harmonization through meta-analysis and annotation
│   ├── Snakefile
│   └── rules/
│       ├── input.smk        # Summary statistic cleaning, liftover, and harmonization (MungeSumstats)
│       ├── ldsc.smk         # LDSC intercept correction (gwaslab)
│       ├── meta.smk         # MAF filter and METAL meta-analysis
│       ├── finemapping.smk  # Lead variants, Bayes factor finemapping, multi-trait colocalization (HyPrColoc)
│       ├── annotation.smk   # MAGMA, PoPS, eQTL colocalization, Ensembl VEP
│       ├── envs/            # Conda environments
│       └── scripts/         # Scripts called by workflow rules
├── src/                     # Downstream analyses run after the workflow
│   ├── meta-analysis/       # Population-specific lead variants
│   ├── mr-mega/             # MR-MEGA lead variants and Manhattan/QQ plots
│   ├── gene-prioritization/ # Annotation of credible sets with eQTL, PoPS, and VEP evidence
│   ├── ldsc/                # SNP heritability
│   ├── genomic-sem/         # Genetic correlation, exploratory and confirmatory factor analysis
│   ├── lava/                # Local genetic correlation
│   ├── mendelian-randomization/  # Risk factor, proteome-wide, and drug target MR; pQTL colocalization
│   ├── figures/             # Manuscript figures
│   └── envs/                # Conda environments for downstream analyses
└── bash/                    # LSF job scripts
    ├── snakemake.sh         # Runs the Snakemake workflow
    ├── mr_mega.sh           # Multi-population meta-regression with MR-MEGA
    └── gwas_inspector.sh    # Quality control metrics for input summary statistics
```

## Analyses

| Analysis | Code | Manuscript |
|---|---|---|
| Harmonization, LDSC correction, MAF filter | `workflow/rules/input.smk`, `ldsc.smk`, `meta.smk` | Methods: Genome-wide association meta-analysis |
| Fixed-effects meta-analysis (METAL) | `workflow/rules/meta.smk` | Fig. 2 |
| Multi-population meta-regression (MR-MEGA) | `bash/mr_mega.sh`, `src/mr-mega/` | Supplementary Fig. 3 |
| Lead variants and shared loci | `workflow/rules/finemapping.smk`, `src/figures/lead_loci_upset.R` | Supplementary Fig. 1 |
| Population-specific lead variants | `src/meta-analysis/population_specific_lead_variants.R` | Supplementary Fig. 2 |
| Bayes factor finemapping | `workflow/rules/finemapping.smk` | Fig. 3 |
| Multi-trait colocalization | `workflow/rules/finemapping.smk` | Fig. 4 |
| Variant-to-gene prioritization (VEP, eQTL colocalization, PoPS) | `workflow/rules/annotation.smk`, `src/gene-prioritization/` | Fig. 3 |
| SNP heritability | `src/ldsc/` | Supplementary Fig. 4 |
| Genetic correlation and Genomic-SEM | `src/genomic-sem/` | Figs. 4 and 5 |
| Local genetic correlation (LAVA) | `src/lava/` | Fig. 4 |
| Mendelian randomization | `src/mendelian-randomization/` | Figs. 6 and 7 |

## Running the analyses

### Requirements

The workflow was developed for an LSF high-performance computing cluster and uses [Snakemake](https://snakemake.readthedocs.io) (≥ 8, with the LSF executor plugin) and conda. Some rules load software through environment modules (METAL, MAGMA, MR-MEGA, Apptainer for VEP). Cluster-specific settings and paths are in `profile/default/config.yaml`, `bash/`, and the `shell` blocks of the workflow rules, and will need to be adapted for other systems.

Input summary statistics and reference data are expected in `data/`, which is not tracked in this repository:

- `data/raw/sumstats/`: input summary statistics named `{study}_{phenotype}_{population}.txt`, matching `config/studies.txt`
- `data/meta/`: reference data, including the hg19-to-hg38 liftover chain file, dbSNP rsID references, Pan-UK Biobank LD scores, MAGMA and PoPS reference files, eQTL Catalogue fibroblast eQTLs, the LAVA reference panel, and Mendelian randomization instruments

### Meta-analysis workflow

```bash
# Submit the workflow to LSF
bsub < bash/snakemake.sh

# Or run directly from the repository root
snakemake --use-conda --profile profile/default
```

### Downstream analyses

Downstream scripts are run from the repository root after the workflow completes, for example `Rscript src/ldsc/heritability.R`. MR-MEGA is run with `bsub < bash/mr_mega.sh`.

## Input data

GWAS summary statistics were obtained from:

| Source | Populations | Access |
|---|---|---|
| FinnGen (release 12) | EUR | https://www.finngen.fi/en/access_results |
| Pan-UK Biobank | EUR, AFR, AMR, CSA, EAS, MID | https://pan.ukbb.broadinstitute.org |
| Million Veteran Program | EUR, AFR, AMR | dbGaP |
| All of Us Research Program | EUR, AFR, AMR | https://www.researchallofus.org |
| BioBank Japan | EAS | https://pheweb.jp |

Case and control counts for each study are listed in `config/studies.txt`. Phenotype definitions are described in the Supplementary Data of the manuscript.

## Data availability

Meta-analysis summary statistics for ventral, inguinal, umbilical, femoral, diaphragmatic, and abdominal hernia will be deposited in the [GWAS Catalog](https://www.ebi.ac.uk/gwas/) upon publication.

| Phenotype | GWAS Catalog accession |
|---|---|
| Ventral hernia | To be added upon publication |
| Inguinal hernia | To be added upon publication |
| Umbilical hernia | To be added upon publication |
| Femoral hernia | To be added upon publication |
| Diaphragmatic hernia | To be added upon publication |
| Abdominal hernia | To be added upon publication |

## Software

| Software | Version |
|---|---|
| R | 4.3 |
| MungeSumstats | Bioconductor |
| gwaslab | 3.6.3 |
| METAL | 2020-05-05 |
| MR-MEGA | 0.2 |
| HyPrColoc | 0.0.2 |
| MAGMA | 1.07a |
| PoPS | [FinucaneLab/pops](https://github.com/FinucaneLab/pops) |
| Ensembl VEP | 113 |
| GenomicSEM | GitHub |
| LAVA | 0.1.0 |
| TwoSampleMR | 0.7.4 |

Conda environment specifications are in `workflow/rules/envs/` and `src/envs/`.

## Citation

If you use this code, please cite:

> Pregnall AM, *et al.* Decoding the genetic architecture of hernia through genome-wide association and multi-trait analyses. *Manuscript under review.*

```bibtex
@misc{pregnallDecodingGeneticArchitecture2026,
  title = {Decoding the Genetic Architecture of Hernia through Genome-Wide Association and Multi-Trait Analyses},
  author = {Pregnall, Andrew M. and Yuan, Shuai and Lawrence, Jeremy M and Abramowitz, Sarah A. and DePaolo, John M and Judy, Renae and Shakt, Gabrielle and Levin, Michael and Damrauer, Scott M. and Wachtel, Heather},
  year = 2026,
  month = may,
  publisher = {{Genetic and Genomic Medicine}},
  doi = {10.64898/2026.05.12.26353033},
  urldate = {2026-05-16},
  archiveprefix = {Genetic and Genomic Medicine},
  copyright = {http://creativecommons.org/licenses/by/4.0/},
  langid = {english}
}

```

## Contact

Questions about the code can be submitted as a [GitHub issue](https://github.com/ampregnall/Hernia-Meta-GWAS/issues).

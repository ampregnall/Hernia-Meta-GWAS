# Workflow configuration

## `studies.txt`

One row per input set of summary statistics. Each row must have a matching file at `data/raw/sumstats/{study}_{phenotype}_{population}.txt`.

| Column | Description |
|---|---|
| `study` | Source study (`FinnGen-R12`, `MVP`, `Pan-UKBB`, `AOU`, `BJP`); selects the study-specific cleaning function in `workflow/rules/scripts/input/process_gwas_utils.R` |
| `phenotype` | Hernia subtype (`Ventral`, `Inguinal`, `Umbilical`, `Femoral`, `Diaphragmatic`, `Abdominal`) |
| `population` | Most similar 1000 Genomes reference population (`EUR`, `AFR`, `AMR`, `CSA`, `EAS`, `MID`); selects the Pan-UK Biobank LD scores used for LDSC correction |
| `genome` | Genome build of the input file (`GrCh37` or `GrCh38`) |
| `cases`, `controls` | Number of cases and controls |
| `n_effective` | Effective sample size, 4 / (1/cases + 1/controls) |

## `config.yaml`

| Key | Description |
|---|---|
| `process_gwas` | Cleaning helpers, liftover chain file, dbSNP SNPlocs tarball for MungeSumstats, and temporary directory for unpacking the tarball |
| `ldsc_scores` | Directory of Pan-UK Biobank LD scores (`{population}/UKBB.{population}`) for LDSC intercept correction |
| `mac_filter` | Minimum expected minor allele count in cases, applied after LDSC correction |
| `meta_analysis` | Trait type, chromosomes (23 = X), and phenotypes and populations to meta-analyze. A population-specific meta-analysis is run when at least two studies are available for that phenotype and population |
| `downstream_phenotypes` | Hernia subtypes carried forward to finemapping, colocalization, gene prioritization, and LDSC-SEG |
| `ldsc_seg` | Population, annotation file, baseline LD scores, weights, and HapMap3 SNP list for tissue and cell type enrichment |
| `hyprcoloc`, `bayesian_finemapping` | Helper functions and parameters for colocalization and finemapping |
| `magma`, `pops`, `eqtl`, `vep` | Reference data and parameters for gene prioritization |
| `modules` | Environment modules loaded by rules (MAGMA, Apptainer) |

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
| `process_gwas` | Helper functions and liftover chain file for summary statistic cleaning |
| `ldsc_correction` | dbSNP reference files used by gwaslab to assign rsIDs |
| `create_meta_analysis_configs` | Phenotypes and populations to meta-analyze. A population-specific meta-analysis is run when at least two studies are available for that phenotype and population |
| `metal_format` | Hernia subtypes carried forward to finemapping, colocalization, and gene prioritization |
| `hyprcoloc`, `bayesian_finemapping` | Helper functions for colocalization and finemapping |

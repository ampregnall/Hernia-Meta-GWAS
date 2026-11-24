# Study metadata, meta-analysis definitions, and helper functions
import os

import pandas as pd

# The tmpdir resource is not created by Snakemake; if it is missing, R falls back to
# /tmp and MungeSumstats unpacks the dbSNP tarball onto the node's local disk.
os.makedirs(config["process_gwas"]["tmpdir"], exist_ok=True)


# --------------
# Input summary statistics
# --------------
studies = pd.read_table("config/studies.txt")
studies["id"] = studies.apply(
    lambda r: f"{r.study}_{r.phenotype}_{r.population}", axis=1
)
metadata = studies.set_index("id").to_dict(orient="index")
ids = list(metadata)


def ldsc_scores(wildcards):
    """Pan-UK Biobank LD scores for the study's population."""
    population = metadata[wildcards.id]["population"]
    return f"{config['ldsc_scores']}/{population}/UKBB.{population}"


# --------------
# Meta-analyses
# --------------
def studies_for(phenotype, population):
    """Summary statistic ids included in one meta-analysis."""
    sel = studies["phenotype"] == phenotype
    if population != "ALL":
        sel &= studies["population"] == population
    return studies.loc[sel, "id"].tolist()


# valid_pairs: multi-population (ALL) and population-specific meta-analyses with >= 2
#              sumstats. These go on to lead variants and finemapping.
# meta_pairs:  valid_pairs plus single-study EAS sumstats, which are passed through the
#              meta-analysis so they share the same output format.
valid_pairs = []
meta_pairs = []

for _phenotype in config["meta_analysis"]["phenotypes"]:
    for _population in config["meta_analysis"]["populations"]:
        _n = len(studies_for(_phenotype, _population))
        if _population == "ALL" or _n >= 2:
            valid_pairs.append((_phenotype, _population))
            meta_pairs.append((_phenotype, _population))
        elif _population == "EAS" and _n == 1:
            meta_pairs.append((_phenotype, _population))


def pair_paths(pattern, pairs):
    """Expand a {phenotype}-{population} pattern over meta-analysis pairs."""
    return [pattern.format(phenotype=ph, population=pop) for ph, pop in pairs]


downstream = config["downstream_phenotypes"]

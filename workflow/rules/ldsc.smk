# Load metadata
df = pd.read_table("config/studies.txt")
df["id"] = df.apply(lambda row: f"{row.study}_{row.phenotype}_{row.population}", axis=1)
df["ldsc"] = df.apply(lambda row: f"data/meta/gwaslab-references/ldsc/{row.population}/UKBB.{row.population}", axis=1)
metadata = {row.id: row._asdict() for row in df.itertuples(index=False)}
ids = list(metadata.keys())

rule ldsc_correction:
    input:
        "data/processed/gwas-summary-statistics/{id}.txt"
    output:
        sumstats = "data/processed/sumstats-ldsc-corrected/{id}_LDSC.txt",
        ldsc = "data/processed/ldsc-stats/{id}.txt"
    params:
        rsid_tsv = config["ldsc_correction"]["rsid_tsv"],
        rsid_vcf = config["ldsc_correction"]["rsid_vcf"],
        ldsc_scores = lambda wildcards : metadata[wildcards.id]["ldsc"]
    conda:
        "gwaslab"
    threads: 2
    resources:
        mem_mb = 16000,
        walltime = 120
    script:
        "scripts/ldsc/ldsc_correction.py"
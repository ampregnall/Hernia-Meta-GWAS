# Load metadata
df = pd.read_table("config/studies.txt")
df["id"] = df.apply(lambda row: f"{row.study}_{row.phenotype}_{row.population}", axis=1)
metadata = {row.id: row._asdict() for row in df.itertuples(index=False)}
ids = list(metadata.keys())

rule process_gwas_sumstats:
    input:
        "data/raw/gwas-summary-statistics/{id}.txt"
    output:
        "data/processed/gwas-summary-statistics/{id}.txt"
    params:
        genome = lambda wildcards: metadata[wildcards.id]["genome"],
        study = lambda wildcards: metadata[wildcards.id]["study"],
        n_cases = lambda wildcards: metadata[wildcards.id]["cases"],
        n_controls = lambda wildcards: metadata[wildcards.id]["controls"],
        utils = config["process_gwas"]["utils"],
        chainfile = config["process_gwas"]["chainfile"],
    script:
        "scripts/input/process_gwas.R"
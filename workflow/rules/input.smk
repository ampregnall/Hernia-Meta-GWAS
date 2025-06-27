rule process_gwas_sumstats:
    input:
        "data/raw/sumstats/{id}.txt"
    output:
        "data/processed/sumstats-harmonized/{id}.txt"
    params:
        genome = lambda wildcards: metadata[wildcards.id]["genome"],
        study = lambda wildcards: metadata[wildcards.id]["study"],
        n_cases = lambda wildcards: metadata[wildcards.id]["cases"],
        n_controls = lambda wildcards: metadata[wildcards.id]["controls"],
        utils = config["process_gwas"]["utils"],
        chainfile = config["process_gwas"]["chainfile"]
    conda:
        "envs/common.yaml"
    script:
        "scripts/input/process_gwas.R"
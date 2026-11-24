rule run_hyprcoloc:
    input:
        sumstats=expand(
            "data/processed/meta-formatted/{phenotype}-ALL_meta_formatted.txt",
            phenotype=downstream,
        ),
        lead=expand(
            "data/processed/lead-variants/meta/{phenotype}-ALL_lead_variants.txt",
            phenotype=downstream,
        ),
    output:
        colocalization="data/processed/finemapping/HyPrColoc_results.txt",
        credset="data/processed/finemapping/HyPrColoc_finemapping.txt",
    log:
        "logs/run_hyprcoloc.log",
    benchmark:
        "benchmarks/run_hyprcoloc.tsv"
    params:
        phenotypes=downstream,
        creds=config["hyprcoloc"]["creds"],
        locus=config["hyprcoloc"]["locus"],
    conda:
        "envs/hyprcoloc.yaml"
    resources:
        mem_mb=120000,
    script:
        "scripts/finemapping/run_hyprcoloc.R"


rule bayesian_finemapping:
    input:
        sumstats="data/processed/meta-formatted/{phenotype}-{population}_meta_formatted.txt",
        lead="data/processed/lead-variants/meta/{phenotype}-{population}_lead_variants.txt",
    output:
        credset="data/processed/finemapping/{phenotype}-{population}_bayesian_finemapping.txt",
    log:
        "logs/bayesian_finemapping/{phenotype}-{population}.log",
    benchmark:
        "benchmarks/bayesian_finemapping/{phenotype}-{population}.tsv"
    params:
        creds=config["bayesian_finemapping"]["creds"],
        locus=config["bayesian_finemapping"]["locus"],
        credible_interval=config["bayesian_finemapping"]["credible_interval"],
        locus_size=config["bayesian_finemapping"]["locus_size"],
    conda:
        "envs/common.yaml"
    script:
        "scripts/finemapping/bayesian_finemapping.R"

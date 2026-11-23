rule identify_lead_variants_metal:
    input:
        "data/processed/metal/{phenotype}_metal-1.txt"
    output:
        lead_vars = "data/processed/lead-variants/metal/{phenotype}_lead_variants.txt",
        fig = "results/mqq-plots/metal/{phenotype}_mqq_plot.png",
    params:
        rsid_tsv = config["ldsc_correction"]["rsid_tsv"],
        rsid_vcf = config["ldsc_correction"]["rsid_vcf"]
    conda:
        "envs/gwaslab.yaml"
    script:
        "scripts/meta-analysis/lead_vars_mqq_metal.py"

rule format_metal:
    input:
        "data/processed/metal/{phenotype}_metal-1.txt"
    output:
        sumstats = "data/processed/metal-formatted/{phenotype}_metal_formatted.txt"
    params:
        rsid_tsv = config["ldsc_correction"]["rsid_tsv"],
        rsid_vcf = config["ldsc_correction"]["rsid_vcf"]
    conda:
        "envs/gwaslab.yaml"
    threads: 4
    script:
        "scripts/finemapping/format_metal.py"

rule run_hyprcoloc:
    input:
        expand("data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt", phenotype=config["metal_format"]["phenotypes"])
    output:
        colocalization = "data/processed/finemapping/HyPrColoc_results.txt",
        credset = "data/processed/finemapping/HyPrColoc_finemapping.txt"
    conda:
        "envs/hyprcoloc.yaml"
    params:
        creds = config["hyprcoloc"]["creds"],
        locus = config["hyprcoloc"]["locus"]
    resources:
        mem_mb = 120000
    script:
        "scripts/finemapping/run_hyprcoloc.R"

rule bayesian_finemapping:
    input:
        sumstats = "data/processed/metal-formatted/{phenotype}_metal_formatted.txt",
        lead = "data/processed/lead-variants/metal/{phenotype}_lead_variants.txt"
    output:
        credset = "data/processed/finemapping/{phenotype}_bayesian_finemapping.txt"
    params:
        creds = config["bayesian_finemapping"]["creds"],
        locus = config["bayesian_finemapping"]["locus"]
    conda:
        "envs/common.yaml"
    script:
        "scripts/finemapping/bayesian_finemapping.R"

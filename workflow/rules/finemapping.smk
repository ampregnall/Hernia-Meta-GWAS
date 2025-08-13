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

rule identify_lead_variants_mrmega:
    input:
        "data/processed/mr-mega/{phenotype}.result"
    output:
        lead_vars = "data/processed/lead-variants/mr-mega/{phenotype}_lead_variants.txt",
        fig = "results/mqq-plots/mr-mega/{phenotype}_mqq_plot.png"
    conda:
        "envs/gwaslab.yaml"
    script:
        "scripts/meta-analysis/lead_vars_mqq_mrmega.py"

rule format_metal:
    input:
        "data/processed/metal/{phenotype}-ALL_metal-1.txt"
    output:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt"
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

rule run_rsparsepro:
    input:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt",
        lead = "data/processed/lead-variants/metal/{phenotype}-ALL_lead_variants.txt",
    output:
        rsparsepro = "data/processed/finemapping/RSparsePro_{phenotype}_finemapping.txt",
        stats = "data/processed/finemapping/RSparsePro_{phenotype}_stats.txt"
    params:
        utils = "workflow/rules/scripts/finemapping/locus_extract.R",
        rsparsepro = "workflow/rules/scripts/finemapping/rsparsepro_ld.py",
        bed = "/project/voltron/Resources/1000-Genomes-Phase3/1kg_all_phase3/1kg_all_phase3.rds"
    resources:
        mem_mb = 64000
    conda:
        "envs/rsparsepro.yaml"
    script:
        "scripts/finemapping/rsparsepro_reticulate.R"

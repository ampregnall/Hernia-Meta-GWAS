rule meta_analysis:
    input:
        lambda wc: expand(
            "data/processed/sumstats-ldsc-corrected/{id}.parquet",
            id=studies_for(wc.phenotype, wc.population),
        ),
    output:
        "data/processed/meta-analysis/{phenotype}-{population}_meta.parquet",
    log:
        "logs/meta_analysis/{phenotype}-{population}.log",
    benchmark:
        "benchmarks/meta_analysis/{phenotype}-{population}.tsv"
    params:
        trait_type=config["meta_analysis"]["trait_type"],
        chromosomes=config["meta_analysis"]["chromosomes"],
    conda:
        "envs/meta.yaml"
    resources:
        mem_mb=32000,
    script:
        "scripts/meta-analysis/meta_analysis.R"


rule format_meta:
    input:
        "data/processed/meta-analysis/{phenotype}-{population}_meta.parquet",
    output:
        sumstats="data/processed/meta-formatted/{phenotype}-{population}_meta_formatted.txt",
    log:
        "logs/format_meta/{phenotype}-{population}.log",
    benchmark:
        "benchmarks/format_meta/{phenotype}-{population}.tsv"
    conda:
        "envs/gwaslab.yaml"
    script:
        "scripts/meta-analysis/format_meta.py"


rule identify_lead_variants:
    input:
        "data/processed/meta-formatted/{phenotype}-{population}_meta_formatted.txt",
    output:
        lead_vars="data/processed/lead-variants/meta/{phenotype}-{population}_lead_variants.txt",
        fig="results/mqq-plots/meta/{phenotype}-{population}_mqq_plot.png",
    log:
        "logs/identify_lead_variants/{phenotype}-{population}.log",
    benchmark:
        "benchmarks/identify_lead_variants/{phenotype}-{population}.tsv"
    conda:
        "envs/gwaslab.yaml"
    script:
        "scripts/meta-analysis/lead_vars_mqq_meta.py"

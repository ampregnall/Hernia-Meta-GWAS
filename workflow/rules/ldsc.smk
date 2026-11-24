rule ldsc_correction:
    input:
        sumstats="data/processed/sumstats-harmonized/{id}.txt",
        counts="data/processed/qc/variant-counts/{id}_harmonization.tsv",
    output:
        sumstats="data/processed/sumstats-ldsc-corrected/{id}.parquet",
        ldsc="data/processed/ldsc-stats/{id}.txt",
        counts="data/processed/qc/variant-counts/{id}.tsv",
    log:
        "logs/ldsc_correction/{id}.log",
    benchmark:
        "benchmarks/ldsc_correction/{id}.tsv"
    params:
        ldsc_scores=ldsc_scores,
        min_case_mac=config["mac_filter"]["min_case_mac"],
    conda:
        "envs/gwaslab.yaml"
    threads: 2
    resources:
        mem_mb=16000,
        walltime=120,
    script:
        "scripts/ldsc/ldsc_correction.py"


rule process_gwas_sumstats:
    input:
        sumstats="data/raw/sumstats/{id}.txt",
        dbsnp_tarball=config["process_gwas"]["dbsnp_tarball"],
    output:
        sumstats="data/processed/sumstats-harmonized/{id}.txt",
        counts="data/processed/qc/variant-counts/{id}_harmonization.tsv",
    log:
        "logs/process_gwas_sumstats/{id}.log",
    benchmark:
        "benchmarks/process_gwas_sumstats/{id}.tsv"
    params:
        study=lambda wc: metadata[wc.id]["study"],
        n_cases=lambda wc: metadata[wc.id]["cases"],
        n_controls=lambda wc: metadata[wc.id]["controls"],
        utils=config["process_gwas"]["utils"],
        chainfile=config["process_gwas"]["chainfile"],
    resources:
        mem_mb=64000,
        tmpdir=config["process_gwas"]["tmpdir"],
    conda:
        "envs/munge.yaml"
    script:
        "scripts/input/process_gwas.R"


rule summarise_variant_counts:
    input:
        expand("data/processed/qc/variant-counts/{id}.tsv", id=ids),
    output:
        "results/qc/variant_counts.tsv",
    log:
        "logs/summarise_variant_counts.log",
    localrule: True
    run:
        counts = pd.concat([pd.read_table(f) for f in input])
        counts.to_csv(output[0], sep="\t", index=False)

rule ldsc_correction:
    input:
        "data/processed/sumstats-harmonized/{id}.txt"
    output:
        sumstats = "data/processed/sumstats-ldsc-corrected/{id}_LDSC.txt",
        ldsc = "data/processed/ldsc-stats/{id}.txt"
    params:
        rsid_tsv = config["ldsc_correction"]["rsid_tsv"],
        rsid_vcf = config["ldsc_correction"]["rsid_vcf"],
        ldsc_scores = lambda wildcards : metadata[wildcards.id]["ldsc"]
    conda:
        "envs/gwaslab.yaml"
    threads: 2
    resources:
        mem_mb = 16000,
        walltime = 120
    script:
        "scripts/ldsc/ldsc_correction.py"
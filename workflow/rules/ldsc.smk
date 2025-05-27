rule munge_sumstats:
    input:
        "data/processed/gwas-summary-statistics/{id}.txt"
    output:
        "data/processed/ldsc-summary-statistics/{id}_LDSC.txt"
    params:
        alleles = lambda wildcards: config["munge_sumstats"]["alleles"],
        chunksize = config["munge_sumstats"]["chunksize"]
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python workflow/rules/scripts/ldsc/munge_sumstats.py \
            --sumstats {input} \
            --out {output} \
            --merge-alleles {params.alleles} \
            --chunksize {params.chunksize}
        """
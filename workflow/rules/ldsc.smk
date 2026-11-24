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


# --------------
# Tissue and cell type enrichment (LDSC-SEG)
# --------------
rule munge_ldsc:
    input:
        sumstats="data/processed/meta-formatted/{phenotype}-{population}_meta_formatted.txt",
        hm3=config["ldsc_seg"]["hm3_snplist"],
    output:
        selected=temp("data/processed/ldsc-seg/munged/{phenotype}-{population}_input.txt"),
        munged="data/processed/ldsc-seg/munged/{phenotype}-{population}.sumstats.gz",
    log:
        "logs/munge_ldsc/{phenotype}-{population}.log",
    benchmark:
        "benchmarks/munge_ldsc/{phenotype}-{population}.tsv"
    params:
        out_prefix=lambda wc: f"data/processed/ldsc-seg/munged/{wc.phenotype}-{wc.population}",
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        # Select columns by header name; LDSC auto-detection would see both SNPID and
        # rsID as SNP columns
        awk -F'\\t' -v OFS='\\t' '
            NR == 1 {{ for (i = 1; i <= NF; i++) col[$i] = i; print "SNP", "A1", "A2", "P", "BETA", "N"; next }}
            $col["rsID"] != "" {{ print $col["rsID"], $col["EA"], $col["NEA"], $col["P"], $col["BETA"], $col["N"] }}
        ' {input.sumstats} > {output.selected}

        munge_sumstats.py \
            --sumstats {output.selected} \
            --merge-alleles {input.hm3} \
            --signed-sumstats BETA,0 \
            --chunksize 500000 \
            --out {params.out_prefix} > {log} 2>&1
        """


rule ldsc_seg:
    input:
        sumstats="data/processed/ldsc-seg/munged/{phenotype}-{population}.sumstats.gz",
        ldcts=config["ldsc_seg"]["ldcts"],
    output:
        "data/processed/ldsc-seg/{phenotype}-{population}.cell_type_results.txt",
    log:
        "logs/ldsc_seg/{phenotype}-{population}.log",
    benchmark:
        "benchmarks/ldsc_seg/{phenotype}-{population}.tsv"
    params:
        baseline=config["ldsc_seg"]["baseline"],
        weights=config["ldsc_seg"]["weights"],
        out_prefix=lambda wc: f"data/processed/ldsc-seg/{wc.phenotype}-{wc.population}",
    conda:
        "envs/ldsc.yaml"
    resources:
        mem_mb=16000,
    shell:
        """
        ldsc.py \
            --h2-cts {input.sumstats} \
            --ref-ld-chr {params.baseline} \
            --ref-ld-chr-cts {input.ldcts} \
            --w-ld-chr {params.weights} \
            --out {params.out_prefix} > {log} 2>&1
        """

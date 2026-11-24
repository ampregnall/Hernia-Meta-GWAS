rule run_magma:
    input:
        sumstats="data/processed/meta-formatted/{phenotype}-ALL_meta_formatted.txt",
    output:
        infile="data/processed/magma/{phenotype}/{phenotype}.txt",
        out="data/processed/magma/{phenotype}/{phenotype}.genes.out",
        raw="data/processed/magma/{phenotype}/{phenotype}.genes.raw",
    log:
        "logs/run_magma/{phenotype}.log",
    benchmark:
        "benchmarks/run_magma/{phenotype}.tsv"
    params:
        bfile=config["magma"]["bfile"],
        gene_annot=config["magma"]["gene_annot"],
        out_prefix=lambda wc: f"data/processed/magma/{wc.phenotype}/{wc.phenotype}",
    envmodules:
        config["modules"]["magma"],
    shell:
        """
        # SNP (rsID), P, and N columns, selected by header name
        awk -F'\\t' -v OFS='\\t' '
            NR == 1 {{ for (i = 1; i <= NF; i++) col[$i] = i; print "SNP", "P", "N"; next }}
            $col["rsID"] != "" {{ print $col["rsID"], $col["P"], $col["N"] }}
        ' {input.sumstats} > {output.infile}

        magma --bfile {params.bfile} \
              --gene-annot {params.gene_annot} \
              --pval {output.infile} use=SNP,P ncol=N \
              --gene-model snp-wise=mean \
              --out {params.out_prefix} > {log} 2>&1
        """


rule run_pops:
    input:
        magma_out="data/processed/magma/{phenotype}/{phenotype}.genes.out",
        magma_raw="data/processed/magma/{phenotype}/{phenotype}.genes.raw",
    output:
        preds="data/processed/pops/{phenotype}.preds",
        coefs="data/processed/pops/{phenotype}.coefs",
        marginals="data/processed/pops/{phenotype}.marginals",
    log:
        "logs/run_pops/{phenotype}.log",
    benchmark:
        "benchmarks/run_pops/{phenotype}.tsv"
    params:
        out_prefix=lambda wc: f"data/processed/pops/{wc.phenotype}",
        magma_prefix=lambda wc: f"data/processed/magma/{wc.phenotype}/{wc.phenotype}",
        gene_annot=config["pops"]["gene_annot"],
        control_features=config["pops"]["control_features"],
        feature_mat_prefix=config["pops"]["feature_mat_prefix"],
        num_feature_chunks=config["pops"]["num_feature_chunks"],
    conda:
        "envs/pops.yaml"
    shell:
        """
        python workflow/rules/scripts/variant-annotation/pops.py \
            --gene_annot_path {params.gene_annot} \
            --feature_mat_prefix {params.feature_mat_prefix} \
            --num_feature_chunks {params.num_feature_chunks} \
            --magma_prefix {params.magma_prefix} \
            --control_features_path {params.control_features} \
            --out_prefix {params.out_prefix} > {log} 2>&1
        """


rule run_hyprcoloc_eqtl:
    input:
        sumstats="data/processed/meta-formatted/{phenotype}-ALL_meta_formatted.txt",
        lead="data/processed/lead-variants/meta/{phenotype}-ALL_lead_variants.txt",
        eqtl=config["eqtl"]["sumstats"],
    output:
        eqtl="data/processed/variant-annotation/eqtl/{phenotype}_eqtl_colocalization.txt",
    log:
        "logs/run_hyprcoloc_eqtl/{phenotype}.log",
    benchmark:
        "benchmarks/run_hyprcoloc_eqtl/{phenotype}.tsv"
    params:
        locus=config["hyprcoloc"]["locus"],
        min_overlap=config["eqtl"]["min_overlap"],
    conda:
        "envs/hyprcoloc.yaml"
    script:
        "scripts/variant-annotation/hyprcoloc_eqtl.R"


rule create_vep_input:
    input:
        sumstats="data/processed/meta-formatted/{phenotype}-{population}_meta_formatted.txt",
        lead="data/processed/lead-variants/meta/{phenotype}-{population}_lead_variants.txt",
    output:
        vep_input="data/processed/vep/input/{phenotype}-{population}.vep",
    log:
        "logs/create_vep_input/{phenotype}-{population}.log",
    benchmark:
        "benchmarks/create_vep_input/{phenotype}-{population}.tsv"
    params:
        vep_prefix=lambda wc: f"data/processed/vep/input/{wc.phenotype}-{wc.population}",
    conda:
        "envs/gwaslab.yaml"
    script:
        "scripts/variant-annotation/create_vep_input.py"


rule run_vep:
    input:
        vep="data/processed/vep/input/{phenotype}-{population}.vep",
    output:
        vep_annotation="data/processed/vep/output/{phenotype}-{population}.txt",
    log:
        "logs/run_vep/{phenotype}-{population}.log",
    benchmark:
        "benchmarks/run_vep/{phenotype}-{population}.tsv"
    params:
        container=config["vep"]["container"],
        cache_dir=config["vep"]["cache_dir"],
        bind=config["vep"]["bind"],
    envmodules:
        config["modules"]["apptainer"],
    shell:
        """
        if [ -s {input.vep} ]; then
            export APPTAINER_BIND="{params.bind}"
            apptainer exec {params.container} vep \
                --cache \
                --offline \
                --pick \
                --dir {params.cache_dir} \
                --input_file {input.vep} \
                --output_file {output.vep_annotation} \
                --species homo_sapiens \
                --assembly GRCh38 \
                --force_overwrite > {log} 2>&1
        else
            echo "No variants to annotate" > {log}
            touch {output.vep_annotation}
        fi
        """

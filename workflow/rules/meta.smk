rule filter_maf:
    input:
        "data/processed/sumstats-ldsc-corrected/{id}_LDSC.txt"    
    output:
        "data/processed/sumstats-maf-filtered/{id}_LDSC.txt"
    conda:
        "envs/common.yaml"
    script:
        "scripts/meta-analysis/filter_maf.R"

rule create_meta_analysis_configs:
    input:
        expand("data/processed/sumstats-maf-filtered/{id}_LDSC.txt", id=ids)
    output:
        metal = metal_config,
        mrmega = mrmega_config
    params:
        utils = config["create_meta_analysis_configs"]["utils"],
        data_dir = "data/processed/sumstats-maf-filtered",
        metal = "data/processed/metal",
        mrmega = "data/processed/mr-mega/config"
    conda:
        "envs/common.yaml"
    script:
        "scripts/meta-analysis/create_meta_analysis_configs.R"

rule run_metal:
    input:
        "data/processed/metal/{meta}_metal-config.txt"
    output:
        "data/processed/metal/{meta}_metal-1.txt"
    shell:
        """
        module load metal/2020.05
        metal {input}
        """

rule run_mrmega:
    input:
        "data/processed/mr-mega/config/{meta}.in"
    output:
        "data/processed/mr-mega/{meta}.result"
    shell:
        """
        module load MR-MEGA/0.2

        NUM_DATASETS=$(wc -l < {input})
        NUM_PCS=$((NUM_DATASETS - 3))

        if [ "$NUM_PCS" -lt 1 ]; then
           NUM_PCS=1
        fi

        MR-MEGA --filelist {input} --out data/processed/mr-mega/{wildcards.meta} \
            --qt \
            --name_pos POS \
            --name_chr CHR \
            --name_n N_CASE \
            --name_beta BETA \
            --name_eaf EAF \
            --name_nea NEA \
            --name_ea EA \
            --name_marker SNPID \
            --name_se SE \
            --pc $NUM_PCS \
            --no_std_names
        """

rule identify_lead_variants_metal:
    input:
        "data/processed/metal/{phenotype}_metal-1.txt"
    output:
        lead_vars = "data/processed/lead-variants/metal/{phenotype}_lead_variants.txt",
        fig = "results/mqq-plots/metal/{phenotype}_mqq_plot.png"
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

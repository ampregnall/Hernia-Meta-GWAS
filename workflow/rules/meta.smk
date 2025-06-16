configfile: "config/config.yaml"

df = pd.read_table("config/studies.txt")
df["id"] = df.apply(lambda row: f"{row.study}_{row.phenotype}_{row.population}", axis=1)
metadata = {row.id: row._asdict() for row in df.itertuples(index=False)}
ids = list(metadata.keys())
id_to_phenotype = {row.id: row.phenotype for row in df.itertuples(index=False)}

rule filter_maf:
    input:
        expand("results/gwas-inspector/{phenotype}-ldsc/QC_{id}_LDSC.txt", zip, phenotype=[id_to_phenotype[i].lower() for i in ids], id=ids)
    output:
        "data/processed/sumstats-maf-filtered/{id}.txt"
    params:
        utils = config["filter_maf"]["utils"]
    conda:
        "envs/common.yaml"
    script:
        "scripts/meta-analysis/filter_maf.R"

rule create_meta_analysis_configs:
    input:
        expand("data/processed/sumstats-maf-filtered/{id}.txt", id=ids)
    output:
        metal = expand("data/processed/metal/{phenotype}-{population}_metal-config.txt", 
            phenotype=config['create_metal_config']['phenotypes'],
            population=config['create_metal_config']['populations']),
        mrmega = expand("data/processed/mr-mega/config/{phenotype}-ALL.in", phenotype=config['create_metal_config']['phenotypes'])
    params:
        utils = config["create_metal_config"]["utils"],
        data_dir = "data/processed/sumstats-maf-filtered",
        metal = "data/processed/metal",
        mrmega = "data/processed/mr-mega/config"
    conda:
        "envs/common.yaml"
    script:
        "scripts/meta-analysis/create_metal_config.R"

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

        #NUM_DATASETS=$(wc -l < {input})
        #NUM_PCS=$((NUM_DATASETS - 3))

        #if [ "$NUM_PCS" -lt 1 ]; then
        #    NUM_PCS=1
        #fi

        MR-MEGA --filelist {input} --out data/processed/mr-mega/{wildcards.config} \
            --qt \
            --name_pos POSITION \
            --name_chr CHR \
            --name_n N_CASE \
            --name_beta EFFECT \
            --name_eaf EFF_ALL_FREQ \
            --name_nea OTHER_ALL \
            --name_ea EFFECT_ALL \
            --name_marker MARKER \
            --name_se STDERR \
            --pc 1 \
            --no_std_names
        """

rule identify_lead_variants:
    input:
        "data/processed/metal/{phenotype}-ALL_metal-1.txt"
    output:
        lead_vars = "data/processed/lead-variants/metal/{phenotype}-ALL_lead_variants.txt",
        fig = "results/mqq-plots/metal/{phenotype}-ALL_mqq_plot.png"
    conda:
        "envs/gwaslab.yaml"
    script:
        "scripts/meta-analysis/lead_vars_mqq.py"

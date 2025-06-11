configfile: "config/config.yaml"# Load metadata

df = pd.read_table("config/studies.txt")
df["id"] = df.apply(lambda row: f"{row.study}_{row.phenotype}_{row.population}", axis=1)
metadata = {row.id: row._asdict() for row in df.itertuples(index=False)}
ids = list(metadata.keys())
id_to_phenotype = {row.id: row.phenotype for row in df.itertuples(index=False)}

rule filter_maf:
    input:
        lambda wildcards: f"results/gwas-inspector/{id_to_phenotype[wildcards.id].lower()}-ldsc/QC_{wildcards.id}.txt"
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
        rules.filter_maf.output
    output:
        metal = expand("data/processed/metal/config/{phenotype}-{population}.txt", 
            phenotype=config['create_metal_configs']['phenotypes'],
            population=config['create_metal_configs']['population']),
        mrmega = expand("data/processed/mr-mega/config/{phenotype}-{population}.in", 
            phenotype=config['create_metal_configs']['phenotypes'],
            population=config['create_metal_configs']['population'])
    params:
        utils = config["create_metal_config"]["utils"],
        data_dir = "data/processed/sumstats-maf-filtered",
        metal = "data/processed/metal/config",
        mrmega = "data/processed/mr-mega/config"
    script:
        "scripts/meta-analysis/create_metal_config.R"

rule run_metal:
    input:
        "data/processed/metal/config/{config}.txt"
    output:
        "data/processed/metal/{config}_metal-1.txt"
    shell:
        """
        module load metal/2020.05
        metal {input}
        """

rule run_mrmega:
    input:
        "data/processed/mr-mega/config/{config}.in"
    output:
        "data/processed/mr-mega/{config}.result"
    shell:
        """
        module load MR-MEGA/0.2

        NUM_DATASETS=$(wc -l < {input})
        NUM_PCS=$((NUM_DATASETS - 3))

        if [ "$NUM_PCS" -lt 1 ]; then
            NUM_PCS=1
        fi

        MR-MEGA --filelist {input} --out {output} \
            --qt \
            --name_pos POS \
            --name_chr CHR \
            --name_n N \
            --name_beta BETA \
            --name_eaf EAF \
            --name_nea NEA \
            --name_ea EA \
            --name_marker SNPID \
            --name_se SE \
            --pc $NUM_PCS \
            --no_std_names
        """
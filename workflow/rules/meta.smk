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
        metal = metal_config
    params:
        utils = config["create_meta_analysis_configs"]["utils"],
        data_dir = "data/processed/sumstats-maf-filtered",
        metal = "data/processed/metal"
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

rule process_gwas_sumstats:
    input:
        "data/raw/gwas-summary-statistics/{study}_{phenotype}_{population}_{genome}.txt"
    output:
        "data/processed/gwas-summary-statistics/{study}_{phenotype}_{population}_GrCh38.txt"
    params:
        study = wildcards.study,
        phenotype = wildcards.phenotype,
        genome = wildcards.genome,
        utils = config["process_gwas"]["utils"]
    script:
        "workflow/rules/scripts/input/process_gwas.R"
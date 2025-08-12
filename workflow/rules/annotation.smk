rule run_hyprcoloc_eqtl:
    input:
        expand("data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt", phenotype=config["metal_format"]["phenotypes"])
    output:
        colocalization = "data/processed/variant-annotation/eqtl/HyPrColoc_eqtl_results.txt",
    conda:
        "envs/hyprcoloc.yaml"
    params:
        locus = config["hyprcoloc"]["locus"]
    resources:
        mem_mb = 200000
    script:
        "scripts/variant-annotation/hyprcoloc_eqtl.R"
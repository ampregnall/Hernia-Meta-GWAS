rule run_magma:
    input:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt"
    output:
        infile = "data/processed/magma/{phenotype}/{phenotype}.txt",
        out = "data/processed/magma/{phenotype}/{phenotype}.genes.out",
        raw = "data/processed/magma/{phenotype}/{phenotype}.genes.raw"
    params:
        bfile = "/project/voltron/Resources/1000-Genomes-Phase3/1kg_all_phase3/all_phase3",
        gene_file = "/project/voltron/Resources/MAGMA/NCBI37/NCBI37.genes.annot",
        out_prefix = lambda wc: f"data/processed/magma/{wc.phenotype}",
        snp_col = "SNP",
        p_col = "P",
        n_col = "N"      
    shell:
        """
        module load magma/1.07a 

        awk 'FNR==1 {{ $15 = "SNP"; print $9, $10, $15 }} FNR > 1 && $15 != "" {{ print $9, $10, $15 }}' {input.sumstats} > {output.infile} 
        
        magma --bfile {params.bfile} \
              --gene-annot {params.gene_file} \
              --pval {output.infile} use={params.snp_col},{params.p_col} ncol={params.n_col} \
              --gene-model snp-wise=mean \
              --out {params.out_prefix}
        """

rule run_pops:  
    input:
        magma_out = "data/processed/magma/{phenotype}/{phenotype}.genes.out",
        magma_raw = "data/processed/magma/{phenotype}/{phenotype}.genes.raw"
    output:
        preds = "data/processed/pops/{phenotype}/{phenotype}.preds",
        coefs = "data/processed/pops/{phenotype}/{phenotype}.coefs",
        marginals = "data/processed/pops/{phenotype}/{phenotype}.marginals"
    params:
        out_prefix = lambda wc: f"data/processed/pops/{wc.phenotype}",
        magma_prefix = lambda wc: f"data/processed/magma/{wc.phenotype}",
        gene_annot_path = "data/meta/pops/utils/gene_annot.txt",
        control_features_path = "data/meta/pops/utils/features_control.txt",
        feature_mat_prefix = "data/meta/pops/features-munged"
    conda:
        "envs/pops.yaml"
    shell:
        """
        python workflow/rules/scripts/annotation/pops.py \
            --gene_annot_path {params.gene_annot_path} \
            --feature_mat_prefix example/data/features_munged/pops_features \
            --num_feature_chunks 115 \
            --magma_prefix {params.magma_prefix} \
            --control_features_path {params.control_features_path} \
            --out_prefix {params.out_prefix}
        """

# TODO: Create rule to run HyPrColo
# TODO Create rule to run open chromatin annotation
# TODO: Create rule to run most severe consequence annotation

# rule run_hyprcoloc_eqtl:
#     input:
#         expand("data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt", phenotype=config["metal_format"]["phenotypes"])
#     output:
#         colocalization = "data/processed/variant-annotation/eqtl/HyPrColoc_eqtl_results.txt",
#     conda:
#         "envs/hyprcoloc.yaml"
#     params:
#         locus = config["hyprcoloc"]["locus"]
#     resources:
#         mem_mb = 200000
#     script:
#         "scripts/variant-annotation/hyprcoloc_eqtl.R"
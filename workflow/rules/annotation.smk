rule run_magma:
    input:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt"
    output:
        infile = "data/processed/magma/{phenotype}/{phenotype}.txt",
        out = "data/processed/magma/{phenotype}/{phenotype}.genes.out",
        raw = "data/processed/magma/{phenotype}/{phenotype}.genes.raw"
    params:
        bfile = "data/meta/magma/1000G_all_hg38",
        gene_file = "data/meta/magma/GRCh38.genes.ensg.genes.annot",
        out_prefix = lambda wc: f"data/processed/magma/{wc.phenotype}/{wc.phenotype}",
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
        preds = "data/processed/pops/{phenotype}.preds",
        coefs = "data/processed/pops/{phenotype}.coefs",
        marginals = "data/processed/pops/{phenotype}.marginals"
    params:
        out_prefix = lambda wc: f"data/processed/pops/{wc.phenotype}",
        magma_prefix = lambda wc: f"data/processed/magma/{wc.phenotype}/{wc.phenotype}",
        gene_annot_path = "data/meta/pops/utils/gene_annot.txt",
        control_features_path = "data/meta/pops/utils/features_control.txt",
        feature_mat_prefix = "data/meta/pops/features-munged/features-munged"
    conda:
        "envs/pops.yaml"
    shell:
        """
        python workflow/rules/scripts/variant-annotation/pops.py \
            --gene_annot_path {params.gene_annot_path} \
            --feature_mat_prefix {params.feature_mat_prefix} \
            --num_feature_chunks 115 \
            --magma_prefix {params.magma_prefix} \
            --control_features_path {params.control_features_path} \
            --out_prefix {params.out_prefix}
        """

rule run_hyprcoloc_eqtl:
    input:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt",
        lead = "data/processed/lead-variants/metal/{phenotype}-ALL_lead_variants.txt"
    output:
        eqtl = "data/processed/variant-annotation/eqtl/{phenotype}_eqtl_colocalization.txt",
    conda:
        "envs/hyprcoloc.yaml"
    params:
        locus = config["hyprcoloc"]["locus"],
        eqtl = "data/meta/eqtl/QTD000216.all.tsv",
        min_overlap = 50
    script:
         "scripts/variant-annotation/hyprcoloc_eqtl.R"

# TODO Create rule to run open chromatin annotation
# TODO: Create rule to run most severe consequence annotation
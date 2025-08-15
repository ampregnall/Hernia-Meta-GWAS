rule run_magma:
    input:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt"
    output:
        infile = "data/processed/magma/{phenotype}.txt",
        out = "data/processed/magma/{phenotype}.genes.out",
        raw = "data/processed/magma/{phenotype}.genes.raw"
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

# TODO: Create rule to run HyPrColo
# TODO: Create rule to run PoPs
# TODO Create rule to run open chromatin annotation
# TODO: Create rule to run most severe consequence annotation
# TODO Create rule to merge all results???

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
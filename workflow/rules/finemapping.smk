rule identify_lead_variants_metal:
    input:
        "data/processed/metal/{phenotype}_metal-1.txt"
    output:
        lead_vars = "data/processed/lead-variants/metal/{phenotype}_lead_variants.txt",
        fig = "results/mqq-plots/metal/{phenotype}_mqq_plot.png",
    params:
        rsid_tsv = config["ldsc_correction"]["rsid_tsv"],
        rsid_vcf = config["ldsc_correction"]["rsid_vcf"]
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

rule format_metal:
    input:
        "data/processed/metal/{phenotype}-ALL_metal-1.txt"
    output:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt"
    params:
        rsid_tsv = config["ldsc_correction"]["rsid_tsv"],
        rsid_vcf = config["ldsc_correction"]["rsid_vcf"]
    conda:
        "envs/gwaslab.yaml"
    threads: 4
    script:
        "scripts/finemapping/format_metal.py"

rule run_hyprcoloc:
    input:
        expand("data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt", phenotype=config["metal_format"]["phenotypes"])
    output:
        colocalization = "data/processed/finemapping/HyPrColoc_results.txt",
        credset = "data/processed/finemapping/HyPrColoc_finemapping.txt"
    conda:
        "envs/hyprcoloc.yaml"
    params:
        creds = config["hyprcoloc"]["creds"],
        locus = config["hyprcoloc"]["locus"]
    resources:
        mem_mb = 120000
    script:
        "scripts/finemapping/run_hyprcoloc.R"

rule bayesian_finemapping:
    input:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt",
        lead = "data/processed/lead-variants/metal/{phenotype}-ALL_lead_variants.txt"
    output:
        credset = "data/processed/finemapping/{phenotype}_bayesian_finemapping.txt"
    params:
        creds = config["bayesian_finemapping"]["creds"],
        locus = config["bayesian_finemapping"]["locus"]
    conda:
        "envs/common.yaml"
    script:
        "scripts/finemapping/bayesian_finemapping.R"
    
checkpoint identify_finemapping_loci:
    input:
        lead = "data/processed/lead-variants/metal/{phenotype}-ALL_lead_variants.txt"
    output:
        directory("data/processed/rsparsepro/{phenotype}/loci")
    run:
        df = pd.read_csv(input.lead, sep="\t", dtype={"CHR": str})
        outdir = Path(output[0])
        outdir.mkdir(parents=True, exist_ok=True)

        for _, row in df.iterrows():
            locus_name = f"{row['CHR']}_{int(row['POS'])}"
            # row is a Series — wrap it so to_csv writes a header + one row
            pd.DataFrame([row]).to_csv(outdir / f"{locus_name}.txt",
                                       sep="\t", index=False)

rule extract_finemapping_loci_snps:
    input:
        locus = "data/processed/rsparsepro/{phenotype}/loci/{locus}.txt",
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt"
    output:
        snps = "data/processed/rsparsepro/{phenotype}/snps/{locus}.snplist",
        a1 = "data/processed/rsparsepro/{phenotype}/snps/{locus}.a1"
    params:
        window = 500000
    run:
        # Create output directory
        Path(output.snps).parent.mkdir(parents=True, exist_ok=True)
        
        # Load data
        df = pd.read_csv(input.sumstats, sep="\t", dtype={"CHR": str, "POS" : int})
        lead = pd.read_csv(input.locus, sep="\t", dtype={"CHR": str, "POS" : int})

        # Extract locus information
        chrom = str(lead.iloc[0]["CHR"])
        pos   = int(lead.iloc[0]["POS"])
        start = max(0, pos - params.window)
        end   = pos + params.window

        snps = df[(df["CHR"] == chrom) & (df["POS"].between(start, end)) & (df["rsID"] != "") & (df["rsID"].notna())]
        snps = snps.drop_duplicates(subset=['rsID'], keep=False)
        snps[["rsID", "EA"]].to_csv(output.a1, sep="\t", index=False, header=False)
        snps[["rsID"]].to_csv(output.snps, sep="\t", index=False, header=False)

rule plink_ld:
    input:
        snps = "data/processed/rsparsepro/{phenotype}/snps/{locus}.snplist",
        a1 = "data/processed/rsparsepro/{phenotype}/snps/{locus}.a1"
    output:
        ld = "data/processed/rsparsepro/{phenotype}/ld/{locus}.ld",
        snplist = "data/processed/rsparsepro/{phenotype}/ld/{locus}.snplist",
        log = "data/processed/rsparsepro/{phenotype}/ld/{locus}.log"    
    shell:
        """
        module load plink/1.9-20210416
        plink --bfile /project/voltron/Resources/1000-Genomes-Phase3/1kg_all_phase3/all_phase3 \
            --extract {input.snps} \
            --a1-allele {input.a1} 2 1 \
            --r square \
            --write-snplist \
            --out data/processed/rsparsepro/{wildcards.phenotype}/ld/{wildcards.locus}
        """

rule extract_rsparsepro_input:
    input:
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt",
        snps = "data/processed/rsparsepro/{phenotype}/ld/{locus}.snplist"
    output:
        zfile = "data/processed/rsparsepro/{phenotype}/z-files/{locus}.txt"
    run:
        df = pd.read_csv(input.sumstats, sep="\t", dtype={"CHR": str, "POS": int})
        snps = pd.read_csv(input.snps, header=None)[0].tolist()  
        df = df[df['rsID'].isin(snps)]
        df = df.set_index('rsID')  
        df = df.reindex(snps) 
        df[["Z"]].to_csv(output.zfile, sep="\t", index=True, index_label="RSID")  

rule run_rsparsepro:
    input:
        zfile = "data/processed/rsparsepro/{phenotype}/z-files/{locus}.txt",
        ld = "data/processed/rsparsepro/{phenotype}/ld/{locus}.ld"
    output:
        finemapping = "data/processed/rsparsepro/{phenotype}/finemapping/{locus}.rsparsepro.txt"
    conda:
        "envs/rsparsepro.yaml"
    resources:
        mem_mb = 4000
    shell:
        """
        python workflow/rules/scripts/finemapping/rsparsepro_ld.py \
            --z {input.zfile} \
            --ld {input.ld} \
            --save "data/processed/rsparsepro/{wildcards.phenotype}/finemapping/{wildcards.locus}"
        """

def rsparsepro_files(wc):
    ck = checkpoints.identify_finemapping_loci.get(phenotype=wc.phenotype)
    # ck.output[0] is the loci directory
    import os
    loci = glob_wildcards(os.path.join(ck.output[0], "{locus}.txt")).locus
    return expand(
        "data/processed/rsparsepro/{phenotype}/finemapping/{locus}.rsparsepro.txt",
        phenotype=wc.phenotype, locus=loci
    )

rule collect_rsparsepro:
    input: rsparsepro_files
    output:
        finemapping = "data/processed/rsparsepro/out/{phenotype}_rsparsepro_finemapping.txt"
    shell:
        """
        awk 'FNR==1 && NR!=1 {{next}} NR==1 {{print; next}} $5 != 0 {{print}}' {input} > {output.finemapping}
        """

rule annotate_rsparsepro:
    input:
        rsparsepro = "data/processed/rsparsepro/out/{phenotype}_rsparsepro_finemapping.txt",
        sumstats = "data/processed/metal-formatted/{phenotype}-ALL_metal_formatted.txt",
        lead = "data/processed/lead-variants/metal/{phenotype}-ALL_lead_variants.txt"
    output:
        finemapping = "data/processed/finemapping/{phenotype}_rsparsepro_finemapping.txt"
    params:
        locus = config["bayesian_finemapping"]["locus"]
    conda:
        "envs/common.yaml"
    script:
        "scripts/finemapping/annotate_rsparsepro.R"
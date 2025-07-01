import gwaslab as gl

ss = gl.Sumstats(snakemake.input[0],
             snpid="MarkerName",
             chrom="Chromosome",
             pos="Position",
             ea="Allele1",
             nea="Allele2",
             eaf="Freq1",
             beta="Effect",
             se="StdErr",
             p="P-value",
             ncase="NCASE",
             ncontrol="NCONTROL",
             n="SAMPLESIZE",
             sep="\t", 
             build="38",
             verbose=True)

# Perform basic check on sumstats
ss.basic_check()

# Extract lead vars
lead_vars = ss.get_lead(anno=True)
lead_vars.to_csv(snakemake.output["lead_vars"], sep="\t", index=False)

var_list = lead_vars.sort_values(by="P").head(10)

# Create MQQ plot with lead variant annotations
ss.plot_mqq(anno="GENENAME", highlight=list(var_list.SNPID), anno_set=list(var_list.SNPID), 
            save=snakemake.output["fig"], save_args={"dpi":300,"facecolor":"white"})

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

# Calculate Z score
ss.data['Z'] = ss.data['BETA'] / ss.data['SE']

# Assign rsID 
ss.assign_rsid(ref_rsid_tsv = snakemake.params["rsid_tsv"],
               ref_rsid_vcf = snakemake.params["rsid_vcf"],
               chr_dict = gl.get_number_to_NC(build="38"), 
               n_cores = snakemake.threads)

# Save formatted sumstats for later steps
ss.data.to_csv(snakemake.output["sumstats"], sep="\t", index=False)

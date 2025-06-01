import gwaslab as gl
import pandas as pd
import numpy as np
from scipy.stats import norm

ss = gl.Sumstats(snakemake.input[0]),
             snpid="MARKER",
             chrom="CHROM",
             pos="POS",
             ea="EFFECT_ALLELE",
             nea="OTHER_ALLELE",
             eaf="EAF",
             beta="BETA",
             se="SE",
             p="P",
             ncase="N_CASE",
             ncontrol="N_CONTROL",
             n="N",
             build="38",      
             sep="\t", 
             verbose=True)

ss.data['SNPID'] = ss.data[['CHR', "POS", "NEA", "EA"]].astype(str).agg(':'.join, axis=1)

ss.assign_rsid(ref_rsid_tsv = "1kg_dbsnp151_hg38_auto.txt.gz",
               chr_dict = gl.get_number_to_NC(build="38"), n_cores = snakemake.threads)

ss_fil = ss.filter_hapmap3(inplace=False)

ss_fil.estimate_h2_by_ldsc(ref_ld = "AMR/UKBB.AMR",  w_ld = "AMR/UKBB.AMR")


if np.float64(ss_fil.ldsc_h2['Intercept'][0]) > 1:
    ss.data['SE'] = ss.data['SE'] * np.sqrt(np.float64(ss_fil.ldsc_h2['Intercept'][0]))
    ss.data['Z'] = ss.data['BETA'] / ss.data['SE']
    ss.data['P'] = 2 * norm.sf(abs(ss.data['Z']))

# Save results
ss.data.to_csv(snakemake.output[0], sep="\t", index=False)
ss_fil.ldsc_h2(snakemake.output[1], sep="\t", index=False)

import gwaslab as gl
import pandas as pd
import numpy as np
from scipy.stats import norm

ss = gl.Sumstats(
    snakemake.input[0],
    snpid="ID",
    rsid="SNP",
    chrom="CHR",
    pos="BP",
    ea="A2",
    nea="A1",
    eaf="FRQ",
    beta="BETA",
    se="SE",
    p="P",
    ncase="N_CAS",
    ncontrol="N_CON",
    n="N",
    build="38",
    sep="\t",
    verbose=True,
)

# ss.assign_rsid(ref_rsid_tsv = snakemake.params["rsid_tsv"],
#                ref_rsid_vcf = snakemake.params["rsid_vcf"],
#                chr_dict = gl.get_number_to_NC(build="38"),
#                n_cores = snakemake.threads)

ss_fil = ss.filter_hapmap3(inplace=False)

ss_fil.estimate_h2_by_ldsc(
    ref_ld=snakemake.params["ldsc_scores"], w_ld=snakemake.params["ldsc_scores"]
)

if np.float64(ss_fil.ldsc_h2["Intercept"][0]) > 1:
    ss.data["SE"] = ss.data["SE"] * np.sqrt(np.float64(ss_fil.ldsc_h2["Intercept"][0]))
    ss.data["Z"] = ss.data["BETA"] / ss.data["SE"]
    ss.data["P"] = 2 * norm.sf(abs(ss.data["Z"]))

# Save results
ss.data.to_csv(snakemake.output["sumstats"], sep="\t", index=False)
ss_fil.ldsc_h2.to_csv(snakemake.output["ldsc"], sep="\t", index=False)

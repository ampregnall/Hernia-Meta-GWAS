import sys

import gwaslab as gl
import numpy as np
import pandas as pd
from scipy.stats import norm

log = open(snakemake.log[0], "w")
sys.stdout = sys.stderr = log

ss = gl.Sumstats(
    snakemake.input["sumstats"],
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

# Convert chromosome labels to numbers (X -> 23). gwaslab keeps CHR as a string on
# load; only fix_chr() converts it.
ss.fix_chr()

# Estimate the LDSC intercept on HapMap3 variants. The LD scores are autosomal, so
# chrX does not contribute to the estimate; the correction is applied to all variants.
ss_fil = ss.filter_hapmap3(inplace=False)

ss_fil.estimate_h2_by_ldsc(
    ref_ld=snakemake.params["ldsc_scores"], w_ld=snakemake.params["ldsc_scores"]
)

intercept = np.float64(ss_fil.ldsc_h2["Intercept"][0])
if intercept > 1:
    ss.data["SE"] = ss.data["SE"] * np.sqrt(intercept)
    ss.data["Z"] = ss.data["BETA"] / ss.data["SE"]
    ss.data["P"] = 2 * norm.sf(abs(ss.data["Z"]))

# Minor allele count filter, applied after the intercept is estimated so that rare
# variants do not change the HapMap3 SNP set used for LDSC
min_case_mac = snakemake.params["min_case_mac"]
maf = np.minimum(ss.data["EAF"], 1 - ss.data["EAF"])
ss.data = ss.data[maf > min_case_mac / (2 * ss.data["N_CASE"])]

counts = pd.read_table(snakemake.input["counts"])
mac_row = pd.DataFrame(
    {
        "id": [snakemake.wildcards["id"]],
        "step": ["mac_filter"],
        "n_variants": [len(ss.data)],
        "n_chrX": [int((ss.data["CHR"] == 23).sum())],
        "n_removed": [counts["n_variants"].iloc[-1] - len(ss.data)],
    }
)
pd.concat([counts, mac_row]).to_csv(snakemake.output["counts"], sep="\t", index=False)

# Write the columns used by the meta-analysis with fixed types, so every study has
# the same parquet schema. CHR is numeric with X coded as 23.
schema = {
    "SNPID": "string",
    "rsID": "string",
    "CHR": "int64",
    "POS": "int64",
    "EA": "string",
    "NEA": "string",
    "EAF": "float64",
    "BETA": "float64",
    "SE": "float64",
    "P": "float64",
    "N_CASE": "int64",
    "N_CONTROL": "int64",
    "N": "int64",
}
out = ss.data[list(schema)].astype(schema)
out.to_parquet(snakemake.output["sumstats"], index=False)
ss_fil.ldsc_h2.to_csv(snakemake.output["ldsc"], sep="\t", index=False)

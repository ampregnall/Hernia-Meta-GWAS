import sys

import gwaslab as gl
import pandas as pd

log = open(snakemake.log[0], "w")
sys.stdout = sys.stderr = log

df = pd.read_parquet(snakemake.input[0])

ss = gl.Sumstats(
    df,
    snpid="SNPID",
    rsid="rsID",
    chrom="CHR",
    pos="POS",
    ea="EA",
    nea="NEA",
    eaf="EAF",
    beta="BETA",
    se="SE",
    p="P",
    z="Z",
    mlog10p="MLOG10P",
    ncase="N_CASE",
    ncontrol="N_CONTROL",
    n="N",
    other=["N_CONTRIBUTIONS", "Q", "Q_DF", "Q_PVAL", "I2"],
    build="38",
    verbose=True,
)

# Perform basic check on sumstats. rsIDs are carried through the meta-analysis from
# harmonization, so they are not reassigned here.
ss.basic_check()

# Save formatted sumstats for later steps
ss.data.to_csv(snakemake.output["sumstats"], sep="\t", index=False)

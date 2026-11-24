import sys

import gwaslab as gl

log = open(snakemake.log[0], "w")
sys.stdout = sys.stderr = log

# Formatted meta-analysis output (already passed basic_check in format_meta)
ss = gl.Sumstats(
    snakemake.input[0],
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
    ncase="N_CASE",
    ncontrol="N_CONTROL",
    n="N",
    sep="\t",
    build="38",
    verbose=True,
)

# Extract lead vars
lead_vars = ss.get_lead(anno=True)
lead_vars.to_csv(snakemake.output["lead_vars"], sep="\t", index=False)

var_list = lead_vars.sort_values(by="P").head(10)

# Create MQQ plot with lead variant annotations
ss.plot_mqq(
    anno="GENENAME",
    highlight=list(var_list.SNPID),
    anno_set=list(var_list.SNPID),
    save=snakemake.output["fig"],
    save_args={"dpi": 300, "facecolor": "white"},
)

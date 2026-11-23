import argparse
import gwaslab as gl
import pandas as pd
import matplotlib.pyplot  as plt

parser = argparse.ArgumentParser(description="Extract lead variants and create MQQ plot from MR-MEGA results")
parser.add_argument("--input", required=True, help="MR-MEGA .result file")
parser.add_argument("--lead-vars", required=True, help="Output path for lead variants")
parser.add_argument("--fig", required=True, help="Output path for MQQ plot")
args = parser.parse_args()

col_names = ["SNPID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE",
             "P", "N", "N_CASE", "N_CONTROL", "STATUS", "LOCATION", "GENE"]

df = pd.read_csv(args.input, sep="\t")
df = df.dropna(subset=["P-value_association"])

if len(df) == 0:
    print("Sample size not large enough to perform meta-analysis for input: ")
    print(args.input)

    lead_vars = pd.DataFrame(columns=col_names)
    lead_vars.to_csv(args.lead_vars, sep="\t", index=False)

    plt.figure(figsize=(8,6))
    plt.text(0.5, 0.5, "Insufficient sample size for meta-analysis", ha="center", va="center", fontsize=16)
    plt.title(args.input)
    plt.axis("off")
    plt.savefig(args.fig, dpi=300, facecolor="white")
    plt.close()
else:
    ss = gl.Sumstats(df,
                snpid="MarkerName",
                chrom="Chromosome",
                pos="Position",
                ea="EA",
                nea="NEA",
                eaf="EAF",
                beta="beta_0",
                se="se_0",
                p="P-value_association",
                ncase="Nsample",
                sep="\t",
                build="38",
                verbose=True)

    # Perform basic check on sumstats
    ss.basic_check()

    # Extract lead vars
    lead_vars = ss.get_lead(anno=True)
    lead_vars.to_csv(args.lead_vars, sep="\t", index=False)

    var_list = lead_vars.sort_values(by="P").head(10)

    # Create MQQ plot with lead variant annotations
    ss.plot_mqq(anno="GENENAME", highlight=list(var_list.SNPID), anno_set=list(var_list.SNPID),
                save=args.fig, save_args={"dpi":300,"facecolor":"white"})

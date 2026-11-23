from pathlib import Path
import gwaslab as gl
import pandas as pd
import polars as pl

LOCUS_SIZE = 1_000_000

df = pl.read_csv(snakemake.input.sumstats, separator="\t")
lead = pl.read_csv(snakemake.input.lead, separator="\t")

if lead.height == 0:
    Path(snakemake.output.vep_input).write_text("")
    raise SystemExit(0)

lead = lead.select(["CHR", "POS"]).with_columns(
    (pl.col("POS") + (LOCUS_SIZE / 2)).alias("END"),
    (pl.col("POS") - (LOCUS_SIZE / 2)).alias("START"),
)

df = df.join(lead, on="CHR", how="inner", suffix="_LEAD", coalesce=False)

df = df.filter(
    (
        (pl.col("CHR") == pl.col("CHR_LEAD"))
        & (pl.col("POS") <= pl.col("END"))
        & (pl.col("POS") >= pl.col("START"))
    )
    | (pl.col("P") < 0.00000005)
)

if df.height == 0:
    Path(snakemake.output.vep_input).write_text("")
    raise SystemExit(0)

df = df.select(pl.all().exclude(["CHR_LEAD", "POS_LEAD", "START", "END"]))

sumstats = gl.Sumstats(
    pd.DataFrame(df, columns=df.columns),
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
    n="N",
    ncase="N_CASE",
    ncontrol="N_CONTROL",
    status="STATUS",
    z="Z",
)

sumstats.to_format(snakemake.params.vep_prefix, fmt="vep")

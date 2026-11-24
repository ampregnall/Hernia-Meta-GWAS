#!/usr/bin/env bash
# Count chrX variants in GWAS summary statistics files.
# Works on plain or gzipped, whitespace- or tab-delimited files.
#
# Usage:  bash count_chrx.sh path/to/sumstats/*  > chrx_counts.tsv
#
# Output columns:
#   file      basename of the input
#   total     number of variants (data rows)
#   chrX      variants coded X / 23 / 25 (PLINK PAR code), with or without a "chr" prefix
#   pct_X     chrX as a percent of total
#   min_pos   smallest chrX position (PAR1 sits below ~2.7 Mb)
#   max_pos   largest chrX position (PAR2 sits above ~155 Mb in hg38)
#
# A file with no recognizable chromosome column is reported as "no chr column"
# (for example, when the only locus info is an ID like 1:12345:A:G).

printf "file\ttotal\tchrX\tpct_X\tmin_pos\tmax_pos\n"

for f in "$@"; do
  gzip -cdf "$f" | awk -v f="$(basename "$f")" '
    BEGIN { FS = "[ \t]+"; OFS = "\t" }
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        h = toupper($i); gsub(/^#/, "", h)
        if (h == "CHR" || h == "CHROM" || h == "CHROMOSOME") c = i
        if (h == "POS" || h == "BP" || h == "POSITION" || h == "GENPOS" || h == "BASE_PAIR_LOCATION") p = i
      }
      if (!c) { print f, "NA", "NA", "no chr column", ".", "."; exit }
      next
    }
    {
      n++
      ch = toupper($c); sub(/^CHR/, "", ch)
      if (ch == "X" || ch == "23" || ch == "25") {
        x++
        if (p) {
          pos = $p + 0
          if (!seen || pos < mn) mn = pos
          if (!seen || pos > mx) mx = pos
          seen = 1
        }
      }
    }
    END {
      if (c) print f, n + 0, x + 0, (n ? sprintf("%.2f", 100 * x / n) : "NA"), (seen ? mn : "."), (seen ? mx : ".")
    }'
done
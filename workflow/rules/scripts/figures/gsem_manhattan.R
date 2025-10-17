source("workflow/rules/scripts/figures/gg_manhattan.R")
library(tidyverse)
library(tidygenomics)
library(ggfastman)

sumstats <- vroom::vroom(
  file = "data/processed/genomic-sem/common-factor-gwas/Hernia_GSEM_common_factor.txt.gz",
  delim = "\t"
)

# Load lead variants
lead_variant_files <- c(
  "data/processed/lead-variants/metal/Ventral-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Inguinal-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Umbilical-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Femoral-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Diaphragmatic-ALL_lead_variants.txt"
)

lead_variants <- vroom::vroom(lead_variant_files, id = "ID") |>
  dplyr::mutate(ID = stringr::str_split_i(basename(ID), pattern = "-", i = 1))

gsem_variants <- vroom::vroom(
  "data/processed/lead-variants/genomic-sem/Hernia_GSEM_common_factor_lead_variants.txt"
) |>
  dplyr::mutate(
    ID = "Common Factor GWAS"
  ) # Hack for function

lead_variants <- dplyr::bind_rows(lead_variants, gsem_variants) |>
  dplyr::mutate(START = POS - 500000, END = POS + 500000) |>
  tidygenomics::genome_cluster(
    by = c("CHR", "START", "END"),
    max_distance = 0
  )

lead_variants <- lead_variants |>
  dplyr::left_join(
    lead_variants |>
      dplyr::group_by(cluster_id) |>
      dplyr::summarise(count = dplyr::n()) |>
      dplyr::mutate(
        shared_locus = dplyr::case_when(count > 1 ~ TRUE, .default = FALSE)
      ) |>
      dplyr::select(!count),
    by = "cluster_id"
  )

annotation_df <- lead_variants |>
  dplyr::filter(ID == "Common Factor GWAS")

plt <- gg_manhattan_df(
  sumstats_df = sumstats,
  annotation_df = annotation_df,
  phenotype = "Common Factor GWAS",
  use_broken_axis = FALSE
)

#
gsem_loci <- lead_variants |>
  dplyr::filter(ID == "Common Factor GWAS") |>
  dplyr::pull(cluster_id)

gsem_variants <- lead_variants |>
  dplyr::filter(cluster_id %in% gsem_loci) |>
  dplyr::mutate(log10p = -log10(P)) |>
  dplyr::mutate(
    ID = factor(
      ID,
      levels = c(
        "Common Factor GWAS",
        "Ventral",
        "Inguinal",
        "Umbilical",
        "Femoral",
        "Diaphragmatic"
      )
    )
  )

gsem_variants <- gsem_variants |>
  dplyr::group_by(cluster_id, GENE) |>
  dplyr::summarise(
    count = n(),
    min_p = min(P), # lowest p-value for that gene in the cluster
    .groups = "drop_last"
  ) |>
  arrange(cluster_id, desc(count), min_p) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(cluster_id, LABEL = GENE) %>%
  right_join(gsem_variants, by = "cluster_id")

geneplt <- gsem_variants |>
  ggplot(aes(x = LABEL, y = forcats::fct_rev(ID), color = ID)) +
  geom_point(aes(size = log10p, shape = factor(shared_locus)), alpha = 0.8) +
  scale_shape_manual(values = c(8, 19)) +
  facet_grid(~CHR, scales = "free_x", space = "free_x") +
  ggsci::scale_color_npg() +
  scale_size_continuous(range = c(1, 5)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.spacing.x = unit(1.5, "mm"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(0, 5, 5, 5)
  )

fig <- cowplot::plot_grid(plt, geneplt, nrow = 2, align = "v", axis = "lr")

ggsave(
  filename = "figures/gsem_manhattan.pdf",
  plot = fig,
  width = 15,
  height = 1.75 * 5,
  units = "in"
)

library(tidyverse)
library(tidygenomics)
library(ggfastman)
source("workflow/rules/scripts/figures/gg_manhattan.R")


# Load summary statistics
ventral <- vroom::vroom(
  "data/processed/metal-formatted/Ventral-ALL_metal_formatted.txt"
)
inguinal <- vroom::vroom(
  "data/processed/metal-formatted/Inguinal-ALL_metal_formatted.txt"
)
umbilical <- vroom::vroom(
  "data/processed/metal-formatted/Umbilical-ALL_metal_formatted.txt"
)
femoral <- vroom::vroom(
  "data/processed/metal-formatted/Femoral-ALL_metal_formatted.txt"
)
diaphragmatic <- vroom::vroom(
  "data/processed/metal-formatted/Diaphragmatic-ALL_metal_formatted.txt"
)

# Load lead variants
lead_files <- c(
  "data/processed/lead-variants/metal/Ventral-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Inguinal-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Umbilical-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Femoral-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Diaphragmatic-ALL_lead_variants.txt"
)

lead_vars <- vroom::vroom(lead_files, id = "ID") |>
  dplyr::bind_rows() |>
  dplyr::mutate(START = POS - 500000, END = POS + 500000) |>
  tidygenomics::genome_cluster(
    by = c("CHR", "START", "END"),
    max_distance = 0
  )

lead_vars <- lead_vars |>
  dplyr::left_join(
    lead_vars |>
      dplyr::group_by(cluster_id) |>
      dplyr::summarise(count = dplyr::n()) |>
      dplyr::mutate(
        shared_locus = dplyr::case_when(count > 1 ~ TRUE, .default = FALSE)
      ) |>
      dplyr::select(!count),
    by = "cluster_id"
  )

ventral_plt <- gg_manhattan_df(
  sumstats_df = ventral,
  annotation_df = lead_vars,
  phenotype = "Ventral",
  use_broken_axis = FALSE
)
inguinal_plt <- gg_manhattan_df(
  sumstats_df = inguinal,
  annotation_df = lead_vars,
  phenotype = "Inguinal",
  use_broken_axis = FALSE
)
umbilical_plt <- gg_manhattan_df(
  sumstats_df = umbilical,
  annotation_df = lead_vars,
  phenotype = "Umbilical",
  use_broken_axis = FALSE
)
femoral_plt <- gg_manhattan_df(
  sumstats_df = femoral,
  annotation_df = lead_vars,
  phenotype = "Femoral",
  use_broken_axis = FALSE
)
diaphragmatic_plt <- gg_manhattan_df(
  sumstats_df = diaphragmatic,
  annotation_df = lead_vars,
  phenotype = "Diaphragmatic",
  use_broken_axis = FALSE
)

fig <- cowplot::plot_grid(
  plotlist = c(
    ventral_plt,
    inguinal_plt,
    umbilical_plt,
    femoral_plt,
    diaphragmatic_plt
  ),
  align = "hv",
  ncol = 1
)

ggsave(
  filename = "figures/manhattan.pdf",
  plot = fig,
  width = 7.5,
  height = 1.75 * 5,
  units = "in"
)

# Load lead variants
metal_files <- c(
  "data/processed/lead-variants/metal/Ventral-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Inguinal-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Umbilical-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Femoral-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Diaphragmatic-ALL_lead_variants.txt"
)

mega_files <- c(
  "data/processed/lead-variants/mr-mega/Ventral-ALL_lead_variants.txt",
  "data/processed/lead-variants/mr-mega/Inguinal-ALL_lead_variants.txt",
  "data/processed/lead-variants/mr-mega/Umbilical-ALL_lead_variants.txt",
  "data/processed/lead-variants/mr-mega/Diaphragmatic-ALL_lead_variants.txt"
)

metal_vars <- vroom::vroom(metal_files, id = "ID") |>
  dplyr::mutate(
    ID = glue::glue(
      "{stringr::str_split_i(basename(ID), pattern = '-', i = 1)} (METAL)"
    )
  ) |>
  dplyr::select(ID, CHR, POS, GENE)

mega_vars <- vroom::vroom(mega_files, id = "ID") |>
  dplyr::mutate(
    ID = glue::glue(
      "{stringr::str_split_i(basename(ID), pattern = '-', i = 1)} (MR)"
    )
  ) |>
  dplyr::select(ID, CHR, POS, GENE)

lead_vars <- dplyr::bind_rows(metal_vars, mega_vars) |>
  dplyr::mutate(
    START = POS - 500000,
    END = POS + 500000
  ) |>
  tidygenomics::genome_cluster(
    by = c("CHR", "START", "END"),
    max_distance = 0
  )

upset_list <- list(
  "Ventral (METAL)" = lead_vars |>
    dplyr::filter(ID == "Ventral (METAL)") |>
    dplyr::pull(cluster_id),
  "Inguinal (METAL)" = lead_vars |>
    dplyr::filter(ID == "Inguinal (METAL)") |>
    dplyr::pull(cluster_id),
  "Umbilical (METAL)" = lead_vars |>
    dplyr::filter(ID == "Umbilical (METAL)") |>
    dplyr::pull(cluster_id),
  "Femoral (METAL)" = lead_vars |>
    dplyr::filter(ID == "Femoral (METAL)") |>
    dplyr::pull(cluster_id),
  "Diaphragmatic (METAL)" = lead_vars |>
    dplyr::filter(ID == "Diaphragmatic (METAL)") |>
    dplyr::pull(cluster_id),
  "Ventral (MR-MEGA)" = lead_vars |>
    dplyr::filter(ID == "Ventral (MR)") |>
    dplyr::pull(cluster_id),
  "Inguinal (MR-MEGA)" = lead_vars |>
    dplyr::filter(ID == "Inguinal (MR)") |>
    dplyr::pull(cluster_id),
  "Umbilical (MR-MEGA)" = lead_vars |>
    dplyr::filter(ID == "Umbilical (MR)") |>
    dplyr::pull(cluster_id),
  "Diaphragmatic (MR-MEGA)" = lead_vars |>
    dplyr::filter(ID == "Diaphragmatic (MR)") |>
    dplyr::pull(cluster_id)
)

# Convert list to the input format expected by UpSetR
upset_mt <- UpSetR::fromList(upset_list)

# visualize
plot <- UpSetR::upset(
  upset_mt,
  sets = names(upset_list),
  order.by = "freq",
  nintersects = NA,
  text.scale = (10 / 7)
)

pdf("figures/metal_mrmega_upset.pdf", width = 7.5, height = 5.5)
print(plot)
dev.off()

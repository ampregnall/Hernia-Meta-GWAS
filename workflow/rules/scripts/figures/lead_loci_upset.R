# Load lead variants
lead_files <- c(
  "data/processed/lead-variants/metal/Ventral-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Inguinal-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Umbilical-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Femoral-ALL_lead_variants.txt",
  "data/processed/lead-variants/metal/Diaphragmatic-ALL_lead_variants.txt"
)

lead_vars <- vroom::vroom(lead_files, id = "ID") |>
  dplyr::mutate(
    ID = stringr::str_split_i(basename(ID), pattern = "-", i = 1),
    START = POS - 500000,
    END = POS + 500000
  ) |>
  tidygenomics::genome_cluster(
    by = c("CHR", "START", "END"),
    max_distance = 0
  )

upset_list <- list(
  "Ventral" = lead_vars |>
    dplyr::filter(ID == "Ventral") |>
    dplyr::pull(cluster_id),
  "Inguinal" = lead_vars |>
    dplyr::filter(ID == "Inguinal") |>
    dplyr::pull(cluster_id),
  "Umbilical" = lead_vars |>
    dplyr::filter(ID == "Umbilical") |>
    dplyr::pull(cluster_id),
  "Femoral" = lead_vars |>
    dplyr::filter(ID == "Femoral") |>
    dplyr::pull(cluster_id),
  "Diaphragmatic" = lead_vars |>
    dplyr::filter(ID == "Diaphragmatic") |>
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

pdf("figures/hernia_subtypes_upset.pdf", width = 7.5, height = 4)
print(plot)
dev.off()

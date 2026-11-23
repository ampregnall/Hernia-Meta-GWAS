proteome_files <- c(
  "data/processed/mendelian-randomization/archive/Diaphragmatic_mr.txt",
  "data/processed/mendelian-randomization/archive/Inguinal_mr.txt",
  "data/processed/mendelian-randomization/archive/Umbilical_mr.txt",
  "data/processed/mendelian-randomization/archive/Femoral_mr.txt",
  "data/processed/mendelian-randomization/archive/Ventral_mr.txt"
)

drug_mr_files <- c(
  "data/processed/mendelian-randomization/Diaphragmatic_DrugTargets_mr.txt",
  "data/processed/mendelian-randomization/Inguinal_DrugTargets_mr.txt",
  "data/processed/mendelian-randomization/Umbilical_DrugTargets_mr.txt",
  "data/processed/mendelian-randomization/Femoral_DrugTargets_mr.txt",
  "data/processed/mendelian-randomization/Ventral_DrugTargets_mr.txt"
)

df_proteome <- vroom(proteome_files, id = "Phenotype") |> 
  filter(stringr::str_detect(exposure, "Fenland")) |> 
  filter(method == "Inverse variance weighted") |> 
  mutate(
    Phenotype = stringr::str_split_i(basename(Phenotype), "_", 1),
    exposure = stringr::str_split_i(exposure, ":", 4),
    z = b / se
  ) |> 
  group_by(Phenotype) |> 
  mutate(fdr = p.adjust(pval, method = "fdr")) |> 
  ungroup() |> 
  filter(fdr < 0.05)

df_pqtl <- vroom("data/processed/mendelian-randomization/MR_PQTL_colocalization.txt") |> 
  filter(posterior_prob > 0.70) |> 
  mutate(PQTL = TRUE) |> 
  select(Phenotype, Gene, PQTL)

# Merge data
df_proteome <- left_join(df_proteome, df_pqtl, by=c("Phenotype", "exposure"="Gene")) |> 
  mutate(group = "1 Proteome MR", PQTL = if_else(is.na(PQTL), FALSE, TRUE)) |> 
  filter(PQTL)


df_drugmr <- vroom(drug_mr_files, id = "Phenotype") |> 
  filter(stringr::str_detect(exposure, "fibroblasts")) |> 
  mutate(
    Phenotype = stringr::str_split_i(basename(Phenotype), "_", 1),
    tissue = case_when(
      str_detect(exposure, "fibroblast") ~ "Fibroblasts",
      .default = "Error"),
    gene = str_split_i(exposure, " ", 1),
    z = b / se
  ) |>
  group_by(Phenotype) |> 
  mutate(fdr = p.adjust(pval, method = "fdr"), group = "2 Drug Repurposing", PQTL = FALSE) |> 
  ungroup() |> 
  filter(fdr < 0.05) |> 
  mutate(exposure = glue::glue("{gene}")) |> 
  select(!c(tissue, gene))

df_plot <- rbind(df_proteome, df_drugmr) |> 
  mutate(Phenotype = factor(Phenotype, levels = c("Ventral", "Inguinal", "Umbilical", "Femoral", "Diaphragmatic")))

p1 <- ggplot(df_plot, aes(x = Phenotype, y = exposure, color = factor(sign(z)))) + 
  geom_point(aes(size = abs(z))) +
  scale_shape_manual(values = c(19, 8)) +
  facet_wrap(~group, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("#011F5B", "#990000"), name = "Z-Score Sign") +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    panel.grid.minor = element_blank(),
    text = element_text(size = 10)
    )

ggsave(
  filename = "figures/mendelian-randomization/drug-repurposing-mr.pdf",
  plot = p1, 
  width = 360,
  height = 340,
  units = "mm",
  dpi = 600
)

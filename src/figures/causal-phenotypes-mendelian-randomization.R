library(ggforestplot)
library(dplyr)
library(vroom)
library(ggplot2)
library(ggstats)
library(cowplot)

mr_files <- c(
  "data/processed/mendelian-randomization/Diaphragmatic_mr.txt",
  "data/processed/mendelian-randomization/Inguinal_mr.txt",
  "data/processed/mendelian-randomization/Femoral_mr.txt",
  "data/processed/mendelian-randomization/Umbilical_mr.txt",
  "data/processed/mendelian-randomization/Ventral_mr.txt"
)

TRAITS_ANTHRO <- c(
  "Height_PMID36224396", 
  "Weight_PMID39024449",
  "BMI_PMID30239722",
  "WaistCircum_PMID40968291",
  "WHRatioBMIAdj_PMID30239722",
  "WBLeanMass_PMID40968291"
)

TRAITS_FAT <- c(
  "ASAT_PMID35773277",
  "ASATBMIAdj_PMID35773277",
  "GFAT_PMID35773277",
  "GFATBMIAdj_PMID35773277",
  "VAT_PMID35773277",
  "VATBMIAdj_PMID35773277"
)

trait_order <- c(
  TRAITS_ANTHRO,
  TRAITS_FAT
)

# Create a named vector for better labels
trait_labels <- c(
  "Height_PMID36224396" = "Height",
  "Weight_PMID39024449" = "Weight",
  "BMI_PMID30239722" = "BMI",
  "WaistCircum_PMID40968291" = "Waist Circumference",
  "WHRatioBMIAdj_PMID30239722" = "WHR (BMI-Adjusted)",
  "WBLeanMass_PMID40968291" = "Whole Body Lean Mass",
  "ASAT_PMID35773277" = "ASAT",
  "ASATBMIAdj_PMID35773277" = "ASAT (BMI-Adjusted)",
  "GFAT_PMID35773277" = "GFAT",
  "GFATBMIAdj_PMID35773277" = "GFAT (BMI-Adjusted)",
  "VAT_PMID35773277" = "VAT",
  "VATBMIAdj_PMID35773277" = "VAT (BMI-Adjusted)"
)

df <- vroom(mr_files, id = "Phenotype") |> 
  dplyr::select(!fdr) |> 
  group_by(Phenotype) |> 
  mutate(fdr = p.adjust(pval, method = "fdr")) |> 
  ungroup()  |> 
  mutate(Phenotype = stringr::str_split_i(basename(Phenotype), "_", 1)) |> 
  filter(exposure %in% c(TRAITS_FAT, TRAITS_ANTHRO)) |> 
  mutate(exposure = factor(exposure, levels = trait_order)) |> 
  group_by(Phenotype, exposure) |> 
  filter(method == "MR Egger" & nsnp == 1 | method == "Inverse variance weighted" & nsnp > 1) |> 
  ungroup() |> 
  mutate(group = case_when(
    exposure %in% TRAITS_ANTHRO ~ "Anthropometric",
    exposure %in% TRAITS_FAT ~ "Fat Distribution"
  )) |> 
  mutate(Phenotype = factor(Phenotype, levels = c("Ventral", "Inguinal", "Umbilical", "Femoral", "Diaphragmatic")))

df <- df |> 
  arrange(exposure) |> 
  mutate(OR = exp(b), ORHIGH = exp(b + 1.96 * se), ORLOW = exp(b - 1.96 * se)) |> 
  mutate(exposure_label = factor(trait_labels[as.character(exposure)], 
                                 levels = trait_labels[trait_order]))

# Define colors
hernia_colors <- c("Ventral" = "#4DBBD5", 
                   "Inguinal" = "#00A087", 
                   "Umbilical" = "#3C5488", 
                   "Femoral" = "#F39B7F", 
                   "Diaphragmatic" = "#8491B1")

# Create Anthropometric plot
p1 <- df |> 
  filter(group == "Anthropometric") |> 
  ggplot(aes(x = OR, y = exposure_label, color = Phenotype)) +
  geom_stripped_rows(odd = "#f5f5f5", even = "white", color = NA) +
  geom_point() +
  geom_errorbarh(aes(xmin = ORLOW, xmax = ORHIGH), height = 0.2, linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = hernia_colors) +
  facet_wrap(~Phenotype, nrow = 1, scales = "free_x") +
  labs(y = NULL, x = "Odds Ratio") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none", panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey30", color = NA),  # Dark grey background
        strip.text = element_text(color = "white")  # White text for contrast
  )

# Create Fat Distribution plot
p2 <- df |> 
  filter(group == "Fat Distribution") |> 
  ggplot(aes(x = OR, y = exposure_label, color = Phenotype)) +
  geom_stripped_rows(odd = "#f5f5f5", even = "white", color = NA) +
  geom_point() +
  geom_errorbarh(aes(xmin = ORLOW, xmax = ORHIGH), height = 0.2, linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = hernia_colors) +
  facet_wrap(~Phenotype, nrow = 1, scales = "free_x") +
  labs(y = NULL, x = "Odds Ratio") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none", panel.grid = element_blank(), ,
        strip.background = element_rect(fill = "grey30", color = NA),  # Dark grey background
        strip.text = element_text(color = "white")  # White text for contrast
  )

# Extract legend
legend <- get_legend(
  df |> 
    filter(group == "Anthropometric") |> 
    ggplot(aes(x = OR, y = exposure_label, color = Phenotype)) +
    geom_point() +
    scale_color_manual(values = hernia_colors) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "bottom")
)

# Combine plots
f1 <- plot_grid(
  plot_grid(p1, p2, ncol = 1, labels = c("A", "B"), label_size = 10),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

ggsave(
  filename = "figures/mendelian-randomization/mendelian-randomization.pdf",
  plot = f1,
  height = 170,
  width = 180,
  dpi = 600,
  units = "mm"
)  
  

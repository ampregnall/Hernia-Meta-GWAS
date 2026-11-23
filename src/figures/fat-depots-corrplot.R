base::load("data/processed/genomic-sem/ldsc/GenomicSEM_hernia_adipose_traits.RData")

trait.names <- c(
  "ASAT",
  "ASATAdjBMI",
  "ASATGFATRatio",
  "BMI",
  "Diaphragmatic",
  "Femoral", 
  "GFAT",
  "GFATAdjBMI",
  "Inguinal",
  "Umbilical",
  "VAT",
  "VATAdjBMI",
  "VATASATRatio",
  "VATGFATRatio",
  "Ventral"
)
S <- ldsc_output$S
R <- cov2cor(S)
colnames(R) <- rownames(R) <- trait.names
colnames(S) <- rownames(S) <- trait.names

p1 <- ggcorrplot(R) +
  theme(
    text = element_text(size = 8)
  )

ggsave(
  filename = "figures/corrplots/hernia-fat-depots-corrplot.pdf",
  plot = p1, 
  width = 170, 
  height = 170, 
  units = "mm",
  dpi = 600
)

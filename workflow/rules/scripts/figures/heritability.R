library(ggplot2)

# Read the data
df <- readr::read_delim("data/processed/heritability/ldsc_heritability.txt")

# Calculate total sample size and select top 2 rows per trait
df_plot <- df %>%
  dplyr::mutate(total_n = n_case + n_control) %>%
  dplyr::group_by(trait) %>%
  dplyr::slice_max(order_by = total_n, n = 2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    trait_prev = paste0(
      trait,
      " ",
      scales::percent(prev, accuracy = 10^(floor(log10(prev)) + 1))
    )
  )

# Create the plot
plt <- ggplot(df_plot, aes(x = trait_prev, y = h2_liability, color = trait)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = h2_liability - 1.96 * h2_liability_se,
      ymax = h2_liability + 1.96 * h2_liability_se
    ),
    width = 0.2
  ) +
  labs(
    x = "Trait (Population Prevalence)",
    y = "Liability-Scale Heritability (h²)",
    title = "Heritability Estimates by Trait and Prevalence"
  ) +
  ggsci::scale_color_npg() +
  coord_flip() +
  theme_minimal() +
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))

ggsave(
  filename = "figures/heritability.pdf",
  plot = plt,
  width = 7.5,
  height = 4,
  units = "in"
)

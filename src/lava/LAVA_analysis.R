
##### LOAD FILES ###############################################################
lava_uni <- readr::read_delim(
  "data/processed/lava/LAVA_hernia_univariate.txt",
  show_col_types = FALSE
)
lava_biv <- readr::read_delim(
  "data/processed/lava/LAVA_hernia_bivariate.txt",
  show_col_types = FALSE
)
lava_pcor <- readr::read_delim(
  "data/processed/lava/LAVA_hernia_pcor_BMIAdj.txt",
  show_col_types = FALSE
)

### HOW MANY SIGNIFICANT UNIVARIATE AND BIVARIATE LOCI FOR BMI
lava_uni_bmi <- lava_uni |> 
  filter(p < 0.05 / 243 & phen == "BMI")

lava_biv_sig_bmi <- lava_biv |> 
  filter(phen2 == "BMI") |> 
  filter(p < 0.05 / 495)

lava_biv_bmi <- lava_biv |> 
  filter(phen2 == "BMI") |> 
  group_by(phen1, phen2) |> 
  summarise(
    n_loci_tested = n(),
    n_nomsig = sum(p < 0.05, na.rm = TRUE),
    n_nomsig_pos = sum(p < 0.05 & rho > 0, na.rm = TRUE),
    n_nomsig_neg = sum(p < 0.05 & rho < 0, na.rm = TRUE),
    median_rho_sig = median(rho[p < 0.05 / 495], na.rm = TRUE),
    median_rho = median(rho, na.rm = TRUE),
    mean_rho_sig = mean(rho[p < 0.05 / 495], na.rm = TRUE),
    mean_rho = mean(rho, na.rm = TRUE)
  ) 



lava_biv_sig_nobmi <- lava_biv |> 
  filter(phen2 != "BMI") |> 
  filter(p < 0.05 / 961)

lava_biv_pairwise_count <- lava_biv_sig_nobmi |> 
  group_by(phen1, phen2) |> 
  summarise(count = n())


lava_biv_hh <- lava_biv %>%
  filter(phen1 != "BMI") %>%
  group_by(phen1, phen2) %>%
  summarise(
    n_loci_tested = n(),
    n_sig_p005 = sum(p < 0.05/ 961, na.rm = TRUE),
    n_sig_pos = sum(p < 0.05 / 961 & rho > 0, na.rm = TRUE),
    n_sig_neg = sum(p < 0.05 / 961 & rho < 0, na.rm = TRUE),
    median_rho_sig = median(rho[p < 0.05 / 961], na.rm = TRUE),
    mean_rho = mean(rho)
  ) %>%
  arrange(desc(n_sig_p005))

pcor_loci <- fread("data/processed/lava/LAVA_hernia_BMIAdj_locus_file.txt") |> 
  select(LOC) |> 
  mutate(locus_id = row_number())

lava_pcor_hh <- lava_pcor %>%
  filter(phen1 != "BMI") %>%
  group_by(phen1, phen2) %>%
  summarise(
    n_loci_tested = n(),
    n_sig_p005 = sum(p < 0.05/ 1494, na.rm = TRUE),
    n_sig_pos = sum(p < 0.05 / 1494 & pcor > 0, na.rm = TRUE),
    n_sig_neg = sum(p < 0.05 / 1494 & pcor < 0, na.rm = TRUE),
    median_rho_sig = median(pcor[p < 0.05 / 1494], na.rm = TRUE),
    mean_rho = mean(pcor, na.rm = TRUE)
  ) %>%
  arrange(desc(n_sig_p005))



lava_pcor <- lava_pcor |> 
  left_join(pcor_loci, by = "locus_id")
  

lava_rho2 <- lava_biv_nobmi %>%
  select(locus, phen1, phen2, rho, p_rho = p)

bmi_comparison <- lava_rho2 %>%
  left_join(
    lava_pcor %>% select(locus = LOC, phen1, phen2, pcor, p_pcor = p),
    by = c("locus", "phen1", "phen2")
  ) %>%
  mutate(delta = pcor - rho)

tbl_bmi_comparison <- bmi_comparison %>%
  filter(!is.na(pcor)) |> 
  group_by(phen1, phen2) %>%
  summarise(
    n_loci = n(),
    mean_rho = mean(rho, na.rm = TRUE),
    mean_pcor = mean(pcor, na.rm = TRUE),
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    pct_delta_pos = mean(delta > 0, na.rm = TRUE)
  )



p1 <- lava_pcor |> select(phen = phen1, r2 = r2.phen1_z)
p2 <- lava_pcor |> select(phen = phen2, r2 = r2.phen2_z)
df <- rbind(p1, p2)

df |> group_by(phen) |> summarise(avg = mean(r2), med = median(r2))

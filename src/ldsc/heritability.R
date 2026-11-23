library(ldscr)
library(data.table)
source("src/ldsc/munge_sumstats.R")

input_dir <- "data/processed/metal-formatted"
output_dir <- "data/processed/sumstats-munged/"
hm3 <- "data/meta/ldsc/w_hm3.snplist"

file_list <- grep(
  pattern = "Abdominal",
  x = list.files(
    input_dir,
    pattern = "_metal_formatted.txt",
    full.names = TRUE
  ),
  value = TRUE,
  invert = TRUE
)

for (file in file_list) {
  data <- vroom::vroom(file) |>
    dplyr::mutate(MAF = dplyr::if_else(EAF < 0.5, EAF, 1 - EAF)) |>
    dplyr::select(SNP = rsID, A1 = EA, A2 = NEA, MAF, BETA, P, N)

  trait_name <- gsub("_metal_formatted.txt", "", basename(file))
  munge_sumstats(
    data,
    hm3 = hm3,
    trait.names = c(trait_name),
    output_path = output_dir
  )
}

calculate_heritability <- function(
  file,
  n_case,
  n_control,
  prev,
  ancestry,
  trait
) {
  data <- data.table::fread(file)
  result <- ldscr::ldsc_h2(
    munged_sumstats = data,
    ancestry = ancestry,
    sample_prev = if (!is.na(n_case) && !is.na(n_control)) {
      n_case / (n_case + n_control)
    } else {
      NA
    },
    population_prev = prev
  )
  result$trait <- trait
  result$n_case <- n_case
  result$n_control <- n_control
  result$prev <- prev
  result$ancestry <- ancestry
  return(result)
}

# List of trait data
traits <- list(
  list(
    file = "data/processed/sumstats-munged/Ventral-ALL.sumstats.gz",
    n_case = 23354,
    n_control = 1596564,
    prev = 0.10,
    ancestry = "EUR",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Ventral-ALL.sumstats.gz",
    n_case = 23354,
    n_control = 1596564,
    prev = 0.30,
    ancestry = "EUR",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-ALL.sumstats.gz",
    n_case = 106421,
    n_control = 1729615,
    prev = 0.27,
    ancestry = "EUR",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-ALL.sumstats.gz",
    n_case = 106421,
    n_control = 1729615,
    prev = 0.43,
    ancestry = "EUR",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Umbilical-ALL.sumstats.gz",
    n_case = 31367,
    n_control = 1581741,
    prev = 0.1,
    ancestry = "EUR",
    trait = "Umbilical"
  ),
  list(
    file = "data/processed/sumstats-munged/Umbilical-ALL.sumstats.gz",
    n_case = 31367,
    n_control = 1581741,
    prev = 0.3,
    ancestry = "EUR",
    trait = "Umbilical"
  ),
  list(
    file = "data/processed/sumstats-munged/Femoral-ALL.sumstats.gz",
    n_case = 2446,
    n_control = 790529,
    prev = 0.016,
    ancestry = "EUR",
    trait = "Femoral"
  ),
  list(
    file = "data/processed/sumstats-munged/Femoral-ALL.sumstats.gz",
    n_case = 2446,
    n_control = 790529,
    prev = 0.073,
    ancestry = "EUR",
    trait = "Femoral"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-ALL.sumstats.gz",
    n_case = 81148,
    n_control = 1545371,
    prev = 0.00013,
    ancestry = "EUR",
    trait = "Diaphragmatic"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-ALL.sumstats.gz",
    n_case = 81148,
    n_control = 1545371,
    prev = 0.00023,
    ancestry = "EUR",
    trait = "Diaphragmatic"
  ),
  ### EUR
  list(
    file = "data/processed/sumstats-munged/Ventral-EUR.sumstats.gz",
    n_case = 21245,
    n_control = 1331234,
    prev = 0.10,
    ancestry = "EUR",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Ventral-EUR.sumstats.gz",
    n_case = 21245,
    n_control = 1331234,
    prev = 0.30,
    ancestry = "EUR",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-EUR.sumstats.gz",
    n_case = 93719,
    n_control = 1297379,
    prev = 0.27,
    ancestry = "EUR",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-EUR.sumstats.gz",
    n_case = 93719,
    n_control = 1297379,
    prev = 0.43,
    ancestry = "EUR",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Umbilical-EUR.sumstats.gz",
    n_case = 27318,
    n_control = 1320619,
    prev = 0.1,
    ancestry = "EUR",
    trait = "Umbilical"
  ),
  list(
    file = "data/processed/sumstats-munged/Umbilical-EUR.sumstats.gz",
    n_case = 27318,
    n_control = 1320619,
    prev = 0.3,
    ancestry = "EUR",
    trait = "Umbilical"
  ),
  list(
    file = "data/processed/sumstats-munged/Femoral-EUR.sumstats.gz",
    n_case = 2446,
    n_control = 790529,
    prev = 0.016,
    ancestry = "EUR",
    trait = "Femoral"
  ),
  list(
    file = "data/processed/sumstats-munged/Femoral-EUR.sumstats.gz",
    n_case = 2446,
    n_control = 790529,
    prev = 0.073,
    ancestry = "EUR",
    trait = "Femoral"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-EUR.sumstats.gz",
    n_case = 73070,
    n_control = 1288859,
    prev = 0.00013,
    ancestry = "EUR",
    trait = "Diaphragmatic"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-EUR.sumstats.gz",
    n_case = 73070,
    n_control = 1288859,
    prev = 0.00023,
    ancestry = "EUR",
    trait = "Diaphragmatic"
  ),
  ### AFR
  list(
    file = "data/processed/sumstats-munged/Ventral-AFR.sumstats.gz",
    n_case = 1258,
    n_control = 166521,
    prev = 0.10,
    ancestry = "AFR",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Ventral-AFR.sumstats.gz",
    n_case = 1258,
    n_control = 166521,
    prev = 0.30,
    ancestry = "AFR",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-AFR.sumstats.gz",
    n_case = 7284,
    n_control = 158627,
    prev = 0.27,
    ancestry = "AFR",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-AFR.sumstats.gz",
    n_case = 7284,
    n_control = 158627,
    prev = 0.43,
    ancestry = "AFR",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Umbilical-AFR.sumstats.gz",
    n_case = 2537,
    n_control = 163798,
    prev = 0.1,
    ancestry = "AFR",
    trait = "Umbilical"
  ),
  list(
    file = "data/processed/sumstats-munged/Umbilical-AFR.sumstats.gz",
    n_case = 2537,
    n_control = 163798,
    prev = 0.3,
    ancestry = "AFR",
    trait = "Umbilical"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-AFR.sumstats.gz",
    n_case = 4635,
    n_control = 157207,
    prev = 0.00013,
    ancestry = "AFR",
    trait = "Diaphragmatic"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-AFR.sumstats.gz",
    n_case = 4635,
    n_control = 157207,
    prev = 0.00023,
    ancestry = "AFR",
    trait = "Diaphragmatic"
  ),
  ### AMR
  list(
    file = "data/processed/sumstats-munged/Ventral-AMR.sumstats.gz",
    n_case = 793,
    n_control = 90972,
    prev = 0.10,
    ancestry = "AMR",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Ventral-AMR.sumstats.gz",
    n_case = 793,
    n_control = 90972,
    prev = 0.30,
    ancestry = "AMR",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-AMR.sumstats.gz",
    n_case = 3387,
    n_control = 87335,
    prev = 0.27,
    ancestry = "AMR",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-AMR.sumstats.gz",
    n_case = 3387,
    n_control = 87335,
    prev = 0.43,
    ancestry = "AMR",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Umbilical-AMR.sumstats.gz",
    n_case = 1422,
    n_control = 89487,
    prev = 0.1,
    ancestry = "AMR",
    trait = "Umbilical"
  ),
  list(
    file = "data/processed/sumstats-munged/Umbilical-AMR.sumstats.gz",
    n_case = 1422,
    n_control = 89487,
    prev = 0.3,
    ancestry = "AMR",
    trait = "Umbilical"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-AMR.sumstats.gz",
    n_case = 2556,
    n_control = 87512,
    prev = 0.00013,
    ancestry = "AMR",
    trait = "Diaphragmatic"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-AMR.sumstats.gz",
    n_case = 2556,
    n_control = 87512,
    prev = 0.00023,
    ancestry = "AMR",
    trait = "Diaphragmatic"
  ),
  ### EAS/CSA
  list(
    file = "data/processed/sumstats-munged/Ventral-AMR.sumstats.gz",
    n_case = 58,
    n_control = 7837,
    prev = 0.10,
    ancestry = "CSA",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Ventral-AMR.sumstats.gz",
    n_case = 58,
    n_control = 7837,
    prev = 0.30,
    ancestry = "CSA",
    trait = "Ventral"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-AMR.sumstats.gz",
    n_case = 1664,
    n_control = 177062,
    prev = 0.27,
    ancestry = "EAS",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Inguinal-AMR.sumstats.gz",
    n_case = 1664,
    n_control = 177062,
    prev = 0.43,
    ancestry = "EAS",
    trait = "Inguinal"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-AMR.sumstats.gz",
    n_case = 88,
    n_control = 2581,
    prev = 0.00013,
    ancestry = "EAS",
    trait = "Diaphragmatic"
  ),
  list(
    file = "data/processed/sumstats-munged/Diaphragmatic-AMR.sumstats.gz",
    n_case = 88,
    n_control = 2581,
    prev = 0.00023,
    ancestry = "EAS",
    trait = "Diaphragmatic"
  )
)

ldsc_heritability <- rbindlist(
  lapply(traits, function(x) {
    calculate_heritability(
      x$file,
      x$n_case,
      x$n_control,
      x$prev,
      x$ancestry,
      x$trait
    )
  }),
  fill = TRUE
)

data.table::fwrite(
  ldsc_heritability,
  file = "data/processed/heritability/ldsc_heritability.txt"
)

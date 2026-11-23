library(LAVA)
library(data.table)
library(vroom)
library(dplyr)
library(tidygenomics)
library(glue)

convert_to_lava_format <- function(input_df, 
                                   output_file = NULL, 
                                   snp_col = "rsID", 
                                   chr_col = "CHR",
                                   bp_col = "POS",
                                   a1_col = "EA",
                                   a2_col = "NEA",
                                   n_col = "N",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   p_col = "P") {
  
  # Read your summary statistics
  sumstats <- input_df
  
  # Calculate Z-score from BETA and SE if not already present
  if ("Z" %in% names(sumstats)) {
    z_score <- sumstats$Z
  } else {
    z_score <- sumstats[[beta_col]] / sumstats[[se_col]]
  }
  
  # Create LAVA-formatted data frame with standardized column names
  lava_sumstats <- data.table(
    SNP = sumstats[[snp_col]],
    CHR = sumstats[[chr_col]],
    BP = sumstats[[bp_col]],
    A1 = sumstats[[a1_col]],
    A2 = sumstats[[a2_col]],
    N = sumstats[[n_col]],
    Z = z_score,
    P = sumstats[[p_col]]
  )
  
  # Create temporary file if output_file not specified
  if (is.null(output_file)) {
    output_file <- tempfile(fileext = ".txt")
    message("Created temporary file: ", output_file)
  }
  
  # Write to file
  fwrite(lava_sumstats, output_file, sep = "\t", quote = FALSE)
  
  return(output_file)
}

# ============================================================================
# STEP 2: Convert sumstats to temporary files
# ============================================================================

ventral <- fread("data/processed/metal-formatted/Ventral-ALL_metal_formatted.txt")
inguinal <- fread("data/processed/metal-formatted/Inguinal-ALL_metal_formatted.txt")
umbilical <- fread("data/processed/metal-formatted/Umbilical-ALL_metal_formatted.txt")
femoral <- fread("data/processed/metal-formatted/Femoral-ALL_metal_formatted.txt")
diaphragmatic <- fread("data/processed/metal-formatted/Diaphragmatic-ALL_metal_formatted.txt")
bmi <- fread("data/meta/lava/BMI_sumstats.txt.gz")
bmi$N <- 806834

ventral_lava <- convert_to_lava_format(input_df = ventral)
inguinal_lava <- convert_to_lava_format(input_df = inguinal)
umbilical_lava <- convert_to_lava_format(input_df = umbilical)
femoral_lava <- convert_to_lava_format(input_df = femoral)
diaphragmatic_lava <- convert_to_lava_format(input_df = diaphragmatic)
bmi_lava <- convert_to_lava_format(
  input_df = bmi, 
  snp_col = "SNP", 
  chr_col = "CHR",
  bp_col = "BP", 
  a1_col = "A1", 
  a2_col = "A2", 
  n_col = "N", 
  beta_col = "BETA", 
  se_col = "SE", 
  p_col = "P"
)

# ============================================================================
# STEP 3: Create input info file
# ============================================================================

input_info <- data.frame(
  phenotype = c("Ventral", "Inguinal", "Umbilical", "Femoral", "Diaphragmatic", "BMI"),
  cases = c(23354, 106421, 31367, 2446, 81148, NA),      
  controls = c(1596564, 1729615, 1581741, 790529, 1545371, NA),   
  prevalence = c(NA, NA, NA, NA, NA, NA), # Can include for local h2 estimates on liability scale
  filename = c(ventral_lava, inguinal_lava, umbilical_lava, femoral_lava, diaphragmatic_lava, bmi_lava),
  stringsAsFactors = FALSE
)

# Save to temporary file
input_info_file <- tempfile(fileext = ".txt")
write.table(input_info, 
            input_info_file, 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

message("Input info file created: ", input_info_file)

# ============================================================================
# STEP 4: Create locus file
# ============================================================================

locus_file <- vroom("data/processed/lava/LAVA_hernia_univariate.txt") |> 
  filter(phen == "BMI" & p < (0.05 / length(unique(locus)))) |> 
  select(LOC = locus, CHR = chr, START = start, STOP = stop) |> 
  as.data.table()

setkey(locus_file, CHR, START, STOP)

get_locus_variants <- function(df) {
  df <- df |> mutate(POS_START = POS, POS_END = POS)
  df <- foverlaps(df, locus_file, by.x=c("CHR", "POS_START", "POS_END"), by.y = c("CHR", "START", "STOP"), nomatch = NULL)
  return(df)
}

ventral <- get_locus_variants(ventral)
inguinal <- get_locus_variants(inguinal)
umbilical <- get_locus_variants(umbilical)
femoral <- get_locus_variants(femoral)
diaphragmatic <- get_locus_variants(diaphragmatic)

variants <- bind_rows(ventral, inguinal, umbilical, femoral, diaphragmatic) |> 
  select(LOC, rsID) |> 
  distinct() |> 
  filter(!is.na(rsID) & rsID != "") |> 
  group_by(LOC) |> 
  summarise(SNPS = stringr::str_flatten(rsID, collapse = ";"))

locus_file <- left_join(locus_file, variants, by="LOC")
fwrite(locus_file, "data/processed/lava/LAVA_hernia_BMIAdj_locus_file.txt")

# ============================================================================
# STEP 5: Process input and run LAVA
# ============================================================================
ref_prefix <- "data/meta/lava/lava-ukb-v1.1"

# Read locus info
loci <- read.loci("data/processed/lava/LAVA_hernia_BMIAdj_locus_file.txt")
n.loc = nrow(loci)

# Process input - this reads and aligns your summary statistics
input <- process.input(
  input.info.file = input_info_file,
  sample.overlap.file = "data/processed/lava/LAVA_sample_overlap.txt",  # Optional: specify if samples overlap
  ref.prefix = ref_prefix,
  phenos = c("Ventral", "Inguinal", "Umbilical", "Femoral", "Diaphragmatic", "BMI")
)

# Define your hernia phenotypes
hernia_phenos <- c("Ventral", "Inguinal", "Umbilical", "Femoral", "Diaphragmatic")
condition_pheno <- "BMI"

# Generate all pairwise combinations
pairwise_combos <- combn(hernia_phenos, 2, simplify = FALSE)

# Initialize results list
pcor_results <- list()

# Loop through each locus (assuming you have a list of processed loci)
for(i in 1:n.loc) {
  locus = process.locus(loci[i,], input)                                          
  print(paste("Starting LAVA analysis for locus", i, "of", n.loc))
  
  # Loop through all pairwise combinations
  for(j in 1:length(pairwise_combos)) {
    pair <- pairwise_combos[[j]]
    
    tryCatch({
      pcor_result <- run.pcor(
        locus, 
        target = pair,
        phenos = condition_pheno
      )
      
      # Store results
      result_name <- paste0("locus_", i, "_", pair[1], "_vs_", pair[2])
      pcor_results[[result_name]] <- pcor_result
      
    }, error = function(e) {
      message(paste("Error in locus", i, "for pair", paste(pair, collapse="-"), ":", e$message))
    })
  }
}

# Combine all results into a single data frame
pcor_all <- do.call(rbind, pcor_results)

# Add locus information
pcor_all$locus_id <- sapply(strsplit(rownames(pcor_all), "_"), function(x) x[2])

# Save results
fwrite(pcor_all, "data/processed/lava/LAVA_hernia_pcor_BMIAdj.txt", sep = "\t", row.names = FALSE)
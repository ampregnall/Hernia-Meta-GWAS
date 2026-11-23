# Extract all variants +/- 500kb of genome-wide significant variant
locus_extract <- function(sumstats_df, # df containing summary statistics
                          sumstats_chr_col = chromosome,
                          sumstats_pos_col = position,
                          locus_df, # df containing gws loci
                          locus_chr_col = chromosome,
                          locus_pos_col = position,
                          locus_gene_col = gene,
                          locus_size = 1e6) {
  
  locus_df <- locus_df %>%
    dplyr::select(chromosome = {{ locus_chr_col }}, lead_pos = {{ locus_pos_col }}, gene = {{ locus_gene_col }}) %>%
    dplyr::distinct(chromosome, lead_pos, .keep_all = TRUE)
  
  sumstats_df <- sumstats_df %>%
    dplyr::rename(chromosome = {{ sumstats_chr_col }}, position = {{ sumstats_pos_col }}) %>%
    dplyr::mutate(chromosome = as.numeric(chromosome), position = as.numeric(position)) %>%
    dplyr::inner_join(locus_df, by = "chromosome", relationship = "many-to-many") %>%
    dplyr::filter(between(position, lead_pos - locus_size / 2, lead_pos + locus_size / 2)) %>%
    dplyr::collect() %>%
    dplyr::mutate(locus_marker = glue::glue("{chromosome}:{lead_pos}")) %>%
    dplyr::select(locus_marker, chromosome, lead_pos, gene, everything())
  
  return(sumstats_df)
}

locus_extract_dtable <- function(sumstats_df, 
                                 locus_df, 
                                 locus_size = 1e6,
                                 sumstats_chr_col, 
                                 sumstats_pos_col,
                                 locus_chr_col, 
                                 locus_pos_col, 
                                 locus_gene_col) {
  
  # Prepare locus_df
  locus_dt <- locus_df %>%
    dplyr::select(chromosome = {{ locus_chr_col }}, lead_pos = {{ locus_pos_col }}, gene = {{ locus_gene_col }}) %>%
    dplyr::distinct(chromosome, lead_pos, .keep_all = TRUE) %>%
    dplyr::mutate(
      start = lead_pos - locus_size / 2,
      end = lead_pos + locus_size / 2
    ) %>%
    data.table::as.data.table()
  
  # Prepare sumstats_df
  sumstats_dt <- sumstats_df %>%
    dplyr::mutate(start = {{ sumstats_pos_col }}, end = {{ sumstats_pos_col }}) %>%
    dplyr::rename(chromosome = {{ sumstats_chr_col }}) %>%
    dplyr::mutate(chromosome = as.numeric(chromosome), start = as.numeric(start), end = as.numeric(end)) %>%
    data.table::as.data.table()
  
  # Set keys for efficient joining
  data.table::setkey(locus_dt, chromosome, start, end)
  data.table::setkey(sumstats_dt, chromosome, start, end)
  
  # Efficient overlap join
  result <- data.table::foverlaps(
    sumstats_dt, locus_dt,
    by.x = c("chromosome", "start", "end"),
    by.y = c("chromosome", "start", "end"),
    type = "within"
  ) %>%
    dplyr::filter(!is.na(lead_pos)) %>%  # Remove non-matching rows
    dplyr::mutate(locus_marker = glue::glue("{chromosome}:{lead_pos}")) %>%
    dplyr::select(locus_marker, chromosome, lead_pos, gene, everything(), -i.start, -i.end)
  
  return(result)
}

locus_extract_coloc <- function(sumstats_df, # df containing summary statistics
                          sumstats_chr_col = chromosome,
                          sumstats_pos_col = position,
                          locus_df, # df containing gws loci
                          locus_chr_col = chromosome,
                          locus_start_col = START_region,
                          locus_end_col = END_region,
                          locus_id_col = ID) {
  
  locus_df <- locus_df %>%
    dplyr::select(chromosome = {{ locus_chr_col }}, start = {{ locus_start_col }}, 
                  end = {{ locus_end_col }}, locus_marker = {{ locus_id_col }})
  
  sumstats_df <- sumstats_df %>%
    dplyr::rename(chromosome = {{ sumstats_chr_col }}, position = {{ sumstats_pos_col }}) %>%
    dplyr::mutate(chromosome = as.numeric(chromosome), position = as.numeric(position)) %>%
    dplyr::inner_join(locus_df, by = "chromosome", relationship = "many-to-many") %>%
    dplyr::filter(between(position, start, end)) %>%
    dplyr::collect() %>%
    dplyr::select(locus_marker, chromosome, everything())
  
  return(sumstats_df)
}

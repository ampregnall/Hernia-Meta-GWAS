metal_config <- function(
  config_name,
  output_dir,
  study_files,
  SCHEME = "STDERR",
  AVERAGEFREQ = "ON",
  MINMAXFREQ = "OFF",
  TRACKPOSITIONS = "ON",
  MARKERLABEL = "MARKER",
  CHROMOSOMELABEL = "CHROM",
  POSITIONLABEL = "POS",
  EFFECT_ALLELE = "EFFECT_ALLELE",
  OTHER_ALLELE = "OTHER_ALLELE",
  EFFECTLABEL = "BETA",
  STDERR = "SE",
  FREQLABEL = "EAF",
  NCASE = NULL,
  NCONTROL = NULL,
  SAMPLESIZE = NULL
) {
  fs::dir_create(output_dir, recurse = TRUE)

  config_outfile <- fs::path(
    output_dir,
    paste0(config_name, "_metal-config.txt")
  )
  meta_outfile <- fs::path(normalizePath(output_dir), config_name)

  study_files <- paste0(
    glue::glue("PROCESS {normalizePath(study_files)}"),
    collapse = "\n"
  )

  config_text <- glue::glue(
    "SCHEME {SCHEME}
    AVERAGEFREQ {AVERAGEFREQ}
    MINMAXFREQ {MINMAXFREQ}
    TRACKPOSITIONS {TRACKPOSITIONS}
    MARKERLABEL {MARKERLABEL}
    CHROMOSOMELABEL {CHROMOSOMELABEL}
    POSITIONLABEL {POSITIONLABEL}
    ALLELELABELS {EFFECT_ALLELE} {OTHER_ALLELE}
    EFFECTLABEL {EFFECTLABEL}
    STDERR {STDERR}
    FREQLABEL {FREQLABEL}

    "
  )

  if (!is.null(NCASE)) {
    config_text <- glue::glue(
      "
    {config_text}\n
    CUSTOMVARIABLE NCASE
    LABEL NCASE as {NCASE}
    "
    )
  }

  if (!is.null(NCONTROL)) {
    config_text <- glue::glue(
      "
    {config_text}\n
    CUSTOMVARIABLE NCONTROL
    LABEL NCONTROL as {NCONTROL}"
    )
  }

  if (!is.null(SAMPLESIZE)) {
    config_text <- glue::glue(
      "
    {config_text}\n
    CUSTOMVARIABLE SAMPLESIZE
    LABEL SAMPLESIZE as {SAMPLESIZE}"
    )
  }

  config_text <- glue::glue(
    "
    {config_text}

    {study_files}

    OUTFILE {meta_outfile}_metal- .txt
    ANALYZE HETEROGENEITY

    QUIT"
  )

  readr::write_lines(x = config_text, file = config_outfile)

  return(config_outfile)
}

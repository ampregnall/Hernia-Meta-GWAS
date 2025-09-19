annotateVariantConsequences <- function(variantIds) {
  
  results <- purrr::map_dfr(variantIds, variantConsequence) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!is.na(MS_GENE))
  
  return(results)
}

variantConsequence <- function(variantId) {
  if (missing(variantId) || is.null(variantId)) {
    
    logger::log_error("The 'variantId' parameter is required. Please provide a valid variant ID.")
    stop()
  }
  
  query <- '
  query variantCredibleSets($variantId: String!) {
  variant(variantId: $variantId) {
    id
    referenceAllele
    alternateAllele
    mostSevereConsequence {id label}
    transcriptConsequences {
      variantConsequences {
        id
        label
      }
    target {
        id
        approvedSymbol
        biotype
      }
    }
  }
}'

  logger::log_info("Sending GraphQL request...")
  
  resp <- httr::POST(
    url = "https://api.platform.opentargets.org/api/v4/graphql",
    httr::add_headers(`Content-Type` = "application/json"),
    body = jsonlite::toJSON(
      list(
        query = query,
        variables = list(
          variantId = variantId
        )
      ),
      auto_unbox = TRUE,
    )
  )
  
  if (httr::status_code(resp) != 200) {
    logger::log_error("GraphQL query failed. HTTP status: ", httr::status_code(resp))
    stop()
  }
  
  output <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
  
  variant <- output$data$variant
  transcript <- output$data$variant$transcriptConsequences
  
  if (is.null(variant) || nrow(tibble::as_tibble(variant$mostSevereConsequence)) == 0) {
    message("No variant consequence data found for given SNP")
    return(dplyr::tibble(SNPID = variantId, MS_CONSEQUENCE = NA, MS_GENE = NA, MS_BIOTYPE = NA))
  }
  
  if (is.null(transcript) || nrow(tibble::as_tibble(transcript)) == 0) {
    message("No transcript consequence data found for given SNP")
    return(dplyr::tibble(SNPID = variantId, MS_CONSEQUENCE = NA, MS_GENE = NA, MS_BIOTYPE = NA))
  }
  
  most_severe_id <- variant$mostSevereConsequence$id
  most_severe_label <- variant$mostSevereConsequence$label
  
  transcript <- transcript %>% 
    tidyr::unnest(variantConsequences) %>%
    dplyr::filter(id %in% most_severe_id & label %in% most_severe_label) %>%
    dplyr::mutate(SNPID = variantId) %>%
    dplyr::select(SNPID, MS_CONSEQUENCE = label, MS_GENE = target.approvedSymbol, MS_BIOTYPE = target.biotype)
  
  return(transcript)
  
}

leadVariantsPhewas <- function(variantIds) {
  
  results <- purrr::map_dfr(variantIds, variantPhewas) %>%
    dplyr::bind_rows()
  
  return(results)
  
}

variantPhewas <- function(variantId) {
  
  if (missing(variantId) || is.null(variantId)) {
    
    logger::log_fatal("The 'variantId' parameter is required. Please provide a valid variant ID.")
    stop()
    }

  
  query <- '
  query variantCredibleSets($variantId: String!) {
  variant(variantId: $variantId) {
    id
    referenceAllele
    alternateAllele
    mostSevereConsequence {id label}
    gwasCredibleSets : credibleSets(studyTypes : [gwas]) {
      count
      rows {
        studyLocusId
        pValueMantissa
        pValueExponent
        beta
        standardError
        finemappingMethod
        confidence
        variant {
          id
          chromosome
          position
          referenceAllele
          alternateAllele
        }
        study {
          traitFromSource
          id
          diseases {
            name
            id
          }
        }
      }
    }
  }
  }'
  
  logger::log_info("Sending GraphQL request...")
  
  resp <- httr::POST(
    url = "https://api.platform.opentargets.org/api/v4/graphql",
    httr::add_headers(`Content-Type` = "application/json"),
    body = jsonlite::toJSON(
      list(
        query = query,
        variables = list(
          variantId = variantId
        )
      ),
      auto_unbox = TRUE,
    )
  )
  
  if (httr::status_code(resp) != 200) {
    stop("GraphQL query failed. HTTP status: ", httr::status_code(resp))
  }
  
  output <- jsonlite::fromJSON(httr::content(resp, as = "text"), flatten = TRUE)
  
  rows <- output$data$variant$gwasCredibleSets$rows
  
  if (is.null(rows) || nrow(tibble::as_tibble(rows)) == 0) {
    message("No GWAS credible sets data found for the given parameters.")
    return(NULL)
  }
  
  # Convert rows to tibble and check for credible set fields
  credible_data <- tibble::as_tibble(rows) %>%
    tidyr::unnest_wider(study.diseases) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(name = stringr::str_flatten(name, collapse = ";"),
                  id = stringr::str_flatten(id, collapse = ";")) %>%
    dplyr::select(SNPID = variant.id,
                  CHR = variant.chromosome,
                  POS = variant.position,
                  REF = variant.referenceAllele,
                  ALT = variant.alternateAllele,
                  BETA = beta,
                  SE = standardError,
                  P_MANTISSA = pValueMantissa,
                  P_EXPONENT = pValueExponent,
                  FINEMAPPING_METHOD = finemappingMethod,
                  STUDY_TRAIT = study.traitFromSource,
                  STUDY_ID = study.id,
                  DISEASE_NAME = name,
                  DISEASE_ID = id)
  
  return(credible_data)
}
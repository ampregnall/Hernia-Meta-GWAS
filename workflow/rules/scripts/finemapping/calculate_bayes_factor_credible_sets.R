calc_credset <- function(
  df,
  locus_marker_col = locus_marker,
  effect_col = effect,
  se_col = std_err,
  samplesize_col = samplesize,
  cred_interval = 0.99
) {
  df %>%
    dplyr::group_by({{ locus_marker_col }}) %>%
    dplyr::mutate(
      bf = exp(
        0.5 * ({{ effect_col }}^2 / {{ se_col }}^2 - log({{ samplesize_col }}))
      )
    ) %>%
    dplyr::mutate(posterior_prob = bf / sum(bf)) %>%
    dplyr::arrange(dplyr::desc(posterior_prob)) %>%
    dplyr::mutate(cum_sum = cumsum(posterior_prob)) %>%
    dplyr::group_by({{ locus_marker_col }}) %>%
    dplyr::filter(cum_sum <= cred_interval | posterior_prob > cred_interval) %>%
    dplyr::select(-cum_sum) %>%
    dplyr::ungroup()
}

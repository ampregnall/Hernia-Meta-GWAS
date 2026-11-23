# ---- Manhattan with a "broken" y-axis using facets (ggrepel-friendly) ----
gg_manhattan_df <- function(
  sumstats_df,
  annotation_df = NULL,
  phenotype = NULL,
  pval_threshold = 0.001,
  log10p_threshold = 3,
  label_col = GENE,
  build = "hg38",
  color1 = "#045ea7",
  color2 = "#82afd3",
  speed = "slow",
  use_broken_axis = TRUE,
  gap_abs_threshold = 25,
  pad_frac = 0.10,
  nonsig_sample_fraction = 0.25,
  label_pad = 12 # headroom (in -log10 p) for upper-facet labels
) {
  by <- dplyr::join_by(CHR, dplyr::between(x$POS, y$START, y$END))

  annotation_df <- annotation_df |>
    dplyr::filter(stringr::str_detect(ID, phenotype)) |>
    dplyr::mutate(log10p = -log10(P))

  ## 1) Find the largest gap between lead peaks
  top_vals <- sort(annotation_df$log10p, decreasing = TRUE)
  need_break <- FALSE
  gap_lower <- gap_upper <- NA_real_

  if (length(top_vals) >= 2) {
    p_val_diff <- top_vals[-length(top_vals)] - top_vals[-1]
    big_i <- which.max(p_val_diff)
    upper <- top_vals[big_i]
    lower <- top_vals[big_i + 1]
    gap_abs <- p_val_diff[big_i]

    # pad the gap so peaks flanking it remain visible
    gap_lower <- lower + pad_frac * (upper - lower)
    gap_upper <- upper - pad_frac * (upper - lower)

    need_break <- isTRUE(use_broken_axis) && (gap_abs >= gap_abs_threshold)
  }

  df <- sumstats_df |>
    dplyr::select(CHR, POS, P) |>
    dplyr::filter(P < pval_threshold) |>
    dplyr::mutate(dplyr::across(.cols = dplyr::everything(), as.numeric)) |>
    tidyr::drop_na() |>
    dplyr::mutate(log10p = -log10(P)) |>
    dplyr::mutate(is_significant = log10p > -log10(5e-8))

  sig_points <- df |> dplyr::filter(is_significant)
  nonsig_points <- df |>
    dplyr::filter(!is_significant) |>
    dplyr::slice_sample(prop = nonsig_sample_fraction)

  df <- dplyr::bind_rows(sig_points, nonsig_points) |>
    dplyr::select(-is_significant) |>
    dplyr::left_join(
      annotation_df |> dplyr::select(CHR, START, END, shared_locus),
      by
    ) |>
    dplyr::left_join(
      annotation_df |>
        dplyr::select(CHR, POS, P, label = {{ label_col }}) |>
        dplyr::slice_min(order_by = P, n = 10) |>
        dplyr::select(!P),
      by = c("CHR", "POS")
    ) |>
    dplyr::mutate(
      highlight = dplyr::case_when(
        log10p > -log10(5e-8) & shared_locus == TRUE ~ "#990000",
        log10p > -log10(5e-8) ~ "#d97373",
        TRUE ~ NA_character_
      ),
      pvalue = 10^(-log10p),
      chr = CHR,
      pos = POS
    ) |>
    dplyr::select(-CHR, -POS)

  max_log10_p <- max(df$log10p, na.rm = TRUE)
  min_log10p <- min(df$log10p, na.rm = TRUE)

  ## 3) If breaking, assign a facet "panel" and drop points in the removed gap
  if (need_break) {
    df <- df |>
      dplyr::mutate(
        panel = dplyr::case_when(
          log10p <= gap_lower ~ "low",
          log10p >= gap_upper ~ "high",
          TRUE ~ NA_character_ # rows inside the removed gap -> dropped
        )
      ) |>
      dplyr::filter(!is.na(panel))
  }

  ## 4) Base Manhattan from ggfastman
  plot <- df |>
    ggfastman::fast_manhattan(
      build = build,
      color1 = color1,
      color2 = color2,
      speed = speed
    ) +
    ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dotted") +
    ggplot2::labs(title = phenotype) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
      axis.title.x = element_blank()
    )

  if (need_break) {
    plot <- plot +
      ggplot2::facet_grid(
        rows = vars(panel),
        scales = "free_y",
        space = "free_y"
      ) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank()
      ) +
      ggh4x::facetted_pos_scales(
        y = list(
          panel == "high" ~
            scale_y_continuous(limits = c(gap_upper, max_log10_p)),
          panel == "low" ~ scale_y_continuous(limits = c(min_log10p, gap_lower))
        )
      ) +
      labs(y = "-log<sub>10</sub>(p-value)") +
      theme(axis.title.y = ggtext::element_markdown())
  } else {
    plot <- plot +
      ggplot2::scale_y_continuous(
        name = "-log<sub>10</sub>(p-value)",
        expand = expansion(mult = c(0.01, 0.07))
      ) +
      theme(
        panel.grid = element_blank(),
        axis.title.y = ggtext::element_markdown()
      )
  }

  ## 6) Labels (work in both modes). No ylim -> let facets handle limits.
  lab_df <- dplyr::filter(plot@data, !is.na(label))
  plot <- plot +
    ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = label),
      size = 8 / ggplot2::.pt,
      color = "black",
      direction = "x",
      nudge_y = 10,
      force_pull = 0,
      angle = 0,
      hjust = 0,
      segment.size = 0.2,
      segment.curvature = -1e-20,
      segment.angle = 175,
      max.overlaps = 50,
      ylim = c(max_log10_p, NA),
      seed = 1
    )

  return(plot)
}

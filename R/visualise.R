
#' @import ggplot2
NULL

#' Plot Time Series with Phenotype-Coloured Anomaly Bands
#'
#' @description
#' The primary summary figure. Draws the full series as a grey line and overlays
#' each detected anomaly segment coloured by phenotype, with a shaded band and
#' confidence-weighted alpha.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param title Character. Plot title. Default: auto-generated.
#' @param show_labels Logical. Show phenotype ID labels above bands. Default: `TRUE`.
#' @param alpha_base Numeric. Base alpha for anomaly lines. Default: 0.9.
#'
#' @return A `ggplot2` object (paper-quality Figure 1).
#' @export
plot_phenotype_timeline <- function(result, title = NULL,
                                     show_labels = TRUE, alpha_base = 0.9) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")
  x        <- result$x
  n        <- length(x)
  seg_len  <- result$seg_len
  inst_tab <- result$instance_table
  n_ph     <- result$mapper$n_phenotypes
  ph_df    <- result$phenotypes

  pal <- grDevices::hcl.colors(max(n_ph, 2L), palette = "Dynamic")
  names(pal) <- as.character(seq_len(n_ph))

  series_df <- data.frame(time = seq_len(n), value = x)

  seg_rows <- lapply(seq_len(nrow(inst_tab)), function(i) {
    st <- inst_tab$anomaly_idx[i]
    en <- min(n, st + seg_len - 1L)
    ph <- as.character(inst_tab$phenotype_id[i])
    conf <- inst_tab$phenotype_confidence[i]
    data.frame(time = st:en, value = x[st:en],
               phenotype = ph, conf = conf, stringsAsFactors = FALSE)
  })
  seg_df <- do.call(rbind, seg_rows)

  band_rows <- lapply(seq_len(nrow(inst_tab)), function(i) {
    st <- inst_tab$anomaly_idx[i]
    en <- min(n, st + seg_len - 1L)
    ph <- as.character(inst_tab$phenotype_id[i])
    data.frame(xmin = st - 0.5, xmax = en + 0.5,
               ymin = -Inf, ymax = Inf, phenotype = ph,
               stringsAsFactors = FALSE)
  })
  band_df <- do.call(rbind, band_rows)

  ph_labels <- paste0("P", seq_len(n_ph), ": ",
                       ph_df$name[order(ph_df$phenotype_id)])

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = band_df,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax,
                   fill = .data$phenotype),
      alpha = 0.18, inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = series_df,
      ggplot2::aes(x = .data$time, y = .data$value),
      colour = "#444444", linewidth = 0.45, alpha = 0.75
    ) +
    ggplot2::geom_line(
      data = seg_df,
      ggplot2::aes(x = .data$time, y = .data$value,
                   colour = .data$phenotype,
                   group  = interaction(.data$phenotype, .data$time),
                   alpha  = .data$conf),
      linewidth = 2.0
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.5, 1.0), guide = "none") +
    ggplot2::scale_fill_manual(values = pal, labels = ph_labels, name = "Phenotype") +
    ggplot2::scale_colour_manual(values = pal, guide  = "none") +
    ggplot2::labs(
      title    = if (is.null(title))
        "AnomalyPheno: Phenotype-Coloured Anomaly Timeline" else title,
      subtitle = sprintf(
        "%d instance(s) . %d phenotype(s)  |  line opacity = phenotype confidence",
        nrow(inst_tab), n_ph),
      x = "Time Index", y = "Value"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle   = ggplot2::element_text(colour = "#666666", size = 10),
      legend.position = "bottom"
    )

  if (show_labels && nrow(inst_tab) > 0) {
    y_top <- max(x, na.rm = TRUE) * 1.04
    lbl_df <- do.call(rbind, lapply(seq_len(nrow(inst_tab)), function(i) {
      st <- inst_tab$anomaly_idx[i]
      en <- min(n, st + seg_len - 1L)
      data.frame(x = (st + en) / 2, y = y_top,
                 label     = paste0("P", inst_tab$phenotype_id[i]),
                 phenotype = as.character(inst_tab$phenotype_id[i]),
                 stringsAsFactors = FALSE)
    }))
    p <- p + ggplot2::geom_text(
      data = lbl_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   label = .data$label, colour = .data$phenotype),
      size = 3.5, fontface = "bold", inherit.aes = FALSE
    )
  }
  p
}

#' Plot Canonical Phenotype Template Gallery
#'
#' @description
#' A faceted panel showing the **medoid segment** (canonical representative) of
#' each phenotype, annotated with severity and mean confidence.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param title Character. Plot title.
#'
#' @return A `ggplot2` object.
#' @export
plot_phenotype_gallery <- function(result, title = NULL) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  n_ph  <- result$mapper$n_phenotypes
  ph_df <- result$phenotypes
  pal   <- grDevices::hcl.colors(max(n_ph, 2L), palette = "Dynamic")

  rows <- lapply(seq_len(n_ph), function(ph) {
    seg <- result$canonical_templates[[ph]]
    conf <- mean(result$instance_table$phenotype_confidence[
      result$instance_table$phenotype_id == ph], na.rm = TRUE)
    sev  <- ph_df$severity[ph_df$phenotype_id == ph]
    data.frame(
      time      = seq_along(seg) - floor(length(seg) / 2),
      value     = seg,
      phenotype = factor(paste0("P", ph, ": ", ph_df$name[ph_df$phenotype_id == ph])),
      ph_id     = as.character(ph),
      sev_conf  = sprintf("%s | conf %.0f%%", sev, conf * 100),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  lvls <- unique(df$phenotype)

  names(pal) <- as.character(seq_len(n_ph))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value,
                                    colour = .data$ph_id)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "#AAAAAA", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::facet_wrap(~ phenotype, scales = "free_y") +
    ggplot2::geom_label(
      data = unique(df[, c("phenotype", "ph_id", "sev_conf")]),
      ggplot2::aes(x = -Inf, y = Inf, label = .data$sev_conf, colour = .data$ph_id),
      hjust = -0.05, vjust = 1.4, size = 3.2, fontface = "italic",
      fill = "white", linewidth = 0, inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::labs(
      title    = if (is.null(title)) "Canonical Templates (Medoid Segments)" else title,
      subtitle = "Centred on peak; severity class and mean confidence shown",
      x = "Steps from Peak", y = "Value"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      strip.text    = ggplot2::element_text(face = "bold", size = 10),
      panel.spacing = ggplot2::unit(1.2, "lines")
    )
}

#' Plot Phenotype Feature Profile (Discriminating Features)
#'
#' @description
#' Bar chart of mean standardised feature values per phenotype, restricted to
#' the **n top features** ranked by F-statistic discriminability (from
#' `result$feature_importance`).
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param n_features Integer. Top features to display. Default: 8L.
#'
#' @return A `ggplot2` object.
#' @export
plot_phenotype_profile <- function(result, n_features = 8L) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  fm        <- result$feat_mat
  inst_tab  <- result$instance_table
  n_ph      <- result$mapper$n_phenotypes
  pal       <- grDevices::hcl.colors(max(n_ph, 2L), palette = "Dynamic")

  col_m  <- colMeans(fm, na.rm = TRUE)
  col_s  <- apply(fm, 2, stats::sd); col_s[col_s < 1e-10] <- 1
  fm_std <- sweep(sweep(fm, 2, col_m, "-"), 2, col_s, "/")

  top_feats <- head(result$feature_importance$feature, n_features)

  rows <- lapply(seq_len(n_ph), function(ph) {
    mems <- which(inst_tab$phenotype_id == ph)
    vals <- colMeans(fm_std[mems, top_feats, drop = FALSE], na.rm = TRUE)
    data.frame(phenotype = factor(ph), feature = top_feats, value = vals,
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df$feature <- factor(df$feature, levels = rev(top_feats))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$feature, y = .data$value,
                                    fill = .data$phenotype)) +
    ggplot2::geom_col(position = "dodge", width = 0.7, alpha = 0.88) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "#999999", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = pal,
      labels = paste0("Phenotype ", seq_len(n_ph)),
      name   = "Phenotype"
    ) +
    ggplot2::labs(
      title    = "Phenotype Feature Profiles",
      subtitle = sprintf("Top %d features by F-statistic discriminability (standardised)",
                         n_features),
      x = NULL, y = "Standardised Mean"
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

#' Plot the Mapper Topological Landscape
#'
#' @description
#' Visualises the Mapper graph where each node is a cluster of anomaly instances,
#' edges connect overlapping nodes, and connected components (phenotypes) are
#' colour-coded. Node size is proportional to member count; node label shows
#' the mean stability of members.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param title Character. Plot title.
#'
#' @return A `ggplot2` object.
#' @export
plot_phenotype_landscape <- function(result, title = NULL) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  mp       <- result$mapper
  nodes    <- mp$nodes
  edges    <- mp$edges
  ph_ids   <- mp$phenotype_id
  inst_tab <- result$instance_table
  n_ph     <- mp$n_phenotypes
  n_nd     <- length(nodes)

  pal <- grDevices::hcl.colors(max(n_ph, 2L), palette = "Dynamic")

  node_ph <- vapply(seq_len(n_nd), function(nd) {
    mems <- nodes[[nd]]
    ph   <- ph_ids[mems]
    as.integer(names(sort(table(ph), decreasing = TRUE))[1])
  }, integer(1))

  fm      <- result$feat_mat
  node_coords <- t(vapply(nodes, function(mems) {
    colMeans(fm[mems, 1:2, drop = FALSE], na.rm = TRUE)
  }, numeric(2)))
  for (j in 1:2) {
    rng <- range(node_coords[, j])
    if (diff(rng) > 1e-10)
      node_coords[, j] <- (node_coords[, j] - rng[1]) / diff(rng)
  }

  node_sizes <- vapply(nodes, length, integer(1))
  node_stab  <- vapply(seq_len(n_nd), function(nd) {
    mems <- nodes[[nd]]
    st   <- inst_tab$stability[mems]
    mean(st, na.rm = TRUE)
  }, numeric(1))

  nd_df <- data.frame(
    x    = node_coords[, 1],
    y    = node_coords[, 2],
    size = node_sizes,
    ph   = as.character(node_ph),
    stab = round(node_stab, 2),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot()

  if (nrow(edges) > 0) {
    edge_df <- do.call(rbind, lapply(seq_len(nrow(edges)), function(i) {
      u <- edges[i, 1]; v <- edges[i, 2]
      data.frame(x = nd_df$x[u], xend = nd_df$x[v],
                 y = nd_df$y[u], yend = nd_df$y[v])
    }))
    p <- p + ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = .data$x, xend = .data$xend,
                   y = .data$y, yend = .data$yend),
      colour = "#AAAAAA", linewidth = 0.8, alpha = 0.7
    )
  }

  p <- p +
    ggplot2::geom_point(
      data = nd_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   size = .data$size, colour = .data$ph),
      alpha = 0.85
    ) +
    ggplot2::geom_text(
      data = nd_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   label = sprintf("%.2f", .data$stab)),
      size = 2.8, colour = "white", fontface = "bold"
    ) +
    ggplot2::scale_colour_manual(
      values = pal,
      labels = paste0("Phenotype ", seq_len(n_ph)),
      name   = "Phenotype"
    ) +
    ggplot2::scale_size_continuous(range = c(5, 18), guide = "none") +
    ggplot2::labs(
      title    = if (is.null(title))
        "AnomalyPheno: Mapper Topological Landscape" else title,
      subtitle = sprintf(
        "%d nodes . %d edges . %d phenotype(s)  |  label = mean stability",
        n_nd, nrow(edges), n_ph),
      x = "PC1 projection", y = "PC2 projection"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right"
    )
  p
}

#' @export
plot.AnomPhenoResult <- function(x, type = c("gallery", "landscape",
                                              "timeline", "profile"), ...) {
  type <- match.arg(type)
  switch(type,
    gallery   = plot_phenotype_gallery(x, ...),
    landscape = plot_phenotype_landscape(x, ...),
    timeline  = plot_phenotype_timeline(x, ...),
    profile   = plot_phenotype_profile(x, ...)
  )
}

#' Plot Bootstrap Stability Analysis Results
#'
#' @description
#' Two-panel visualisation of [pheno_bootstrap()] output:
#' \enumerate{
#'   \item **n_phenotypes distribution** -- bar chart showing how often each
#'     phenotype count was recovered across B bootstrap resamples. A tall bar
#'     at the original count indicates structural robustness.
#'   \item **Per-phenotype stability** -- horizontal bar chart showing the
#'     fraction of bootstrap runs where each phenotype's core members stayed
#'     together. Scores above 70% indicate stable phenotypes.
#' }
#'
#' @param boot An `AnomPhenoBootstrap` from [pheno_bootstrap()].
#' @param title Character. Overall plot title. Default: auto-generated.
#'
#' @return A `ggplot2` object (two facets side by side via `patchwork` if
#'   available, otherwise the stability panel alone).
#' @export
plot_bootstrap_stability <- function(boot, title = NULL) {
  if (!inherits(boot, "AnomPhenoBootstrap"))
    stop("'boot' must be an AnomPhenoBootstrap from pheno_bootstrap().")

  dist_df <- data.frame(
    n_ph  = as.integer(names(boot$n_pheno_distribution)),
    count = as.integer(boot$n_pheno_distribution),
    is_original = as.integer(names(boot$n_pheno_distribution)) == boot$n_pheno_original
  )

  p_dist <- ggplot2::ggplot(dist_df,
    ggplot2::aes(x = factor(.data$n_ph), y = .data$count,
                 fill = .data$is_original)) +
    ggplot2::geom_col(width = 0.65, alpha = 0.88) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#1565c0", "FALSE" = "#90a4ae"),
      labels = c("TRUE" = "Original", "FALSE" = "Other"),
      name   = ""
    ) +
    ggplot2::labs(
      title    = "Bootstrap: n_phenotypes Distribution",
      subtitle = sprintf("Mode = %d  |  Original = %d",
                         boot$n_pheno_mode, boot$n_pheno_original),
      x = "Number of Phenotypes", y = "Bootstrap Count"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )

  stab_df <- data.frame(
    phenotype = factor(paste0("P", seq_along(boot$phenotype_stability))),
    stability = boot$phenotype_stability * 100,
    tier      = cut(boot$phenotype_stability * 100,
                    breaks = c(-Inf, 50, 70, 90, Inf),
                    labels = c("Unstable", "Moderate", "Stable", "Very Stable"))
  )

  tier_pal <- c(
    "Unstable"    = "#ef5350",
    "Moderate"    = "#ffa726",
    "Stable"      = "#66bb6a",
    "Very Stable" = "#1565c0"
  )

  p_stab <- ggplot2::ggplot(stab_df,
    ggplot2::aes(x = .data$phenotype, y = .data$stability, fill = .data$tier)) +
    ggplot2::geom_col(width = 0.6, alpha = 0.88) +
    ggplot2::geom_hline(yintercept = 70, linetype = "dashed",
                        colour = "#555555", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", .data$stability)),
                       hjust = -0.2, size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = tier_pal, name = "Stability Tier") +
    ggplot2::scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = "Per-Phenotype Stability",
      subtitle = sprintf("Mean Jaccard = %.3f  |  dashed line = 70%% threshold",
                         boot$mean_jaccard),
      x = NULL, y = "Stability (%)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(p_dist, p_stab, ncol = 2) +
      patchwork::plot_annotation(
        title = if (is.null(title))
          "AnomalyPheno: Bootstrap Stability Analysis" else title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
      )
    return(combined)
  }
  p_stab
}

#' Plot Continuous Anomaly Score Time Series
#'
#' @description
#' Draws the raw time series and overlays a continuous anomaly severity score
#' from `pheno_score()`, coloured by the detected phenotype.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param type Character. Score type for `pheno_score()`. Default: `"combined"`.
#' @param title Character. Plot title.
#'
#' @return A `ggplot2` object.
#' @export
plot_score_timeline <- function(result,
                                 type  = c("combined", "zscore", "confidence"),
                                 title = NULL) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")
  type    <- match.arg(type)
  df      <- as.data.frame(result, score_type = type)
  n_ph    <- result$mapper$n_phenotypes
  ph_keys <- as.character(seq_len(n_ph))
  pal     <- setNames(grDevices::hcl.colors(max(n_ph, 2L), palette = "Dynamic"),
                      ph_keys)

  df$ph_col <- ifelse(is.na(df$phenotype_id), NA_character_,
                      as.character(df$phenotype_id))
  anom_df <- df[df$is_anomaly & !is.na(df$ph_col), ]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$value),
                       colour = "#888888", linewidth = 0.4, alpha = 0.7) +
    ggplot2::geom_col(
      data = anom_df,
      ggplot2::aes(x = .data$time, y = .data$anomaly_score,
                   fill = .data$ph_col),
      width = 1, alpha = 0.75, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = pal,
      labels = paste0("Phenotype ", ph_keys),
      name   = "Phenotype",
      na.value = "transparent"
    ) +
    ggplot2::labs(
      title    = if (is.null(title))
        sprintf("Anomaly Score Timeline [%s]", type) else title,
      subtitle = "Grey line = original series . Bars = anomaly score per segment",
      x = "Time Index", y = sprintf("Anomaly Score (%s)", type)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}


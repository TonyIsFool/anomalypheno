
#' Merge Weak or Unstable Phenotypes
#'
#' @description
#' Post-processes an `AnomPhenoResult` to **merge over-split phenotypes** into
#' their nearest neighbours in the 27-dimensional standardised feature space.
#' A phenotype is considered "weak" if it has fewer than `min_instances`
#' members, or if its mean multi-resolution stability score is below
#' `min_stability`.
#'
#' Weak phenotype instances are reassigned to the closest *strong* phenotype
#' (by Euclidean centroid distance). The result is then re-characterised and
#' re-named using the discriminant-based naming algorithm.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param min_instances Integer. Minimum number of instances for a phenotype
#'   to be retained. Default: `2L`.
#' @param min_stability Numeric in \[0,1\]. Minimum mean stability score.
#'   Default: `0.60`.
#' @param verbose Logical. Print merge summary. Default: `TRUE`.
#'
#' @return A pruned `AnomPhenoResult` with `$pruned = TRUE` and
#'   `$n_merged` reporting how many phenotypes were absorbed.
#'
#' @examples
#' data(anomaly_benchmark)
#' result  <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result) && result$mapper$n_phenotypes > 2) {
#'   pruned <- pheno_prune(result, min_instances = 2L, min_stability = 0.7)
#'   print(pruned)
#' }
#'
#' @export
pheno_prune <- function(result,
                         min_instances = 2L,
                         min_stability = 0.60,
                         verbose       = TRUE) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  inst_tab  <- result$instance_table
  ph_df     <- result$phenotypes
  feat_mat  <- result$feat_mat
  n_ph      <- result$mapper$n_phenotypes
  n_inst    <- nrow(inst_tab)

  cm  <- colMeans(feat_mat, na.rm = TRUE)
  cs  <- apply(feat_mat, 2, stats::sd); cs[cs < 1e-10] <- 1
  fms <- sweep(sweep(feat_mat, 2, cm, "-"), 2, cs, "/")

  cents <- do.call(rbind, lapply(seq_len(n_ph), function(ph) {
    mems <- which(inst_tab$phenotype_id == ph)
    if (length(mems) == 0) return(rep(NA_real_, ncol(fms)))
    colMeans(fms[mems, , drop = FALSE], na.rm = TRUE)
  }))

  ph_cnt  <- tabulate(inst_tab$phenotype_id, nbins = n_ph)
  ph_stab <- vapply(seq_len(n_ph), function(ph) {
    mems <- which(inst_tab$phenotype_id == ph)
    if (length(mems) == 0) return(0)
    mean(inst_tab$stability[mems], na.rm = TRUE)
  }, numeric(1))

  is_weak  <- (ph_cnt < min_instances) | (ph_stab < min_stability)
  is_empty <- ph_cnt == 0L
  is_weak  <- is_weak | is_empty
  strong   <- which(!is_weak)

  if (length(strong) == 0) {
    if (verbose) message("pheno_prune: all phenotypes are weak -- no pruning applied.")
    return(result)
  }
  if (!any(is_weak)) {
    if (verbose) message("pheno_prune: no weak phenotypes found -- result unchanged.")
    result$pruned  <- FALSE
    result$n_merged <- 0L
    return(result)
  }

  n_merged <- sum(is_weak)

  remap <- seq_len(n_ph)
  for (ph in which(is_weak)) {
    dists <- vapply(strong, function(s) {
      cv <- cents[ph, ] - cents[s, ]
      sqrt(sum(cv^2, na.rm = TRUE))
    }, numeric(1))
    remap[ph] <- strong[which.min(dists)]
  }

  new_ids <- remap[inst_tab$phenotype_id]

  ulabs  <- sort(unique(new_ids))
  K_new  <- length(ulabs)
  renum  <- setNames(seq_len(K_new), as.character(ulabs))
  new_ids <- as.integer(renum[as.character(new_ids)])
  inst_tab$phenotype_id <- new_ids

  new_ph_df <- do.call(rbind, lapply(seq_len(K_new), function(ph) {
    mems <- which(new_ids == ph)
    orig_strong <- ulabs[ph]
    old_row     <- ph_df[ph_df$phenotype_id == orig_strong, ]
    if (nrow(old_row) == 0) old_row <- ph_df[1, ]
    data.frame(
      phenotype_id   = ph,
      name           = old_row$name[1],
      count          = length(mems),
      mean_amplitude = mean(abs(inst_tab$peak_abs[mems]), na.rm = TRUE),
      mean_duration  = mean(inst_tab$duration[mems],      na.rm = TRUE),
      mean_zscore    = mean(inst_tab$global_zscore[mems], na.rm = TRUE),
      mean_hurst     = if ("mean_hurst" %in% names(old_row)) old_row$mean_hurst[1] else NA_real_,
      mean_apen      = if ("mean_apen"  %in% names(old_row)) old_row$mean_apen[1]  else NA_real_,
      severity       = old_row$severity[1],
      stringsAsFactors = FALSE
    )
  }))

  new_tmpl <- lapply(seq_len(K_new), function(ph) {
    mems <- which(new_ids == ph)
    sub  <- fms[mems, , drop = FALSE]
    if (nrow(sub) == 1) {
      med_g <- mems[1]
    } else {
      dm      <- as.matrix(stats::dist(sub))
      med_loc <- which.min(rowSums(dm))
      med_g   <- mems[med_loc]
    }
    st <- inst_tab$anomaly_idx[med_g]
    en <- min(length(result$x), st + result$seg_len - 1L)
    result$x[st:en]
  })

  new_ph_df$name <- .discriminant_names(fms, new_ids, K_new)

  result$instance_table      <- inst_tab
  result$phenotypes          <- new_ph_df
  result$canonical_templates <- new_tmpl
  result$mapper$n_phenotypes <- K_new
  result$pruned              <- TRUE
  result$n_merged            <- n_merged

  if (verbose)
    message(sprintf(
      "pheno_prune: %d \u2192 %d phenotypes (%d absorbed).",
      n_ph, K_new, n_merged
    ))
  result
}

#' Rename Phenotypes Using Maximum Discriminant Feature Deviations
#'
#' @description
#' Replaces the phenotype names in an `AnomPhenoResult` with interpretable,
#' **automatically-generated labels** derived from each phenotype's most
#' distinctive features relative to the global centroid.
#'
#' Each phenotype is named by its top two deviating features:
#' \itemize{
#'   \item \strong{Primary}: the feature where the phenotype centroid departs
#'     most (in standardised units) from the cross-phenotype mean.
#'   \item \strong{Secondary}: the next most distinctive feature.
#' }
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()] or `pheno_prune()`.
#'
#' @return The same `AnomPhenoResult` with updated `$phenotypes$name` column.
#'
#' @examples
#' data(anomaly_benchmark)
#' result  <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result)) {
#'   result2 <- pheno_rename(result)
#'   print(result2$phenotypes[, c("phenotype_id", "name")])
#' }
#'
#' @export
pheno_rename <- function(result) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")
  fms    <- result$feat_mat
  ids    <- result$instance_table$phenotype_id
  K      <- result$mapper$n_phenotypes
  cm     <- colMeans(fms, na.rm = TRUE)
  cs     <- apply(fms, 2, stats::sd); cs[cs < 1e-10] <- 1
  fms_st <- sweep(sweep(fms, 2, cm, "-"), 2, cs, "/")
  result$phenotypes$name <- .discriminant_names(fms_st, ids, K)
  result
}

#' Compute Internal Clustering Validation Metrics
#'
#' @description
#' Quantifies the **statistical quality** of the phenotype partition using
#' three complementary internal clustering criteria:
#'
#' \describe{
#'   \item{Silhouette}{Mean s(i) \in \[-1,1\]: measures how well each instance
#'     fits its assigned phenotype relative to the nearest alternative.
#'     Values > 0.5 indicate a strong structure.}
#'   \item{Calinski-Harabász (CH)}{Between-cluster scatter / within-cluster
#'     scatter, corrected for k. Higher is better; no upper bound.}
#'   \item{Davies-Bouldin (DB)}{Mean ratio of intra-cluster spread to
#'     inter-cluster centroid distance. Lower is better; 0 is ideal.}
#' }
#'
#' These three metrics provide triangulated evidence for phenotype quality
#' and are commonly requested by reviewers of unsupervised clustering papers.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()] or `pheno_prune()`.
#'
#' @return An `AnomPhenoQuality` object. Print and plot methods available.
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result)) {
#'   qlt <- pheno_cluster_quality(result)
#'   print(qlt)
#' }
#'
#' @export
pheno_cluster_quality <- function(result) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  fms    <- result$feat_mat
  ids    <- result$instance_table$phenotype_id
  n      <- nrow(fms)
  k      <- result$mapper$n_phenotypes

  if (k < 2L || n < 3L) {
    warning("pheno_cluster_quality: need k >= 2 and n >= 3; returning NULL.")
    return(NULL)
  }

  cm     <- colMeans(fms, na.rm = TRUE)
  cs     <- apply(fms, 2, stats::sd); cs[cs < 1e-10] <- 1
  fms_st <- sweep(sweep(fms, 2, cm, "-"), 2, cs, "/")
  dm     <- as.matrix(stats::dist(fms_st))

  cents <- lapply(seq_len(k), function(ph) {
    mems <- which(ids == ph)
    if (length(mems) == 0) return(rep(NA_real_, ncol(fms_st)))
    colMeans(fms_st[mems, , drop = FALSE], na.rm = TRUE)
  })

  sil <- vapply(seq_len(n), function(i) {
    ph_i  <- ids[i]
    same  <- which(ids == ph_i & seq_len(n) != i)
    if (length(same) == 0) return(0)
    a     <- mean(dm[i, same])
    other <- unique(ids[ids != ph_i])
    if (length(other) == 0) return(0)
    b <- min(vapply(other, function(o)
      mean(dm[i, which(ids == o)]), numeric(1)))
    (b - a) / max(a, b)
  }, numeric(1))

  gm    <- colMeans(fms_st, na.rm = TRUE)
  n_k   <- tabulate(ids, nbins = k)
  SSB   <- sum(vapply(seq_len(k), function(j)
    n_k[j] * sum((cents[[j]] - gm)^2, na.rm = TRUE), numeric(1)))
  SSW   <- sum(vapply(seq_len(n), function(i) {
    c_i <- cents[[ids[i]]]
    sum((fms_st[i, ] - c_i)^2, na.rm = TRUE)
  }, numeric(1)))
  CH <- if (SSW > 0 && k > 1) (SSB / (k - 1)) / (SSW / (n - k)) else NA_real_

  sigma <- vapply(seq_len(k), function(j) {
    mems <- which(ids == j)
    if (length(mems) < 2) return(0)
    mean(vapply(mems, function(i)
      sqrt(sum((fms_st[i, ] - cents[[j]])^2, na.rm = TRUE)), numeric(1)))
  }, numeric(1))

  DB <- mean(vapply(seq_len(k), function(i) {
    vals <- vapply(seq_len(k), function(j) {
      if (i == j) return(-Inf)
      dij <- sqrt(sum((cents[[i]] - cents[[j]])^2, na.rm = TRUE))
      if (dij < 1e-10) return(Inf)
      (sigma[i] + sigma[j]) / dij
    }, numeric(1))
    max(vals[is.finite(vals)], 0)
  }, numeric(1)))

  structure(
    list(
      silhouette_mean         = mean(sil),
      silhouette_per_instance = sil,
      calinski_harabasz       = CH,
      davies_bouldin          = DB,
      n_phenotypes            = k,
      n_instances             = n,
      labels                  = ids
    ),
    class = "AnomPhenoQuality"
  )
}

#' @export
print.AnomPhenoQuality <- function(x, ...) {
  cat("\u2500\u2500 AnomalyPheno: Cluster Quality Report ",
      strrep("\u2500", 30), "\n", sep = "")
  cat(sprintf("  Phenotypes    : %d\n", x$n_phenotypes))
  cat(sprintf("  Instances     : %d\n", x$n_instances))
  cat(sprintf("  Silhouette    : %.4f   [> 0.5 = strong structure]\n",
              x$silhouette_mean))
  cat(sprintf("  Calinski-Harab\u00e1sz : %.2f   [higher is better]\n",
              x$calinski_harabasz))
  cat(sprintf("  Davies-Bouldin: %.4f   [lower is better; 0 = ideal]\n",
              x$davies_bouldin))
  invisible(x)
}

#' @export
plot.AnomPhenoQuality <- function(x, ...) {
  n_ph   <- x$n_phenotypes
  sil    <- x$silhouette_per_instance
  ids    <- x$labels
  pal    <- grDevices::hcl.colors(max(n_ph, 2L), palette = "Dynamic")
  df     <- data.frame(
    instance  = order(ids, sil),
    sil_value = sort(sil[order(ids, sil)]),
    phenotype = factor(ids[order(ids, sil)])
  )
  df$x_pos <- seq_len(nrow(df))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x_pos,
                                    y = .data$sil_value,
                                    fill = .data$phenotype)) +
    ggplot2::geom_col(width = 1, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5,
                        colour = "#333333") +
    ggplot2::geom_hline(yintercept = x$silhouette_mean,
                        linetype = "dashed", colour = "#CC0000", linewidth = 0.7) +
    ggplot2::scale_fill_manual(
      values = setNames(pal, as.character(seq_len(n_ph))),
      name   = "Phenotype",
      labels = paste0("P", seq_len(n_ph))
    ) +
    ggplot2::annotate("text", x = nrow(df) * 0.02,
                      y = x$silhouette_mean + 0.03,
                      label = sprintf("mean = %.3f", x$silhouette_mean),
                      colour = "#CC0000", size = 3.5, hjust = 0) +
    ggplot2::labs(
      title    = "Silhouette Plot: Phenotype Cluster Quality",
      subtitle = sprintf(
        "CH = %.2f  |  DB = %.4f  |  n = %d instances",
        x$calinski_harabasz, x$davies_bouldin, x$n_instances),
      x = "Instances (sorted by phenotype, then silhouette)",
      y = "Silhouette Width s(i)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

#' Track Phenotype Prevalence Over Sliding Windows
#'
#' @description
#' Computes the **fractional prevalence of each phenotype** within overlapping
#' windows of the time series, revealing how the distribution of anomaly
#' types evolves over time.
#'
#' This is particularly useful for:
#' \itemize{
#'   \item Detecting **regime shifts** -- sustained changes in which phenotype
#'     dominates.
#'   \item Monitoring **treatment effects** -- e.g., does the spike phenotype
#'     diminish after a clinical intervention?
#'   \item **Longitudinal studies** -- comparing phenotype distributions across
#'     temporal segments.
#' }
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()] or `pheno_prune()`.
#' @param window_size Integer. Window size in time steps. Default: `100L`.
#' @param step Integer. Step size between windows. Default: `20L`.
#'
#' @return An `AnomPhenoTrack` object with:
#'   \describe{
#'     \item{`track_matrix`}{Numeric matrix (\eqn{W \times K}): phenotype
#'       prevalence proportions per window.}
#'     \item{`window_starts`}{Integer vector of window start indices.}
#'     \item{`window_midpoints`}{Numeric vector of window midpoints.}
#'     \item{`n_phenotypes`}{Number of phenotypes.}
#'   }
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result)) {
#'   trk <- pheno_track(result, window_size = 120L, step = 20L)
#'   plot(trk)
#' }
#'
#' @export
pheno_track <- function(result, window_size = 100L, step = 20L) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")
  window_size <- as.integer(window_size)
  step        <- as.integer(step)
  n           <- length(result$x)
  inst_tab    <- result$instance_table
  n_ph        <- result$mapper$n_phenotypes

  if (window_size >= n)
    stop("'window_size' must be smaller than the series length.")

  starts <- seq(1L, n - window_size + 1L, by = step)
  mids   <- starts + window_size / 2

  track_mat <- do.call(rbind, lapply(starts, function(st) {
    en   <- st + window_size - 1L
    mask <- inst_tab$anomaly_idx >= st & inst_tab$anomaly_idx <= en
    cnts <- tabulate(inst_tab$phenotype_id[mask], nbins = n_ph)
    if (sum(cnts) == 0) return(rep(0, n_ph))
    cnts / sum(cnts)
  }))

  colnames(track_mat) <- paste0("P", seq_len(n_ph))

  structure(
    list(
      track_matrix    = track_mat,
      window_starts   = starts,
      window_midpoints = mids,
      window_size     = window_size,
      step            = step,
      n_phenotypes    = n_ph,
      phenotypes      = result$phenotypes,
      n               = n
    ),
    class = "AnomPhenoTrack"
  )
}

#' @export
print.AnomPhenoTrack <- function(x, ...) {
  cat("\u2500\u2500 AnomalyPheno: Phenotype Evolution Track ",
      strrep("\u2500", 28), "\n", sep = "")
  cat(sprintf("  Windows       : %d  (size = %d, step = %d)\n",
              nrow(x$track_matrix), x$window_size, x$step))
  cat(sprintf("  Phenotypes    : %d\n", x$n_phenotypes))
  cat(sprintf("  Time range    : [1, %d]\n", x$n))
  cat("\nPhenotype prevalence (window x phenotype):\n")
  pm <- round(x$track_matrix * 100, 1)
  rownames(pm) <- sprintf("t=%d", x$window_midpoints)
  print(pm)
  invisible(x)
}

#' @export
plot.AnomPhenoTrack <- function(x, title = NULL, ...) {
  n_ph <- x$n_phenotypes
  pal  <- grDevices::hcl.colors(max(n_ph, 2L), palette = "Dynamic")

  rows <- lapply(seq_len(n_ph), function(ph) {
    data.frame(
      midpoint  = x$window_midpoints,
      prevalence = x$track_matrix[, ph] * 100,
      phenotype = factor(ph),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$midpoint,
                                    y = .data$prevalence,
                                    fill = .data$phenotype)) +
    ggplot2::geom_area(position = "stack", alpha = 0.82) +
    ggplot2::scale_fill_manual(
      values = setNames(pal, as.character(seq_len(n_ph))),
      labels = if (!is.null(x$phenotypes))
        paste0("P", seq_len(n_ph), ": ", x$phenotypes$name) else
        paste0("Phenotype ", seq_len(n_ph)),
      name = "Phenotype"
    ) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(v, "%"),
                                expand  = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(
      title    = if (is.null(title))
        "Phenotype Prevalence Over Time" else title,
      subtitle = sprintf(
        "Sliding window: size = %d, step = %d | Stacked = 100%% of anomalies per window",
        x$window_size, x$step),
      x = "Time Index (window midpoint)", y = "Prevalence (%)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
}

.FEAT_LABEL <- c(
  peak_signed      = "Peak Direction",
  peak_abs         = "Amplitude",
  sharpness        = "Sharpness",
  onset_slope      = "Onset Rate",
  recovery_slope   = "Recovery Rate",
  auc_baseline     = "Area-Under-Curve",
  duration         = "Duration",
  asymmetry        = "Wave Asymmetry",
  dom_freq         = "Dominant Frequency",
  spec_entropy     = "Spectral Entropy",
  hf_energy_ratio  = "HF Energy",
  spec_centroid    = "Spectral Centroid",
  spec_flatness    = "Spectral Flatness",
  pre_mean         = "Pre-Event Level",
  pre_var          = "Pre-Event Variability",
  post_mean        = "Post-Event Level",
  post_var         = "Post-Event Variability",
  local_zscore     = "Local Z-Score",
  global_zscore    = "Global Z-Score",
  pct_rank         = "Severity Rank",
  dev_median       = "Median Deviation",
  iqr_score        = "IQR Score",
  approx_entropy   = "Complexity",
  turning_pt_rate  = "Oscillation Rate",
  trend_slope      = "Trend Slope",
  hurst_approx     = "Persistence",
  zero_cross_rate  = "Zero-Crossing Rate"
)

.discriminant_names <- function(fms_std, ids, K) {
  cents <- do.call(rbind, lapply(seq_len(K), function(ph) {
    mems <- which(ids == ph)
    if (length(mems) == 0) return(rep(0, ncol(fms_std)))
    colMeans(fms_std[mems, , drop = FALSE], na.rm = TRUE)
  }))
  grand <- colMeans(cents, na.rm = TRUE)

  vapply(seq_len(K), function(ph) {
    dev   <- cents[ph, ] - grand
    adev  <- abs(dev)
    ranked <- order(adev, decreasing = TRUE)

    parts <- character(0)
    for (r in seq_len(min(2, length(ranked)))) {
      feat <- colnames(fms_std)[ranked[r]]
      sign <- if (dev[ranked[r]] >= 0) "High" else "Low"
      lbl  <- if (feat %in% names(.FEAT_LABEL)) .FEAT_LABEL[[feat]] else feat
      parts <- c(parts, sprintf("%s %s", sign, lbl))
    }
    paste(parts, collapse = " | ")
  }, character(1))
}

#' Auto-Tune Detection Threshold via Grid Search
#'
#' @description
#' Searches a grid of `q` values, evaluates each using internal clustering
#' quality (Silhouette, CH, DB), and returns the configuration that best
#' separates phenotypes. Eliminates manual parameter selection.
#'
#' @param x Numeric vector. Time series.
#' @param q_grid Numeric vector. Candidate threshold values. Default:
#'   `c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0)`.
#' @param engine Character. Detection engine: `"sr"` or `"lof"`.
#' @param criterion Character. Optimisation criterion: `"silhouette"`,
#'   `"ch"` (Calinski-Harabász), or `"composite"` (average of all three).
#' @param min_instances Integer. Minimum anomaly instances required for a `q`
#'   to be scored. Default: `3L`.
#' @param verbose Logical. Show progress. Default: `TRUE`.
#'
#' @return An `AnomPhenoTune` object with fields `best_q`, `best_result`,
#'   `metrics` (data frame), and `criterion`.
#'
#' @examples
#' data(anomaly_benchmark)
#' tune <- suppressMessages(
#'   pheno_auto_tune(anomaly_benchmark$x,
#'                   q_grid = c(2.0, 2.5, 3.0), criterion = "composite"))
#' print(tune)
#'
#' @export
pheno_auto_tune <- function(x,
                              q_grid        = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0),
                              engine        = c("sr", "lof"),
                              criterion     = c("composite", "silhouette", "ch"),
                              min_instances = 3L,
                              verbose       = TRUE) {
  engine    <- match.arg(engine)
  criterion <- match.arg(criterion)

  metrics <- data.frame(
    q                 = q_grid,
    n_detected        = NA_integer_,
    n_phenotypes      = NA_integer_,
    silhouette        = NA_real_,
    calinski_harabasz = NA_real_,
    davies_bouldin    = NA_real_,
    stringsAsFactors  = FALSE
  )
  stored <- vector("list", length(q_grid))

  for (i in seq_along(q_grid)) {
    qi <- q_grid[i]
    if (verbose) message(sprintf("  pheno_auto_tune: q = %.2f ...", qi))
    ri <- suppressMessages(tryCatch(
      pheno_diagnose(x, q = qi, engine = engine),
      error = function(e) NULL
    ))
    n_det <- if (is.null(ri)) 0L else nrow(ri$instance_table)
    n_ph  <- if (is.null(ri)) 0L else ri$mapper$n_phenotypes
    metrics$n_detected[i]   <- n_det
    metrics$n_phenotypes[i] <- n_ph
    if (is.null(ri) || n_det < min_instances || n_ph < 2L) next
    qlt <- pheno_cluster_quality(ri)
    if (is.null(qlt)) next
    metrics$silhouette[i]        <- qlt$silhouette_mean
    metrics$calinski_harabasz[i] <- qlt$calinski_harabasz
    metrics$davies_bouldin[i]    <- qlt$davies_bouldin
    stored[[i]] <- ri
  }

  .safe_norm <- function(v) {
    rng <- range(v, na.rm = TRUE)
    if (diff(rng) < 1e-10) return(ifelse(is.na(v), NA_real_, 0.5))
    (v - rng[1]) / diff(rng)
  }

  score <- switch(criterion,
    silhouette = metrics$silhouette,
    ch         = metrics$calinski_harabasz,
    composite  = {
      s  <- .safe_norm(metrics$silhouette)
      ch <- .safe_norm(metrics$calinski_harabasz)
      db <- 1 - .safe_norm(metrics$davies_bouldin)
      rowMeans(cbind(s, ch, db), na.rm = TRUE)
    }
  )

  valid_idx <- which(!is.na(score))
  if (length(valid_idx) == 0) {
    warning("pheno_auto_tune: no valid q found; returning first non-NULL result.")
    non_null <- which(!vapply(stored, is.null, logical(1)))
    best_idx <- if (length(non_null) > 0) non_null[1] else 1L
  } else {
    best_idx <- valid_idx[which.max(score[valid_idx])]
  }

  if (verbose) message(sprintf("  Best q = %.2f (criterion: %s)",
                                q_grid[best_idx], criterion))
  structure(
    list(metrics     = metrics,
         best_q      = q_grid[best_idx],
         best_result = stored[[best_idx]],
         criterion   = criterion,
         q_grid      = q_grid,
         engine      = engine),
    class = "AnomPhenoTune"
  )
}

#' @export
print.AnomPhenoTune <- function(x, ...) {
  cat("\u2500\u2500 AnomalyPheno Auto-Tune Report ", strrep("\u2500", 32), "\n", sep = "")
  cat(sprintf("  Best q        : %.2f\n", x$best_q))
  cat(sprintf("  Criterion     : %s\n",   x$criterion))
  cat(sprintf("  Engine        : %s\n",   x$engine))
  cat("\nMetrics across q grid:\n")
  m <- x$metrics
  m[, 4:6] <- round(m[, 4:6], 4)
  print(m, row.names = FALSE)
  invisible(x)
}

#' @export
plot.AnomPhenoTune <- function(x, ...) {
  df <- x$metrics
  df_long <- do.call(rbind, list(
    data.frame(q = df$q, value = df$silhouette,
               metric = "Silhouette", stringsAsFactors = FALSE),
    data.frame(q = df$q, value = df$calinski_harabasz / max(df$calinski_harabasz, na.rm=TRUE),
               metric = "CH (scaled)", stringsAsFactors = FALSE),
    data.frame(q = df$q, value = 1 - df$davies_bouldin / max(df$davies_bouldin, na.rm=TRUE),
               metric = "1 - DB (scaled)", stringsAsFactors = FALSE)
  ))
  ggplot2::ggplot(df_long, ggplot2::aes(x = .data$q, y = .data$value,
                                         colour = .data$metric, group = .data$metric)) +
    ggplot2::geom_line(linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = x$best_q, linetype = "dashed",
                        colour = "#CC0000", linewidth = 0.8) +
    ggplot2::annotate("text", x = x$best_q, y = Inf,
                      label = sprintf(" q* = %.2f", x$best_q),
                      colour = "#CC0000", vjust = 1.5, hjust = 0, size = 3.5) +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Metric") +
    ggplot2::labs(
      title    = "Auto-Tune: Cluster Quality vs Detection Threshold",
      subtitle = sprintf("Criterion: %s  |  Engine: %s", x$criterion, x$engine),
      x = "Detection Threshold q", y = "Quality Score"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}

#' Compute Phenotype-to-Phenotype Markov Transition Matrix
#'
#' @description
#' Treats the temporally ordered sequence of anomaly phenotypes as a
#' discrete-time Markov chain and estimates the **transition probability
#' matrix** P where P\[i,j\] = Prob(next anomaly is phenotype j | current
#' anomaly is phenotype i).
#'
#' Also computes the **stationary distribution** (long-run prevalence of each
#' phenotype) via power iteration, revealing which phenotype the system
#' gravitates towards over time.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#'
#' @return An `AnomPhenoTransition` object with:
#'   \describe{
#'     \item{`transition_matrix`}{K \eqn{\times} K row-stochastic probability matrix.}
#'     \item{`transition_counts`}{Raw count matrix.}
#'     \item{`stationary_dist`}{Stationary distribution vector (length K).}
#'     \item{`sequence`}{Time-ordered phenotype ID sequence.}
#'   }
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result)) {
#'   tr <- pheno_transition(result)
#'   print(tr)
#' }
#'
#' @export
pheno_transition <- function(result) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")
  inst_tab <- result$instance_table
  n_ph     <- result$mapper$n_phenotypes

  ord  <- order(inst_tab$anomaly_idx)
  seqs <- inst_tab$phenotype_id[ord]

  if (length(seqs) < 2L) {
    warning("pheno_transition: need at least 2 anomaly instances.")
    return(NULL)
  }

  cnt <- matrix(0L, nrow = n_ph, ncol = n_ph,
                dimnames = list(paste0("P", seq_len(n_ph)),
                                paste0("P", seq_len(n_ph))))
  for (i in seq_len(length(seqs) - 1L))
    cnt[seqs[i], seqs[i + 1L]] <- cnt[seqs[i], seqs[i + 1L]] + 1L

  row_s <- pmax(rowSums(cnt), 1L)
  prob  <- cnt / row_s

  pi_v <- rep(1 / n_ph, n_ph)
  for (it in seq_len(1000L)) {
    pi_n <- as.numeric(pi_v %*% prob)
    if (max(abs(pi_n - pi_v)) < 1e-10) break
    pi_v <- pi_n
  }

  structure(
    list(transition_matrix = prob,
         transition_counts = cnt,
         stationary_dist   = pi_v,
         sequence          = seqs,
         n_phenotypes      = n_ph,
         phenotypes        = result$phenotypes),
    class = "AnomPhenoTransition"
  )
}

#' @export
print.AnomPhenoTransition <- function(x, ...) {
  cat("\u2500\u2500 AnomalyPheno Phenotype Transition Matrix ",
      strrep("\u2500", 22), "\n", sep = "")
  cat(sprintf("  Phenotypes  : %d   |   Transitions observed: %d\n",
              x$n_phenotypes, length(x$sequence) - 1L))
  cat("\nTransition probabilities (row = From, col = To):\n")
  print(round(x$transition_matrix, 3))
  cat("\nStationary distribution:\n")
  sv <- setNames(round(x$stationary_dist, 4), paste0("P", seq_len(x$n_phenotypes)))
  print(sv)
  invisible(x)
}

#' @export
plot.AnomPhenoTransition <- function(x, ...) {
  n_ph <- x$n_phenotypes
  mat  <- x$transition_matrix
  df   <- expand.grid(From = paste0("P", seq_len(n_ph)),
                      To   = paste0("P", seq_len(n_ph)),
                      stringsAsFactors = FALSE)
  df$prob <- as.vector(t(mat))

  pal <- grDevices::hcl.colors(9, palette = "YlOrRd")

  sd_df <- data.frame(
    phenotype = paste0("P", seq_len(n_ph)),
    stat      = round(x$stationary_dist * 100, 1)
  )

  p_heat <- ggplot2::ggplot(df, ggplot2::aes(x = .data$To, y = .data$From,
                                               fill = .data$prob)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$prob)),
                       size = 3.5, fontface = "bold",
                       colour = ifelse(df$prob > 0.5, "white", "#333333")) +
    ggplot2::scale_fill_gradientn(
      colours = pal, limits = c(0, 1),
      name = "Transition\nProbability"
    ) +
    ggplot2::labs(
      title    = "Phenotype Markov Transition Matrix",
      subtitle = paste0("Stationary dist: ",
                        paste(sprintf("P%d=%.1f%%", seq_len(n_ph),
                                      x$stationary_dist * 100),
                              collapse = "  ")),
      x = "To Phenotype", y = "From Phenotype"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      axis.text   = ggplot2::element_text(face = "bold"),
      panel.grid  = ggplot2::element_blank()
    )
  p_heat
}

#' Detect Phenotype Regime Changes via CUSUM
#'
#' @description
#' Applies a **tabular CUSUM test** to the binary phenotype indicator series
#' to detect sustained shifts in the prevalence of a target phenotype.
#' Returns the time indices of detected regime changes.
#'
#' The CUSUM statistic accumulates deviations from the reference mean:
#' \deqn{C^+_i = \max(0,\ C^+_{i-1} + (y_i - \mu_0) - k)}
#' \deqn{C^-_i = \max(0,\ C^-_{i-1} - (y_i - \mu_0) - k)}
#' A regime change is signalled when \eqn{C^+ > h} or \eqn{C^- > h}.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param target_ph Integer. Phenotype ID to monitor. `NULL` = most frequent.
#' @param h Numeric. CUSUM decision threshold. Default: `3`.
#' @param k Numeric. Reference slack (typically `0.5`). Default: `0.5`.
#'
#' @return An `AnomPhenoRegime` object.
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result)) {
#'   reg <- pheno_regime(result, h = 2)
#'   print(reg)
#' }
#'
#' @export
pheno_regime <- function(result, target_ph = NULL, h = 3, k = 0.5) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")
  inst_tab <- result$instance_table
  n_ph     <- result$mapper$n_phenotypes

  ord   <- order(inst_tab$anomaly_idx)
  times <- inst_tab$anomaly_idx[ord]
  ids   <- inst_tab$phenotype_id[ord]

  if (length(ids) < 3L) {
    warning("pheno_regime: need at least 3 anomaly instances.")
    return(NULL)
  }

  if (is.null(target_ph))
    target_ph <- as.integer(names(sort(table(ids), decreasing = TRUE))[1])

  y   <- as.integer(ids == target_ph)
  mu0 <- mean(y)
  n   <- length(y)

  C_pos <- numeric(n)
  C_neg <- numeric(n)
  for (i in seq_len(n)) {
    z        <- y[i] - mu0
    prev_p   <- if (i == 1L) 0 else C_pos[i - 1L]
    prev_n   <- if (i == 1L) 0 else C_neg[i - 1L]
    C_pos[i] <- max(0, prev_p + z - k)
    C_neg[i] <- max(0, prev_n - z - k)
  }

  sigs <- which(C_pos > h | C_neg > h)
  changepoints <- integer(0)
  if (length(sigs) > 0L) {
    grps <- cumsum(c(1L, diff(sigs) > 2L))
    for (g in unique(grps))
      changepoints <- c(changepoints, sigs[grps == g][1L])
    changepoints <- times[changepoints]
  }

  structure(
    list(changepoints = changepoints,
         cusum_pos    = C_pos,
         cusum_neg    = C_neg,
         threshold    = h,
         target_ph    = target_ph,
         times        = times,
         y            = y,
         n_phenotypes = n_ph,
         phenotypes   = result$phenotypes),
    class = "AnomPhenoRegime"
  )
}

#' @export
print.AnomPhenoRegime <- function(x, ...) {
  cat("\u2500\u2500 AnomalyPheno CUSUM Regime Detection ",
      strrep("\u2500", 28), "\n", sep = "")
  cat(sprintf("  Target phenotype  : P%d\n", x$target_ph))
  cat(sprintf("  CUSUM threshold h : %.2f\n", x$threshold))
  cat(sprintf("  Instances monitored: %d\n", length(x$y)))
  if (length(x$changepoints) == 0) {
    cat("  Regime changes    : none detected\n")
  } else {
    cat(sprintf("  Regime changes    : %d detected at time index(es): %s\n",
                length(x$changepoints),
                paste(x$changepoints, collapse = ", ")))
  }
  invisible(x)
}

#' @export
plot.AnomPhenoRegime <- function(x, ...) {
  n  <- length(x$times)
  df <- data.frame(
    idx     = seq_len(n),
    time    = x$times,
    y       = x$y,
    C_pos   = x$cusum_pos,
    C_neg   = x$cusum_neg
  )

  cp_df <- if (length(x$changepoints) > 0)
    data.frame(time = x$changepoints) else NULL

  pA <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_step(ggplot2::aes(y = .data$y), colour = "#1565c0",
                       linewidth = 1.0) +
    ggplot2::scale_y_continuous(breaks = c(0, 1),
                                labels = c("Other", paste0("P", x$target_ph))) +
    ggplot2::labs(title = sprintf("CUSUM Regime Detection  (target = P%d)", x$target_ph),
                  x = NULL, y = "Phenotype") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  pB <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$C_pos), colour = "#e53935",
                       linewidth = 0.9) +
    ggplot2::geom_line(ggplot2::aes(y = .data$C_neg), colour = "#1e88e5",
                       linewidth = 0.9, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = x$threshold, linetype = "dotted",
                        colour = "#333333", linewidth = 0.7) +
    ggplot2::annotate("text", x = min(df$time), y = x$threshold,
                      label = sprintf(" h = %.1f", x$threshold),
                      colour = "#333333", size = 3, vjust = -0.5, hjust = 0) +
    ggplot2::labs(x = "Time Index (anomaly position)",
                  y = "CUSUM Statistic",
                  subtitle = "Red = C\u207a (increase)  |  Blue dashed = C\u207b (decrease)") +
    ggplot2::theme_minimal(base_size = 10)

  if (!is.null(cp_df)) {
    pA <- pA + ggplot2::geom_vline(data = cp_df,
                                    ggplot2::aes(xintercept = .data$time),
                                    colour = "#e65100", linetype = "dashed",
                                    linewidth = 0.8)
    pB <- pB + ggplot2::geom_vline(data = cp_df,
                                    ggplot2::aes(xintercept = .data$time),
                                    colour = "#e65100", linetype = "dashed",
                                    linewidth = 0.8)
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(pA, pB, ncol = 1, heights = c(1, 2))
  } else pB
}

#' Plot Phenotype Feature Fingerprint as Radar Chart
#'
#' @description
#' Draws a polar coordinate radar (spider) chart showing the **standardised
#' mean value of the top features** for each phenotype. Unlike the bar chart
#' in [plot_phenotype_profile()], the radar chart makes it visually
#' intuitive to compare the overall *shape* of each phenotype's profile.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param n_features Integer. Number of features to include. Default: `8L`.
#' @param title Character. Plot title.
#'
#' @return A `ggplot2` object.
#' @export
plot_phenotype_radar <- function(result, n_features = 8L, title = NULL) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  fm       <- result$feat_mat
  inst_tab <- result$instance_table
  n_ph     <- result$mapper$n_phenotypes
  pal      <- grDevices::hcl.colors(max(n_ph, 2L), palette = "Dynamic")

  cm  <- colMeans(fm, na.rm = TRUE)
  cs  <- apply(fm, 2, stats::sd); cs[cs < 1e-10] <- 1
  fms <- sweep(sweep(fm, 2, cm, "-"), 2, cs, "/")

  top_feats <- head(result$feature_importance$feature, n_features)
  n_f       <- length(top_feats)

  rows <- lapply(seq_len(n_ph), function(ph) {
    mems <- which(inst_tab$phenotype_id == ph)
    vals <- colMeans(fms[mems, top_feats, drop = FALSE], na.rm = TRUE)
    vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals) + 1e-10)
    angs  <- seq(0, 2 * pi, length.out = n_f + 1)[-(n_f + 1)]
    feats <- c(top_feats, top_feats[1])
    vs    <- c(vals_scaled, vals_scaled[1])
    as_   <- c(angs, angs[1])
    data.frame(
      x         = vs * cos(as_),
      y         = vs * sin(as_),
      feature   = feats,
      angle     = as_,
      phenotype = factor(ph),
      stringsAsFactors = FALSE
    )
  })
  df_poly <- do.call(rbind, rows)

  grid_df <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1.0), function(r) {
    th <- seq(0, 2 * pi, length.out = 100)
    data.frame(x = r * cos(th), y = r * sin(th), r = r)
  }))

  angs  <- seq(0, 2 * pi, length.out = n_f + 1)[-(n_f + 1)]
  spoke_df <- data.frame(
    x0 = 0, y0 = 0,
    x1 = cos(angs), y1 = sin(angs),
    label = top_feats,
    lx = 1.22 * cos(angs), ly = 1.22 * sin(angs)
  )

  ggplot2::ggplot() +
    ggplot2::geom_path(data = grid_df,
                       ggplot2::aes(x = .data$x, y = .data$y,
                                    group = .data$r),
                       colour = "#CCCCCC", linewidth = 0.4) +
    ggplot2::geom_segment(data = spoke_df,
                          ggplot2::aes(x = .data$x0, y = .data$y0,
                                       xend = .data$x1, yend = .data$y1),
                          colour = "#AAAAAA", linewidth = 0.5) +
    ggplot2::geom_text(data = spoke_df,
                       ggplot2::aes(x = .data$lx, y = .data$ly,
                                    label = .data$label),
                       size = 2.8, colour = "#444444") +
    ggplot2::geom_polygon(data = df_poly,
                          ggplot2::aes(x = .data$x, y = .data$y,
                                       fill = .data$phenotype,
                                       colour = .data$phenotype,
                                       group = .data$phenotype),
                          alpha = 0.25, linewidth = 1.1) +
    ggplot2::scale_fill_manual(
      values = setNames(pal, as.character(seq_len(n_ph))),
      labels = paste0("Phenotype ", seq_len(n_ph)), name = "Phenotype"
    ) +
    ggplot2::scale_colour_manual(
      values = setNames(pal, as.character(seq_len(n_ph))), guide = "none"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title    = if (is.null(title)) "Phenotype Feature Fingerprints (Radar)" else title,
      subtitle = sprintf("Top %d discriminating features (standardised, scaled to [0,1])", n_f)
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "#666666"),
      legend.position = "right"
    )
}

#' Gaussian Mixture Model Phenotyping with Automatic Model Selection
#'
#' @description
#' An alternative to the Mapper-based phenotyping engine that fits a
#' **Gaussian Mixture Model** (diagonal covariance) to the standardised
#' 27-dimensional feature space via the EM algorithm. The number of
#' components K is chosen automatically using the **Bayesian Information
#' Criterion** (BIC), balancing goodness-of-fit against model complexity.
#'
#' Unlike the hard Mapper assignments, GMM returns **soft posterior
#' probabilities** P(phenotype k | instance i), enabling uncertainty-aware
#' downstream analysis.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param k_max Integer. Maximum number of mixture components. Default: `6L`.
#' @param n_init Integer. Random initialisations per k. Default: `5L`.
#' @param max_iter Integer. Maximum EM iterations per fit. Default: `150L`.
#' @param tol Numeric. Log-likelihood convergence tolerance. Default: `1e-6`.
#' @param seed Integer. Random seed for reproducibility.
#' @param verbose Logical. Print BIC table. Default: `TRUE`.
#'
#' @return An `AnomPhenoGMM` object with fields:
#'   \describe{
#'     \item{`labels`}{Hard phenotype assignments (MAP estimate).}
#'     \item{`posterior`}{n \eqn{\times} K soft probability matrix.}
#'     \item{`mu`}{Component mean vectors (K \eqn{\times} p).}
#'     \item{`best_k`}{Selected number of components.}
#'     \item{`bic_per_k`}{BIC value for each k in 1..k_max.}
#'   }
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result) && nrow(result$instance_table) >= 4) {
#'   gmm <- suppressMessages(pheno_gmm(result, k_max = 4L, n_init = 2L, seed = 42))
#'   print(gmm)
#' }
#'
#' @export
pheno_gmm <- function(result, k_max = 6L, n_init = 5L, max_iter = 150L,
                       tol = 1e-6, seed = NULL, verbose = TRUE) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  fm  <- result$feat_mat
  cm  <- colMeans(fm, na.rm = TRUE)
  cs  <- apply(fm, 2, stats::sd); cs[cs < 1e-10] <- 1
  fms <- sweep(sweep(fm, 2, cm, "-"), 2, cs, "/")
  n   <- nrow(fms); p <- ncol(fms)

  if (!is.null(seed)) set.seed(seed)

  .log_dens <- function(mat, mu_j, sig_j) {
    sig_j[sig_j < 1e-6] <- 1e-6
    -0.5 * p * log(2 * pi) - 0.5 * sum(log(sig_j)) -
      0.5 * rowSums(sweep(mat, 2, mu_j)^2 /
                      matrix(sig_j, nrow = n, ncol = p, byrow = TRUE))
  }

  best_bic <- Inf; best_gmm <- NULL; best_k <- 1L
  bic_vals <- rep(Inf, k_max)

  for (k in seq_len(k_max)) {
    best_ll_k <- -Inf; best_fit_k <- NULL

    for (init in seq_len(n_init)) {
      km <- tryCatch(
        stats::kmeans(fms, centers = k, nstart = 1L, iter.max = 20L),
        error = function(e) NULL)
      if (is.null(km)) next

      mu    <- km$centers
      sigma <- lapply(seq_len(k), function(j) rep(1, p))
      pi_w  <- tabulate(km$cluster, nbins = k) / n
      pi_w  <- pmax(pi_w, 1e-4); pi_w <- pi_w / sum(pi_w)
      ll_prev <- -Inf; gamma <- matrix(0, n, k)

      for (iter in seq_len(max_iter)) {
        lmat <- vapply(seq_len(k), function(j)
          log(pi_w[j]) + .log_dens(fms, mu[j, ], sigma[[j]]),
          numeric(n))
        lse  <- apply(lmat, 1, function(r) { mx <- max(r); mx + log(sum(exp(r - mx))) })
        ll   <- sum(lse)
        gamma <- exp(lmat - lse)
        if (abs(ll - ll_prev) < tol) break
        ll_prev <- ll
        nk <- colSums(gamma); nk[nk < 1e-6] <- 1e-6
        pi_w <- nk / n
        for (j in seq_len(k)) {
          mu[j, ]   <- colSums(gamma[, j] * fms) / nk[j]
          v         <- colSums(gamma[, j] * sweep(fms, 2, mu[j, ])^2) / nk[j]
          sigma[[j]] <- pmax(v, 1e-6)
        }
      }
      if (ll_prev > best_ll_k) {
        best_ll_k  <- ll_prev
        best_fit_k <- list(mu = mu, sigma = sigma, pi = pi_w,
                           gamma = gamma, ll = ll_prev)
      }
    }
    if (is.null(best_fit_k)) next
    n_par    <- k * p * 2 + (k - 1)
    bic      <- -2 * best_fit_k$ll + n_par * log(n)
    bic_vals[k] <- bic
    if (verbose) message(sprintf("  GMM k=%d: BIC=%.2f", k, bic))
    if (bic < best_bic) { best_bic <- bic; best_gmm <- best_fit_k; best_k <- k }
  }

  if (is.null(best_gmm)) stop("GMM fitting failed for all k.")
  labels <- apply(best_gmm$gamma, 1, which.max)

  structure(
    list(labels       = labels,
         posterior    = best_gmm$gamma,
         mu           = best_gmm$mu,
         sigma        = best_gmm$sigma,
         pi           = best_gmm$pi,
         log_lik      = best_gmm$ll,
         bic          = best_bic,
         bic_per_k    = bic_vals,
         best_k       = best_k,
         k_max        = k_max,
         feat_mat_std = fms,
         n_instances  = n),
    class = "AnomPhenoGMM"
  )
}

#' @export
print.AnomPhenoGMM <- function(x, ...) {
  cat("\u2500\u2500 AnomalyPheno Gaussian Mixture Model ",
      strrep("\u2500", 28), "\n", sep = "")
  cat(sprintf("  Optimal K     : %d components (BIC-selected)\n", x$best_k))
  cat(sprintf("  Log-likelihood: %.4f\n", x$log_lik))
  cat(sprintf("  BIC           : %.4f\n", x$bic))
  cat("\nBIC per k:\n")
  bdf <- data.frame(k = seq_len(x$k_max), BIC = round(x$bic_per_k, 3))
  bdf$selected <- ifelse(bdf$k == x$best_k, "<--", "")
  print(bdf, row.names = FALSE)
  cat("\nMixing weights:\n")
  mw <- setNames(round(x$pi, 4), paste0("C", seq_len(x$best_k)))
  print(mw)
  cat("\nPosterior uncertainty (entropy per instance):\n")
  ent <- -rowSums(x$posterior * log(x$posterior + 1e-10))
  cat(sprintf("  Mean: %.3f  |  Max: %.3f  |  Min: %.3f\n",
              mean(ent), max(ent), min(ent)))
  invisible(x)
}

#' @export
plot.AnomPhenoGMM <- function(x, ...) {
  pal <- grDevices::hcl.colors(max(x$best_k, 2L), palette = "Dynamic")
  k   <- x$best_k

  bdf <- data.frame(k = seq_len(x$k_max),
                    BIC = x$bic_per_k,
                    is_best = seq_len(x$k_max) == x$best_k)
  pA <- ggplot2::ggplot(bdf, ggplot2::aes(x = .data$k, y = .data$BIC)) +
    ggplot2::geom_line(colour = "#5c85d6", linewidth = 1.1, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$is_best), size = 4,
                        na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "#5c85d6", "TRUE" = "#e53935"),
                                 guide = "none") +
    ggplot2::annotate("text", x = x$best_k, y = min(x$bic_per_k, na.rm = TRUE),
                      label = sprintf(" K*=%d", x$best_k),
                      colour = "#e53935", size = 3.5, hjust = 0, vjust = 1.5) +
    ggplot2::labs(title  = "GMM Model Selection (BIC)",
                  x = "Number of Components K", y = "BIC") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ent <- -rowSums(x$posterior * log(x$posterior + 1e-10))
  edf <- data.frame(instance  = seq_len(x$n_instances),
                    entropy   = ent,
                    component = factor(x$labels))
  pB <- ggplot2::ggplot(edf, ggplot2::aes(x = .data$instance, y = .data$entropy,
                                           fill = .data$component)) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::scale_fill_manual(
      values = setNames(pal, as.character(seq_len(k))),
      name   = "Component"
    ) +
    ggplot2::labs(x = "Instance Index", y = "Assignment Entropy",
                  subtitle = "Lower = more certain assignment") +
    ggplot2::theme_minimal(base_size = 10)

  if (requireNamespace("patchwork", quietly = TRUE))
    patchwork::wrap_plots(pA, pB, ncol = 2)
  else pA
}

#' Functional Principal Component Analysis of Anomaly Waveforms
#'
#' @description
#' Treats each detected anomaly segment as a **functional object** and
#' decomposes the space of anomaly shapes via Functional Principal Component
#' Analysis (FPCA).
#'
#' The k-th FPC captures the k-th dominant **mode of shape variation** across
#' all anomaly instances. FPC scores (projections) may be used as additional
#' features for downstream clustering or regression.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param n_comp Integer. Number of functional principal components. Default: `3L`.
#' @param align_length Integer. Common segment length after linear interpolation.
#'   Default: `result$seg_len`.
#'
#' @return An `AnomPhenoFPCA` object.
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result)) {
#'   fp <- pheno_fpca(result, n_comp = 3L)
#'   print(fp)
#'   plot(fp)
#' }
#'
#' @export
pheno_fpca <- function(result, n_comp = 3L, align_length = NULL) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  inst_tab <- result$instance_table
  n_inst   <- nrow(inst_tab)
  if (n_inst < 3L) {
    warning("pheno_fpca: need at least 3 instances.")
    return(NULL)
  }

  seg_len <- result$seg_len
  if (is.null(align_length)) align_length <- seg_len

  segs <- do.call(rbind, lapply(inst_tab$anomaly_idx, function(st) {
    en  <- min(length(result$x), st + seg_len - 1L)
    seg <- result$x[st:en]
    if (length(seg) != align_length)
      seg <- stats::approx(seq_along(seg), seg, n = align_length)$y
    seg
  }))

  mean_fn <- colMeans(segs, na.rm = TRUE)
  segs_c  <- sweep(segs, 2, mean_fn, "-")

  sv     <- svd(segs_c / sqrt(max(n_inst - 1, 1)))
  n_comp <- min(n_comp, nrow(sv$u), length(sv$d))

  scores      <- segs_c %*% sv$v[, seq_len(n_comp), drop = FALSE]
  var_exp     <- sv$d^2 / sum(sv$d^2)

  structure(
    list(eigenfunctions = sv$v[, seq_len(n_comp), drop = FALSE],
         scores         = scores[, seq_len(n_comp), drop = FALSE],
         eigenvalues    = sv$d[seq_len(n_comp)]^2,
         var_explained  = var_exp[seq_len(n_comp)],
         mean_function  = mean_fn,
         align_length   = align_length,
         n_comp         = n_comp,
         n_instances    = n_inst,
         labels         = inst_tab$phenotype_id,
         phenotypes     = result$phenotypes),
    class = "AnomPhenoFPCA"
  )
}

#' @export
print.AnomPhenoFPCA <- function(x, ...) {
  cat("\u2500\u2500 AnomalyPheno Functional PCA ",
      strrep("\u2500", 38), "\n", sep = "")
  cat(sprintf("  Instances     : %d\n", x$n_instances))
  cat(sprintf("  Align length  : %d time steps\n", x$align_length))
  cat(sprintf("  Components    : %d\n", x$n_comp))
  cat("\nVariance explained:\n")
  for (j in seq_len(x$n_comp))
    cat(sprintf("  FPC%d: %.2f%%\n", j, x$var_explained[j] * 100))
  cat(sprintf("  Cumulative: %.2f%%\n", sum(x$var_explained) * 100))
  invisible(x)
}

#' @export
plot.AnomPhenoFPCA <- function(x, comp = 1L, scale = 1.5, ...) {
  comp <- min(comp, x$n_comp)
  t_   <- seq_len(x$align_length)
  ef   <- x$eigenfunctions[, comp]
  mf   <- x$mean_function

  ampl <- scale * sqrt(x$eigenvalues[comp])
  df_shape <- data.frame(
    t     = rep(t_, 3),
    value = c(mf, mf + ampl * ef, mf - ampl * ef),
    curve = rep(c("Mean", "+FPC", "-FPC"), each = x$align_length)
  )
  df_shape$curve <- factor(df_shape$curve, levels = c("+FPC", "Mean", "-FPC"))

  pal <- grDevices::hcl.colors(max(length(unique(x$labels)), 2L),
                                palette = "Dynamic")
  n_ph <- length(unique(x$labels))

  pA <- ggplot2::ggplot(df_shape,
                        ggplot2::aes(x = .data$t, y = .data$value,
                                     colour = .data$curve,
                                     linetype = .data$curve,
                                     linewidth = .data$curve)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_manual(
      values = c("Mean" = "#333333", "+FPC" = "#e53935", "-FPC" = "#1565c0"),
      name = NULL) +
    ggplot2::scale_linetype_manual(
      values = c("Mean" = "solid", "+FPC" = "dashed", "-FPC" = "dashed"),
      name = NULL) +
    ggplot2::scale_linewidth_manual(
      values = c("Mean" = 1.2, "+FPC" = 0.9, "-FPC" = 0.9), name = NULL) +
    ggplot2::labs(
      title    = sprintf("FPC%d  (%.1f%% variance explained)", comp,
                         x$var_explained[comp] * 100),
      subtitle = sprintf("Mean \u00b1 %.1f \u00d7 \u03c3(FPC%d)", scale, comp),
      x = "Steps from anomaly onset", y = "Value") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")

  if (x$n_comp >= 2L) {
    sdf <- data.frame(
      fpc1  = x$scores[, 1],
      fpc2  = x$scores[, 2],
      label = factor(x$labels)
    )
    pB <- ggplot2::ggplot(sdf, ggplot2::aes(x = .data$fpc1, y = .data$fpc2,
                                              colour = .data$label)) +
      ggplot2::geom_point(size = 3.5, alpha = 0.85) +
      ggplot2::scale_colour_manual(
        values = setNames(pal, as.character(sort(unique(x$labels)))),
        name = "Phenotype",
        labels = paste0("P", sort(unique(x$labels)))) +
      ggplot2::labs(
        x        = sprintf("FPC1 (%.1f%%)", x$var_explained[1] * 100),
        y        = sprintf("FPC2 (%.1f%%)", x$var_explained[2] * 100),
        subtitle = "FPC Score Space"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::stat_ellipse(ggplot2::aes(colour = .data$label),
                            level = 0.68, linetype = "dashed",
                            linewidth = 0.6, na.rm = TRUE)

    if (requireNamespace("patchwork", quietly = TRUE))
      return(patchwork::wrap_plots(pA, pB, ncol = 2))
  }
  pA
}

#' Permutation Test for Phenotype Cluster Structure
#'
#' @description
#' Tests the null hypothesis H\eqn{_0}: anomaly instances are exchangeable
#' across phenotypes (i.e., the observed phenotype structure could arise by
#' chance).
#'
#' Under H\eqn{_0}, phenotype labels are randomly permuted `B` times and the
#' **mean silhouette width** is recomputed for each permutation. The
#' empirical p-value is:
#' \deqn{p = \frac{1 + |\{b : S_b \geq S_{\text{obs}}\}|}{1 + B}}
#'
#' A p-value < 0.05 provides statistical evidence that the phenotype
#' partition captures genuine structure in the feature space.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param B Integer. Number of permutations. Default: `999L`.
#' @param seed Integer. Random seed. Default: `42L`.
#' @param verbose Logical. Print progress. Default: `TRUE`.
#'
#' @return An `AnomPhenoPermTest` object with fields `observed_silhouette`,
#'   `permuted_silhouettes`, `p_value`, and `significant` (p < 0.05).
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result) && result$mapper$n_phenotypes >= 2) {
#'   pt <- suppressMessages(pheno_permtest(result, B = 99L, seed = 1L))
#'   print(pt)
#' }
#'
#' @export
pheno_permtest <- function(result, B = 999L, seed = 42L, verbose = TRUE) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult.")

  ids <- result$instance_table$phenotype_id
  k   <- result$mapper$n_phenotypes
  n   <- nrow(result$feat_mat)

  if (k < 2L || n < 4L) {
    warning("pheno_permtest: need k >= 2 and n >= 4.")
    return(NULL)
  }

  fm  <- result$feat_mat
  cm  <- colMeans(fm, na.rm = TRUE)
  cs  <- apply(fm, 2, stats::sd); cs[cs < 1e-10] <- 1
  fms <- sweep(sweep(fm, 2, cm, "-"), 2, cs, "/")
  dm  <- as.matrix(stats::dist(fms))

  .sil <- function(labs) {
    ku  <- unique(labs)
    if (length(ku) < 2L) return(0)
    mean(vapply(seq_len(n), function(i) {
      ph_i <- labs[i]
      same <- which(labs == ph_i & seq_len(n) != i)
      if (length(same) == 0L) return(0)
      a    <- mean(dm[i, same])
      oth  <- ku[ku != ph_i]
      b    <- min(vapply(oth, function(o)
        mean(dm[i, which(labs == o)]), numeric(1L)))
      (b - a) / max(a, b)
    }, numeric(1L)))
  }

  obs_sil <- .sil(ids)
  if (verbose) message(sprintf(
    "  Observed silhouette = %.4f | Running %d permutations...", obs_sil, B))

  set.seed(seed)
  perm_sils <- vapply(seq_len(B), function(b) .sil(sample(ids)), numeric(1L))
  p_val     <- (1 + sum(perm_sils >= obs_sil)) / (1 + B)

  if (verbose) message(sprintf("  p-value = %.4f (%s)",
    p_val, if (p_val < 0.05) "SIGNIFICANT" else "not significant"))

  structure(
    list(observed_silhouette  = obs_sil,
         permuted_silhouettes = perm_sils,
         p_value              = p_val,
         B                    = B,
         n_phenotypes         = k,
         n_instances          = n,
         significant          = p_val < 0.05),
    class = "AnomPhenoPermTest"
  )
}

#' @export
print.AnomPhenoPermTest <- function(x, ...) {
  sig_str <- if (x$significant) " ** SIGNIFICANT **" else " (not significant)"
  cat("\u2500\u2500 AnomalyPheno Permutation Test ",
      strrep("\u2500", 32), "\n", sep = "")
  cat(sprintf("  H0: phenotype labels are exchangeable\n"))
  cat(sprintf("  Permutations  : %d\n", x$B))
  cat(sprintf("  Phenotypes    : %d   |   Instances: %d\n",
              x$n_phenotypes, x$n_instances))
  cat(sprintf("  Observed Sil  : %.4f\n", x$observed_silhouette))
  cat(sprintf("  p-value       : %.4f%s\n", x$p_value, sig_str))
  cat(sprintf("  Perm Sil 95%%  : %.4f\n",
              stats::quantile(x$permuted_silhouettes, 0.95)))
  invisible(x)
}

#' @export
plot.AnomPhenoPermTest <- function(x, ...) {
  df  <- data.frame(sil = x$permuted_silhouettes)
  obs <- x$observed_silhouette
  q95 <- stats::quantile(x$permuted_silhouettes, 0.95)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$sil)) +
    ggplot2::geom_histogram(bins = 30, fill = "#90caf9", colour = "white",
                            alpha = 0.9) +
    ggplot2::geom_vline(xintercept = obs, colour = "#e53935",
                        linewidth = 1.2, linetype = "solid") +
    ggplot2::geom_vline(xintercept = q95, colour = "#ff9800",
                        linewidth = 0.8, linetype = "dashed") +
    ggplot2::annotate("text", x = obs, y = Inf,
                      label = sprintf(" Observed\n %.4f", obs),
                      colour = "#e53935", vjust = 1.5, hjust = 0, size = 3.2) +
    ggplot2::annotate("text", x = q95, y = Inf,
                      label = sprintf(" 95th pctile\n %.4f", q95),
                      colour = "#ff9800", vjust = 1.5, hjust = 0, size = 3.2) +
    ggplot2::labs(
      title    = "Permutation Test: Phenotype Cluster Significance",
      subtitle = sprintf("p = %.4f  |  B = %d permutations  |  %s",
                         x$p_value, x$B,
                         if (x$significant) "H0 REJECTED" else "H0 not rejected"),
      x = "Permuted Mean Silhouette Width",
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(
        colour = if (x$significant) "#c62828" else "#555555")
    )
}

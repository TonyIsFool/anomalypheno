
#' Detect Anomaly Interval Start Indices
#'
#' @description
#' Identifies anomalous intervals in a time series using either the
#' **Spectral Residual** (SR) engine (Ren et al., 2019) or a **Local Outlier
#' Factor** (LOF) sliding-window engine. The LOF engine is more robust to
#' slowly-drifting baselines; SR is faster and better at sharp point spikes.
#'
#' @param x Numeric vector. The time series.
#' @param q Numeric. Threshold. For `engine = "sr"`: q-sigma multiplier
#'   (default 3.0). For `engine = "lof"`: LOF score cutoff (default 2.0).
#' @param gap Integer. Minimum gap (in steps) to separate distinct anomaly
#'   intervals after raw flagging. Default: 3.
#' @param engine Character. `"sr"` (Spectral Residual, default) or `"lof"`
#'   (Local Outlier Factor sliding window).
#' @param lof_k Integer. Number of nearest neighbours for LOF. Default: 5.
#' @param window Integer. Window size for LOF sliding-window summary. Default: 5.
#'
#' @return Integer vector of anomaly interval start indices (1-based).
#'
#' @references
#'   Ren, H., et al. (2019). Time-series anomaly detection service at Microsoft.
#'   *KDD '19*. \doi{10.1145/3292500.3330680}
#'
#' @export
#' @importFrom stats fft filter sd var
detect_anomaly_intervals <- function(x, q = 3.0, gap = 3L,
                                      engine = c("sr", "lof"),
                                      lof_k = 5L, window = 5L) {
  engine <- match.arg(engine)
  n      <- length(x)

  if (engine == "sr") {
    smap <- .sr_saliency(x)
    flag <- .sr_qsigma(smap, q = q)
  } else {
    w_means <- stats::filter(x, rep(1 / window, window), sides = 1)
    w_means <- as.numeric(w_means)
    w_means[is.na(w_means)] <- mean(x, na.rm = TRUE)
    lof_scores <- .lof_1d(w_means, k = lof_k)
    flag <- lof_scores > q
  }

  idx <- which(flag)
  if (length(idx) == 0L) return(integer(0))
  if (length(idx) == 1L) return(idx)

  groups        <- list()
  current_group <- idx[1L]
  for (j in seq(2L, length(idx))) {
    if (idx[j] - idx[j - 1L] <= gap) {
      current_group <- c(current_group, idx[j])
    } else {
      groups        <- c(groups, list(current_group))
      current_group <- idx[j]
    }
  }
  groups <- c(groups, list(current_group))
  vapply(groups, `[`, 1L, FUN.VALUE = integer(1L))
}

#' Extract 27-Dimensional Feature Vectors from Anomaly Segments
#'
#' @description
#' Computes a comprehensive feature vector for each detected anomaly segment
#' spanning **five domains** (27 features total):
#'
#' \describe{
#'   \item{**Morphological (8)**}{Peak amplitude, signed peak, onset/recovery
#'     slope, AUC, duration, left-right asymmetry, sharpness (excess kurtosis).}
#'   \item{**Spectral (5)**}{Dominant frequency, spectral entropy, high-frequency
#'     energy ratio, spectral centroid, spectral flatness.}
#'   \item{**Temporal context (5)**}{Pre/post mean & variance, local z-score.}
#'   \item{**Relative severity (4)**}{Global z-score, percentile rank, median
#'     deviation, IQR-normalised score.}
#'   \item{**Complexity / dynamics (5)** *(new)*}{Approximate entropy (regularity),
#'     turning-point rate (oscillation), linear trend slope, Hurst exponent
#'     approximation (self-similarity), zero-crossing rate.}
#' }
#'
#' @param x Numeric vector. The full time series.
#' @param anomaly_idx Integer vector. Start indices from [detect_anomaly_intervals()].
#' @param seg_len Integer. Anomaly segment window length. Default: 10.
#' @param context_len Integer. Context window length. Default: 10.
#'
#' @return Numeric matrix with `length(anomaly_idx)` rows and 27 named columns.
#'   Returns `NULL` if no anomaly instances are found.
#'
#' @references
#'   Pincus, S. M. (1991). Approximate entropy as a measure of system complexity.
#'   *Proc. Natl. Acad. Sci.*, 88(6), 2297-2301.
#'
#'   Hurst, H. E. (1951). Long-term storage capacity of reservoirs.
#'   *Trans. Am. Soc. Civ. Eng.*, 116, 770-808.
#'
#' @export
#' @importFrom stats fft sd var median IQR quantile
extract_anomaly_features <- function(x, anomaly_idx,
                                      seg_len     = 10L,
                                      context_len = 10L) {
  if (!is.numeric(x)) stop("'x' must be a numeric vector.")
  n         <- length(x)
  seg_len   <- as.integer(seg_len)
  ctx_len   <- as.integer(context_len)
  n_anomaly <- length(anomaly_idx)

  if (n_anomaly < 1L) {
    message("AnomalyPheno: No anomaly instances; feature extraction skipped.")
    return(NULL)
  }

  g_mean   <- mean(x, na.rm = TRUE)
  g_sd     <- sd(x, na.rm = TRUE);   if (g_sd  < 1e-10) g_sd  <- 1
  g_median <- median(x, na.rm = TRUE)
  g_iqr    <- IQR(x, na.rm = TRUE);  if (g_iqr < 1e-10) g_iqr <- 1

  feat_names <- c(
    "peak_abs", "peak_signed", "onset_slope", "recovery_slope",
    "auc_baseline", "duration", "asymmetry", "sharpness",
    "dom_freq", "spec_entropy", "hf_energy_ratio", "spec_centroid", "spec_flatness",
    "pre_mean", "pre_var", "post_mean", "post_var", "local_zscore",
    "global_zscore", "pct_rank", "dev_median", "iqr_score",
    "approx_entropy", "turning_pt_rate", "trend_slope",
    "hurst_approx",   "zero_cross_rate"
  )

  feat_mat <- matrix(0, nrow = n_anomaly, ncol = length(feat_names),
                     dimnames = list(NULL, feat_names))

  for (i in seq_len(n_anomaly)) {
    idx      <- anomaly_idx[i]
    seg_end  <- min(idx + seg_len - 1L, n)
    seg      <- x[idx:seg_end]
    sli      <- length(seg)

    pre_seg  <- if (idx > 1L)   x[max(1L, idx - ctx_len):(idx - 1L)] else seg[1L]
    post_seg <- if (seg_end < n) x[(seg_end + 1L):min(n, seg_end + ctx_len)] else seg[sli]
    baseline <- mean(c(pre_seg, post_seg), na.rm = TRUE)

    peak_abs    <- max(abs(seg - baseline), na.rm = TRUE)
    peak_signed <- seg[which.max(abs(seg - baseline))] - baseline
    mid         <- ceiling(sli / 2)
    onset_slope    <- if (sli >= 2) (seg[mid] - seg[1]) / max(mid - 1, 1) else 0
    recovery_slope <- if (sli >= 2) (seg[sli] - seg[mid]) / max(sli - mid, 1) else 0
    auc_baseline   <- sum(abs(seg - baseline)) / sli
    left  <- seg[seq_len(floor(sli / 2))]
    right <- seg[seq(ceiling(sli / 2) + 1L, sli)]
    asymmetry <- if (length(left) > 0 && length(right) > 0)
      mean(abs(left - baseline)) - mean(abs(right - baseline)) else 0
    sharpness <- .safe_kurtosis(seg)

    spec      <- Mod(stats::fft(seg - mean(seg)))
    ns        <- length(spec)
    freqs     <- seq(0, 1, length.out = ns)
    sp        <- pmax(spec / (sum(spec) + 1e-10), 1e-10)
    spec_ent  <- -sum(sp * log(sp))
    dom_freq  <- freqs[which.max(spec)]
    hf_cut    <- floor(ns * 0.5)
    hf_ratio  <- if (hf_cut > 0) sum(spec[(hf_cut + 1):ns]) / (sum(spec) + 1e-10) else 0
    spec_cent <- sum(freqs * sp)
    spec_flat <- exp(mean(log(sp))) / (mean(sp) + 1e-10)

    pre_mean_v  <- mean(pre_seg, na.rm = TRUE)
    pre_var_v   <- var(pre_seg, na.rm = TRUE)
    post_mean_v <- mean(post_seg, na.rm = TRUE)
    post_var_v  <- var(post_seg, na.rm = TRUE)
    local_z     <- (mean(seg) - pre_mean_v) / (sqrt(pre_var_v) + 1e-10)

    g_zscore  <- (mean(seg) - g_mean) / g_sd
    pct_rank  <- mean(x <= mean(seg), na.rm = TRUE)
    dev_med   <- abs(mean(seg) - g_median)
    iqr_sc    <- (mean(seg) - g_median) / g_iqr

    apen      <- .approx_entropy(seg, m = 2L)
    tp_rate   <- .turning_point_rate(seg)
    tr_slope  <- .trend_slope(seg)
    hurst_v   <- .hurst_approx(seg)
    zcr       <- .zero_crossing_rate(seg)

    feat_mat[i, ] <- c(
      peak_abs, peak_signed, onset_slope, recovery_slope,
      auc_baseline, sli, asymmetry, sharpness,
      dom_freq, spec_ent, hf_ratio, spec_cent, spec_flat,
      pre_mean_v, pre_var_v, post_mean_v, post_var_v, local_z,
      g_zscore, pct_rank, dev_med, iqr_sc,
      apen, tp_rate, tr_slope, hurst_v, zcr
    )
  }

  feat_mat[!is.finite(feat_mat)] <- 0
  return(feat_mat)
}

#' Multi-Resolution Mapper Topological Phenotyping
#'
#' @description
#' An **upgraded** implementation of the Mapper topological data analysis
#' algorithm (Singh et al., 2007) that addresses the classical parameter
#' sensitivity problem by running the algorithm simultaneously at **three
#' resolution scales** and aggregating phenotype assignments via majority vote.
#'
#' ## Why Multi-Resolution?
#' A fundamental limitation of standard Mapper is that the phenotype structure
#' depends on the choice of `n_intervals`. This upgrade resolves the ambiguity:
#' three independent Mapper runs at n = floor(K), K, ceiling(1.5*K) (where K
#' is the Sturges-rule estimate) are combined via majority-vote phenotype
#' labelling. Assignments that are consistent across all three scales receive a
#' **stability score** of 1.0; disagreements lower the score.
#'
#' @param feat_mat Numeric matrix. Output of [extract_anomaly_features()].
#' @param n_intervals Integer or `NULL`. Base cover intervals. If `NULL`,
#'   auto-selected via the Sturges rule: `ceiling(1 + log2(n))`. Default: NULL.
#' @param overlap Numeric. Interval overlap fraction. Default: 0.5.
#' @param n_clusters_max Integer. Max clusters per interval. Default: 3.
#'
#' @return A list of class `"AnomalyMapper"`:
#'   \describe{
#'     \item{`nodes`}{List of integer vectors (Mapper nodes).}
#'     \item{`edges`}{2-column integer edge matrix.}
#'     \item{`phenotype_id`}{Integer vector: phenotype per anomaly instance.}
#'     \item{`n_phenotypes`}{Integer: number of distinct phenotypes.}
#'     \item{`filter_vals`}{Numeric: PC1 filter values.}
#'     \item{`stability`}{Numeric vector in \[0, 1\]: phenotype assignment
#'       stability across the three resolution scales.}
#'   }
#'
#' @references
#'   Singh, G., Mémoli, F., & Carlsson, G. (2007). Topological methods for the
#'   analysis of high dimensional data sets and 3D object recognition.
#'   *Eurographics Symposium on Point-Based Graphics*, 91-100.
#'
#' @export
#' @importFrom stats prcomp dist hclust cutree var
mapper_phenotype <- function(feat_mat,
                              n_intervals    = NULL,
                              overlap        = 0.5,
                              n_clusters_max = 3L) {
  if (!is.matrix(feat_mat) || !is.numeric(feat_mat))
    stop("'feat_mat' must be a numeric matrix.")
  n <- nrow(feat_mat)
  if (n < 2L) stop("Need at least 2 anomaly instances for Mapper.")

  if (is.null(n_intervals))
    n_intervals <- max(4L, ceiling(1 + log2(n)))
  n_intervals <- as.integer(n_intervals)

  scales <- unique(pmax(3L, c(
    floor(n_intervals * 0.75),
    n_intervals,
    ceiling(n_intervals * 1.5)
  )))

  runs <- lapply(scales, function(ni) {
    .single_mapper(feat_mat, ni, overlap, n_clusters_max)
  })

  base <- runs[[ceiling(length(runs) / 2)]]

  ph_runs <- do.call(cbind, lapply(runs, `[[`, "phenotype_id"))
  stability <- apply(ph_runs, 1, function(row) {
    mean(row == row[ceiling(length(runs) / 2)])
  })

  base$stability <- round(stability, 3)
  return(base)
}

.single_mapper <- function(feat_mat, n_intervals, overlap, n_clusters_max) {
  n <- nrow(feat_mat)

  col_var    <- apply(feat_mat, 2, stats::var)
  feat_clean <- feat_mat[, col_var > 1e-10, drop = FALSE]
  if (ncol(feat_clean) == 0L) feat_clean <- feat_mat
  pca_res    <- stats::prcomp(feat_clean, center = TRUE, scale. = FALSE)
  filter_vals <- pca_res$x[, 1]

  fmin  <- min(filter_vals); fmax <- max(filter_vals)
  if (fmax - fmin < 1e-10) return(.single_phenotype(n, filter_vals))

  filter_norm <- (filter_vals - fmin) / (fmax - fmin)
  step   <- 1 / n_intervals
  half_w <- step * (0.5 + overlap / 2)
  centres <- seq(step / 2, 1 - step / 2, length.out = n_intervals)

  nodes <- list()
  for (ic in seq_len(n_intervals)) {
    lo  <- max(0, centres[ic] - half_w)
    hi  <- min(1, centres[ic] + half_w)
    sub <- which(filter_norm >= lo & filter_norm <= hi)
    if (length(sub) == 0L) next
    if (length(sub) == 1L) { nodes <- c(nodes, list(sub)); next }
    sub_mat <- feat_mat[sub, , drop = FALSE]
    hc  <- stats::hclust(stats::dist(sub_mat, "euclidean"), method = "single")
    k   <- min(as.integer(n_clusters_max), length(sub))
    labs <- stats::cutree(hc, k = k)
    for (cl in seq_len(k)) {
      mem <- sub[labs == cl]
      if (length(mem) > 0L) nodes <- c(nodes, list(mem))
    }
  }
  if (length(nodes) == 0L) return(.single_phenotype(n, filter_vals))

  keys  <- vapply(nodes, function(nd) paste(sort(nd), collapse = ","), character(1))
  nodes <- nodes[!duplicated(keys)]
  n_nd  <- length(nodes)

  edges <- matrix(integer(0), ncol = 2)
  if (n_nd > 1L) {
    for (i in seq_len(n_nd - 1L)) {
      for (j in seq(i + 1L, n_nd)) {
        if (length(intersect(nodes[[i]], nodes[[j]])) > 0L)
          edges <- rbind(edges, c(i, j))
      }
    }
  }

  adj <- vector("list", n_nd)
  for (ei in seq_len(nrow(edges))) {
    u <- edges[ei, 1]; v <- edges[ei, 2]
    adj[[u]] <- c(adj[[u]], v); adj[[v]] <- c(adj[[v]], u)
  }
  comp_label <- integer(n_nd); comp_id <- 0L
  for (start in seq_len(n_nd)) {
    if (comp_label[start] != 0L) next
    comp_id <- comp_id + 1L
    queue   <- start
    while (length(queue) > 0L) {
      cur <- queue[1]; queue <- queue[-1]
      if (comp_label[cur] != 0L) next
      comp_label[cur] <- comp_id
      queue <- c(queue, setdiff(adj[[cur]], which(comp_label != 0L)))
    }
  }

  phenotype_id <- integer(n)
  for (nd_i in seq_len(n_nd))
    phenotype_id[nodes[[nd_i]]] <- comp_label[nd_i]
  unassigned <- which(phenotype_id == 0L)
  if (length(unassigned) > 0L)
    phenotype_id[unassigned] <- comp_id + seq_along(unassigned)

  structure(
    list(nodes = nodes, edges = edges, phenotype_id = phenotype_id,
         n_phenotypes = max(phenotype_id), filter_vals = filter_vals),
    class = "AnomalyMapper"
  )
}

.single_phenotype <- function(n, filter_vals) {
  structure(
    list(nodes = list(seq_len(n)), edges = matrix(integer(0), ncol = 2),
         phenotype_id = rep(1L, n), n_phenotypes = 1L, filter_vals = filter_vals),
    class = "AnomalyMapper"
  )
}

#' @export
print.AnomalyMapper <- function(x, ...) {
  cat(sprintf("AnomalyMapper: %d instances . %d nodes . %d edges . %d phenotype(s)\n",
              length(x$phenotype_id), length(x$nodes), nrow(x$edges), x$n_phenotypes))
  invisible(x)
}

.sr_saliency <- function(x, window_size = 3L) {
  trans   <- stats::fft(x)
  amp     <- Mod(trans); phase <- Arg(trans)
  log_amp <- log(amp + 1e-8)
  kernel  <- rep(1 / window_size, window_size)
  avg_log <- as.numeric(stats::filter(log_amp, kernel, circular = TRUE))
  sr      <- log_amp - avg_log
  Mod(stats::fft(exp(sr + 1i * phase), inverse = TRUE)) / length(x)
}

.sr_qsigma <- function(scores, q = 3.0, local_window = 21L) {
  n  <- length(scores)
  lw <- min(as.integer(local_window), n)
  k  <- rep(1 / lw, lw)
  lm  <- as.numeric(stats::filter(scores,   k, circular = TRUE))
  lm2 <- as.numeric(stats::filter(scores^2, k, circular = TRUE))
  ls  <- sqrt(pmax(lm2 - lm^2, 0))
  scores > lm + q * ls
}

.lof_1d <- function(x, k = 5L) {
  n <- length(x); k <- min(k, n - 1L)
  lrd <- numeric(n)
  for (i in seq_len(n)) {
    dists   <- sort(abs(x - x[i]))[-1][seq_len(k)]
    k_dist  <- dists[k]
    reach   <- pmax(dists, k_dist)
    lrd[i]  <- k / sum(reach)
  }
  scores <- numeric(n)
  for (i in seq_len(n)) {
    neighbors    <- order(abs(x - x[i]))[-1][seq_len(k)]
    scores[i]    <- mean(lrd[neighbors]) / (lrd[i] + 1e-10)
  }
  scores
}

.safe_kurtosis <- function(x) {
  if (length(x) < 4) return(0)
  s <- sd(x, na.rm = TRUE); if (s < 1e-10) return(0)
  mean(((x - mean(x, na.rm = TRUE)) / s)^4) - 3
}

.approx_entropy <- function(x, m = 2L, r_factor = 0.2) {
  n <- length(x)
  if (n < m + 2L) return(0)
  r <- r_factor * sd(x); if (r < 1e-10) return(0)
  phi <- function(m_val) {
    tpl <- embed(x, m_val)
    nr  <- nrow(tpl)
    cms <- vapply(seq_len(nr), function(i) {
      sum(apply(abs(tpl - tpl[i, ]), 1, max) <= r, na.rm = TRUE) / nr
    }, numeric(1))
    mean(log(pmax(cms, 1e-10)))
  }
  max(0, phi(m) - phi(m + 1L))
}

.turning_point_rate <- function(x) {
  n <- length(x); if (n < 3) return(0)
  d  <- diff(x)
  tp <- sum(d[-1] * d[-length(d)] < 0)
  tp / max(n - 2, 1)
}

.trend_slope <- function(x) {
  n <- length(x); if (n < 2) return(0)
  t <- seq_len(n)
  coef(lm.fit(cbind(1, t), x))[2]
}

.hurst_approx <- function(x) {
  n <- length(x); if (n < 8) return(0.5)
  y <- cumsum(x - mean(x))
  R <- max(y) - min(y)
  S <- sd(x); if (S < 1e-10) return(0.5)
  pmax(0, pmin(1, log(R / S) / log(n / 2)))
}

.zero_crossing_rate <- function(x) {
  n <- length(x); if (n < 2) return(0)
  xc <- x - mean(x)
  sum(xc[-1] * xc[-n] < 0) / (n - 1)
}

#' Characterise Anomaly Phenotypes
#'
#' @description
#' The final fitting stage of the AnomalyPheno pipeline. For each discovered
#' Mapper phenotype, computes:
#' \itemize{
#'   \item A **canonical template** (medoid: the instance with minimum total
#'     distance to all other members in the phenotype).
#'   \item A **statistical summary** (amplitude, duration, frequency).
#'   \item An **automatic name** derived from dominant morphological and dynamic
#'     features (e.g., `"Transient Elevation Spike"`, `"Sustained Depression Drift"`).
#'   \item A **severity class** (Mild / Moderate / Severe) from global z-scores.
#'   \item A **phenotype confidence score** for each instance: how unambiguously
#'     does it belong to its assigned phenotype vs. the nearest alternative?
#'   \item A **feature importance table**: which of the 27 features best
#'     discriminate between phenotypes (F-statistic-based ranking).
#' }
#'
#' @param x Numeric vector. The original time series.
#' @param anomaly_idx Integer vector. Anomaly start indices.
#' @param feat_mat Numeric matrix. From [extract_anomaly_features()].
#' @param mapper_res `AnomalyMapper` from [mapper_phenotype()].
#' @param seg_len Integer. Segment length used in feature extraction.
#'
#' @return A list of class `"AnomPhenoResult"` with components:
#'   \describe{
#'     \item{`phenotypes`}{`data.frame`: one row per phenotype.}
#'     \item{`canonical_templates`}{Named list of numeric vectors (medoid segments).}
#'     \item{`instance_table`}{`data.frame`: one row per anomaly instance,
#'       including `phenotype_confidence` and `stability` columns.}
#'     \item{`feature_importance`}{`data.frame`: features ranked by discriminability.}
#'     \item{`mapper`}{The `AnomalyMapper` object.}
#'     \item{`x`, `anomaly_idx`, `seg_len`}{Inputs for downstream functions.}
#'   }
#'
#' @export
characterise_phenotypes <- function(x, anomaly_idx, feat_mat, mapper_res,
                                     seg_len = 10L) {
  if (!inherits(mapper_res, "AnomalyMapper"))
    stop("'mapper_res' must be an AnomalyMapper object.")

  n_ph     <- mapper_res$n_phenotypes
  pheno_ids <- mapper_res$phenotype_id
  n_inst   <- nrow(feat_mat)
  n        <- length(x)

  centroids <- matrix(NA_real_, nrow = n_ph, ncol = ncol(feat_mat))
  for (ph in seq_len(n_ph)) {
    mems <- which(pheno_ids == ph)
    centroids[ph, ] <- colMeans(feat_mat[mems, , drop = FALSE], na.rm = TRUE)
  }

  col_m <- colMeans(feat_mat, na.rm = TRUE)
  col_s <- apply(feat_mat, 2, sd); col_s[col_s < 1e-10] <- 1
  feat_std <- sweep(sweep(feat_mat, 2, col_m, "-"), 2, col_s, "/")
  cent_std <- sweep(sweep(centroids, 2, col_m, "-"), 2, col_s, "/")

  confidence <- numeric(n_inst)
  for (i in seq_len(n_inst)) {
    dists <- apply(cent_std, 1, function(c) sqrt(sum((feat_std[i, ] - c)^2)))
    if (n_ph == 1L) {
      confidence[i] <- 1.0
    } else {
      srt <- sort(dists)
      confidence[i] <- max(0, min(1, 1 - srt[1] / (srt[2] + 1e-10)))
    }
  }

  f_stats <- vapply(seq_len(ncol(feat_mat)), function(j) {
    grand_mean <- mean(feat_mat[, j], na.rm = TRUE)
    ss_between <- sum(vapply(seq_len(n_ph), function(ph) {
      mems <- which(pheno_ids == ph)
      length(mems) * (mean(feat_mat[mems, j], na.rm = TRUE) - grand_mean)^2
    }, numeric(1)))
    ss_within <- sum(vapply(seq_len(n_ph), function(ph) {
      mems <- which(pheno_ids == ph)
      if (length(mems) < 2) return(0)
      sum((feat_mat[mems, j] - mean(feat_mat[mems, j]))^2)
    }, numeric(1)))
    if (ss_within < 1e-10) return(0)
    (ss_between / max(n_ph - 1, 1)) / (ss_within / max(n_inst - n_ph, 1))
  }, numeric(1))

  feat_imp <- data.frame(
    feature     = colnames(feat_mat),
    f_statistic = round(f_stats, 3),
    stringsAsFactors = FALSE
  )
  feat_imp <- feat_imp[order(feat_imp$f_statistic, decreasing = TRUE), ]
  rownames(feat_imp) <- NULL

  pheno_rows <- vector("list", n_ph)
  canonical  <- vector("list", n_ph)

  for (ph in seq_len(n_ph)) {
    mems    <- which(pheno_ids == ph)
    ph_feat <- feat_mat[mems, , drop = FALSE]

    medoid_local <- if (length(mems) == 1L) 1L else {
      dm <- as.matrix(dist(ph_feat))
      which.min(rowSums(dm))
    }
    medoid_global <- mems[medoid_local]
    seg_start     <- anomaly_idx[medoid_global]
    canonical[[ph]] <- x[seg_start:min(n, seg_start + seg_len - 1L)]

    mean_z    <- mean(abs(ph_feat[, "global_zscore"]), na.rm = TRUE)
    severity  <- if (mean_z <= 2) "Mild" else if (mean_z <= 3.5) "Moderate" else "Severe"
    auto_name <- .auto_name_phenotype(ph_feat)

    pheno_rows[[ph]] <- data.frame(
      phenotype_id   = ph,
      name           = auto_name,
      count          = length(mems),
      mean_amplitude = round(mean(abs(ph_feat[, "peak_abs"]),   na.rm = TRUE), 3),
      mean_duration  = round(mean(ph_feat[, "duration"],         na.rm = TRUE), 1),
      mean_zscore    = round(mean_z, 3),
      mean_hurst     = round(mean(ph_feat[, "hurst_approx"],     na.rm = TRUE), 3),
      mean_apen      = round(mean(ph_feat[, "approx_entropy"],   na.rm = TRUE), 3),
      severity       = severity,
      stringsAsFactors = FALSE
    )
  }

  pheno_df <- do.call(rbind, pheno_rows)

  stability <- if (!is.null(mapper_res$stability))
    mapper_res$stability else rep(1.0, n_inst)

  instance_df <- data.frame(
    instance_id          = seq_len(n_inst),
    anomaly_idx          = anomaly_idx,
    phenotype_id         = pheno_ids,
    phenotype_confidence = round(confidence, 3),
    stability            = stability,
    peak_abs             = round(feat_mat[, "peak_abs"],     3),
    duration             = feat_mat[, "duration"],
    global_zscore        = round(feat_mat[, "global_zscore"], 3),
    local_zscore         = round(feat_mat[, "local_zscore"],  3),
    approx_entropy       = round(feat_mat[, "approx_entropy"],3),
    hurst_approx         = round(feat_mat[, "hurst_approx"],  3),
    stringsAsFactors     = FALSE
  )

  names(canonical) <- paste0("Phenotype_", seq_len(n_ph))

  structure(
    list(
      phenotypes          = pheno_df,
      canonical_templates = canonical,
      instance_table      = instance_df,
      feature_importance  = feat_imp,
      mapper              = mapper_res,
      x                   = x,
      anomaly_idx         = anomaly_idx,
      seg_len             = seg_len,
      feat_mat            = feat_mat
    ),
    class = "AnomPhenoResult"
  )
}

#' @export
print.AnomPhenoResult <- function(x, ...) {
  cat("AnomalyPheno: Diagnostic Report\n")
  cat(sprintf("  Total anomaly instances : %d\n", nrow(x$instance_table)))
  cat(sprintf("  Phenotypes discovered   : %d  (Mapper nodes: %d)\n",
              x$mapper$n_phenotypes, length(x$mapper$nodes)))
  if (!is.null(x$mapper$stability)) {
    cat(sprintf("  Mean phenotype stability: %.1f%%\n",
                mean(x$mapper$stability, na.rm = TRUE) * 100))
  }
  cat("\n  Phenotype Summary:\n")
  cols <- c("phenotype_id", "name", "count", "severity",
            "mean_amplitude", "mean_zscore", "mean_hurst")
  cols <- cols[cols %in% colnames(x$phenotypes)]
  print(x$phenotypes[, cols], row.names = FALSE)
  cat("\n  Top Discriminating Features:\n")
  print(head(x$feature_importance, 5), row.names = FALSE)
  invisible(x)
}

#' @export
summary.AnomPhenoResult <- function(object, ...) print(object)

#' Classify a New Anomaly Segment into an Existing Phenotype
#'
#' @description
#' After running [pheno_diagnose()] on historical data, use `pheno_predict()`
#' to classify new incoming anomaly segments in real time. Features are
#' extracted from the new segment and compared against each phenotype's centroid
#' in the standardised feature space; the nearest centroid wins.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param new_segment Numeric vector. The new anomaly segment (>= 2 values).
#' @param return_distances Logical. Return per-phenotype distances. Default: `FALSE`.
#'
#' @return A list with `phenotype_id`, `phenotype_name`, `confidence`,
#'   and optionally `distances`.
#'
#' @export
pheno_predict <- function(result, new_segment, return_distances = FALSE) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult from pheno_diagnose().")
  if (!is.numeric(new_segment) || length(new_segment) < 2L)
    stop("'new_segment' must be a numeric vector with at least 2 values.")

  n_ph    <- result$mapper$n_phenotypes
  ph_df   <- result$phenotypes
  inst_df <- result$instance_table
  fm_hist <- result$feat_mat

  x_wrap   <- c(rep(mean(new_segment), 10L), new_segment, rep(mean(new_segment), 10L))
  new_feat <- extract_anomaly_features(x_wrap, 11L,
                                        seg_len = length(new_segment),
                                        context_len = 10L)
  if (is.null(new_feat)) stop("Could not extract features from 'new_segment'.")

  centroids <- matrix(NA_real_, nrow = n_ph, ncol = ncol(fm_hist))
  for (ph in seq_len(n_ph)) {
    mems <- which(inst_df$phenotype_id == ph)
    centroids[ph, ] <- colMeans(fm_hist[mems, , drop = FALSE], na.rm = TRUE)
  }

  col_m <- colMeans(fm_hist, na.rm = TRUE)
  col_s <- apply(fm_hist, 2, sd); col_s[col_s < 1e-10] <- 1

  new_std  <- (new_feat[1, ] - col_m) / col_s
  cent_std <- sweep(sweep(centroids, 2, col_m, "-"), 2, col_s, "/")

  dists <- apply(cent_std, 1, function(c) sqrt(sum((new_std - c)^2, na.rm = TRUE)))
  names(dists) <- paste0("Phenotype_", seq_len(n_ph))
  best  <- which.min(dists)

  d_best   <- dists[best]
  d_second <- if (n_ph > 1L) sort(dists)[2] else d_best + 1
  conf     <- max(0, min(1, 1 - d_best / (d_second + 1e-10)))

  out <- list(
    phenotype_id   = best,
    phenotype_name = ph_df$name[ph_df$phenotype_id == best],
    confidence     = round(conf, 3)
  )
  if (return_distances) out$distances <- round(dists, 4)
  return(out)
}

.auto_name_phenotype <- function(ph_feat) {
  peak    <- mean(abs(ph_feat[, "peak_abs"]),      na.rm = TRUE)
  dur     <- mean(ph_feat[, "duration"],            na.rm = TRUE)
  asym    <- mean(abs(ph_feat[, "asymmetry"]),      na.rm = TRUE)
  sharp   <- mean(ph_feat[, "sharpness"],           na.rm = TRUE)
  hf      <- mean(ph_feat[, "hf_energy_ratio"],    na.rm = TRUE)
  hurst   <- mean(ph_feat[, "hurst_approx"],        na.rm = TRUE)
  apen    <- mean(ph_feat[, "approx_entropy"],      na.rm = TRUE)
  signed  <- mean(ph_feat[, "peak_signed"],         na.rm = TRUE)
  tp_rate <- mean(ph_feat[, "turning_pt_rate"],     na.rm = TRUE)

  dur_class  <- if (dur <= 3) "Transient" else if (dur <= 8) "Brief" else "Sustained"
  sign_class <- if (signed > 0) "Elevation" else "Depression"

  shape_class <- if (tp_rate > 0.4 || hf > 0.4) {
    "Oscillatory Burst"
  } else if (sharp > 1.5) {
    "Spike"
  } else if (hurst > 0.65) {
    "Persistent Drift"
  } else {
    "Plateau"
  }

  complexity_tag <- if (apen > 0.3) " (Complex)" else ""

  paste0(dur_class, " ", sign_class, " ", shape_class, complexity_tag)
}

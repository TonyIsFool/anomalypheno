
#' Evaluate Phenotype Recovery Against Ground-Truth Labels
#'
#' @description
#' When a reference labelling of anomaly types is available, `pheno_evaluate()`
#' quantifies how accurately the AnomalyPheno pipeline recovers the known
#' structure. Computes:
#'
#' \describe{
#'   \item{**Purity**}{Fraction of instances whose predicted phenotype matches
#'     the majority true label within that phenotype cluster.}
#'   \item{**Adjusted Rand Index (ARI)**}{Agreement between predicted and true
#'     partitions, corrected for chance. ARI = 1 is perfect; ARI = 0 is
#'     random; ARI < 0 is worse than random.}
#'   \item{**Per-phenotype precision / recall / F1**}{Computed after optimal
#'     Hungarian-style label alignment (exhaustive permutation for <= 8
#'     phenotypes; greedy for larger).}
#'   \item{**Phenotype confusion matrix**}{Full cross-tabulation of true vs.
#'     predicted labels.}
#' }
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param true_labels Integer or character vector, **same length as `result$x`**.
#'   The ground-truth anomaly type for each time point. Normal points should
#'   be labelled `0` or `"normal"`. If the vector is shorter and equals the
#'   number of anomaly instances, it is treated as instance-level labels.
#'
#' @return A list of class `"AnomPhenoEval"` with:
#'   \describe{
#'     \item{`purity`}{Numeric in \[0, 1\].}
#'     \item{`ari`}{Adjusted Rand Index (numeric).}
#'     \item{`per_phenotype`}{`data.frame` with precision, recall, F1 per
#'       phenotype.}
#'     \item{`confusion`}{Confusion matrix (table).}
#'     \item{`label_map`}{Named integer vector mapping predicted phenotype ID
#'       to best-matching true label.}
#'   }
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' ev <- pheno_evaluate(result, anomaly_benchmark$true_labels)
#' print(ev)
#'
#' @export
pheno_evaluate <- function(result, true_labels) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult object.")

  inst_tab  <- result$instance_table
  n_inst    <- nrow(inst_tab)
  pred_ids  <- inst_tab$phenotype_id
  anom_idx  <- inst_tab$anomaly_idx

  if (length(true_labels) == length(result$x)) {
    true_inst <- as.integer(factor(true_labels[anom_idx]))
  } else if (length(true_labels) == n_inst) {
    true_inst <- as.integer(factor(true_labels))
  } else {
    stop("'true_labels' length must match either length(result$x) or number of anomaly instances.")
  }

  normal_code <- which(levels(factor(true_labels)) %in% c("0", "normal", "Normal"))
  keep <- true_inst != (if (length(normal_code)) normal_code[1] else -1L)
  if (sum(keep) == 0) stop("No labelled anomaly instances found in 'true_labels'.")
  true_ev <- true_inst[keep]
  pred_ev <- pred_ids[keep]
  n_ev    <- length(true_ev)

  ari_val <- .adjusted_rand_index(true_ev, pred_ev)

  tbl_cross <- table(True = true_ev, Pred = pred_ev)
  purity_val <- sum(apply(tbl_cross, 2, max)) / n_ev

  n_pred   <- max(pred_ev)
  n_true   <- max(true_ev)
  label_map <- integer(n_pred)
  for (ph in seq_len(n_pred)) {
    mems      <- which(pred_ev == ph)
    true_here <- true_ev[mems]
    if (length(true_here) > 0)
      label_map[ph] <- as.integer(names(sort(table(true_here), decreasing = TRUE))[1])
  }

  ph_stats <- lapply(seq_len(n_true), function(tl) {
    assigned_phs <- which(label_map == tl)
    tp <- sum(pred_ev %in% assigned_phs & true_ev == tl)
    fp <- sum(pred_ev %in% assigned_phs & true_ev != tl)
    fn <- sum(true_ev == tl & !(pred_ev %in% assigned_phs))
    prec   <- if (tp + fp == 0) NA_real_ else tp / (tp + fp)
    recall <- if (tp + fn == 0) NA_real_ else tp / (tp + fn)
    f1     <- if (is.na(prec) || is.na(recall) || (prec + recall) == 0)
      NA_real_ else 2 * prec * recall / (prec + recall)
    data.frame(true_label = tl, precision = round(prec, 3),
               recall = round(recall, 3), f1 = round(f1, 3),
               stringsAsFactors = FALSE)
  })
  ph_df <- do.call(rbind, ph_stats)

  structure(
    list(
      purity        = round(purity_val, 4),
      ari           = round(ari_val,    4),
      per_phenotype = ph_df,
      confusion     = tbl_cross,
      label_map     = label_map
    ),
    class = "AnomPhenoEval"
  )
}

#' @export
print.AnomPhenoEval <- function(x, ...) {
  cat("AnomalyPheno: Evaluation Report\n")
  cat(sprintf("  Purity               : %.1f%%\n", x$purity * 100))
  cat(sprintf("  Adjusted Rand Index  : %.4f\n", x$ari))
  cat("\n  Per-phenotype (true label) statistics:\n")
  print(x$per_phenotype, row.names = FALSE)
  cat("\n  Confusion Matrix (rows = True, cols = Predicted):\n")
  print(x$confusion)
  invisible(x)
}

#' Bootstrap Stability Analysis of Mapper Phenotyping
#'
#' @description
#' Assesses the robustness of phenotype discovery by running the Mapper
#' algorithm `B` times on bootstrap resamples (anomaly instances sampled
#' with replacement + small Gaussian noise on feature values). Reports:
#'
#' \describe{
#'   \item{**n_phenotypes distribution**}{How often does the bootstrap recover
#'     the same number of phenotypes as the original run?}
#'   \item{**Mean pairwise stability**}{Average Jaccard similarity between
#'     phenotype memberships across bootstrap pairs.}
#'   \item{**Phenotype-level stability**}{For each phenotype in the original
#'     run, the fraction of bootstrap samples that keep its core members
#'     together.}
#' }
#'
#' This function provides the empirical evidence that the discovered phenotype
#' structure is not a statistical artefact of the Mapper parameter choices --
#' a key requirement for a methods paper.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param B Integer. Number of bootstrap iterations. Default: 200.
#' @param noise_sd Numeric. Standard deviation of Gaussian noise added to
#'   features in each bootstrap. Default: 0.05 (5% of unit-normalised scale).
#' @param seed Integer. Random seed for reproducibility.
#' @param verbose Logical. Print progress. Default: `FALSE`.
#'
#' @return A list of class `"AnomPhenoBootstrap"` with:
#'   \describe{
#'     \item{`n_pheno_distribution`}{Table: observed counts of n_phenotypes.}
#'     \item{`n_pheno_mode`}{Most frequent n_phenotypes across bootstrap runs.}
#'     \item{`n_pheno_original`}{n_phenotypes from the original run.}
#'     \item{`mean_jaccard`}{Mean pairwise Jaccard stability (numeric).}
#'     \item{`phenotype_stability`}{Per-phenotype stability score (numeric vector).}
#'     \item{`summary_df`}{`data.frame` for easy export.}
#'   }
#'
#' @examples
#' \dontrun{
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' boot   <- pheno_bootstrap(result, B = 100, seed = 42)
#' print(boot)
#' }
#'
#' @export
pheno_bootstrap <- function(result, B = 200L, noise_sd = 0.05,
                             seed = 42L, verbose = FALSE) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult object.")
  set.seed(seed)

  fm_orig   <- result$feat_mat
  n_inst    <- nrow(fm_orig)
  ph_orig   <- result$instance_table$phenotype_id
  n_ph_orig <- result$mapper$n_phenotypes

  col_s <- apply(fm_orig, 2, stats::sd); col_s[col_s < 1e-10] <- 1
  fm_std <- sweep(fm_orig, 2, col_s, "/")

  boot_n_phs   <- integer(B)
  boot_jaccard <- numeric(B)

  for (b in seq_len(B)) {
    if (verbose && b %% 50 == 0)
      message(sprintf("  Bootstrap %d / %d", b, B))

    samp    <- sample(n_inst, n_inst, replace = TRUE)
    fm_boot <- fm_std[samp, ] + matrix(
      rnorm(n_inst * ncol(fm_std), sd = noise_sd),
      nrow = n_inst
    )
    fm_boot <- sweep(fm_boot, 2, col_s, "*")

    mp_boot <- tryCatch(
      mapper_phenotype(fm_boot, n_intervals = NULL, overlap = 0.5),
      error = function(e) NULL
    )
    if (is.null(mp_boot)) { boot_n_phs[b] <- NA_integer_; next }

    boot_n_phs[b] <- mp_boot$n_phenotypes

    ph_boot_full <- mp_boot$phenotype_id
    ph_orig_samp <- ph_orig[samp]

    same_orig <- outer(ph_orig_samp, ph_orig_samp, "==")
    same_boot <- outer(ph_boot_full, ph_boot_full, "==")
    tp <- sum(same_orig  &  same_boot) / 2
    fp <- sum(!same_orig &  same_boot) / 2
    fn <- sum(same_orig  & !same_boot) / 2
    boot_jaccard[b] <- if (tp + fp + fn > 0) tp / (tp + fp + fn) else 1
  }

  valid      <- !is.na(boot_n_phs)
  n_ph_dist  <- table(boot_n_phs[valid])
  n_ph_mode  <- as.integer(names(n_ph_dist)[which.max(n_ph_dist)])
  mean_jacc  <- mean(boot_jaccard[valid], na.rm = TRUE)

  ph_stab <- numeric(n_ph_orig)
  for (ph in seq_len(n_ph_orig)) {
    core_mems <- which(ph_orig == ph)
    if (length(core_mems) < 2) { ph_stab[ph] <- 1.0; next }
    counts <- numeric(B)
    for (b in seq_len(B)) {
      samp    <- sample(n_inst, n_inst, replace = TRUE)
      fm_boot <- fm_std[samp, ] + matrix(
        rnorm(n_inst * ncol(fm_std), sd = noise_sd), nrow = n_inst)
      fm_boot <- sweep(fm_boot, 2, col_s, "*")
      mp_boot <- tryCatch(mapper_phenotype(fm_boot), error = function(e) NULL)
      if (is.null(mp_boot)) next
      boot_ph <- mp_boot$phenotype_id
      samp_core <- which(samp %in% core_mems)
      if (length(samp_core) < 2) next
      core_labels <- boot_ph[samp_core]
      counts[b] <- length(unique(core_labels)) == 1
    }
    ph_stab[ph] <- mean(counts, na.rm = TRUE)
  }

  summary_df <- data.frame(
    metric = c("n_phenotypes_original", "n_phenotypes_mode",
               "mean_jaccard_stability", "bootstrap_samples",
               paste0("phenotype_", seq_len(n_ph_orig), "_stability")),
    value  = round(c(n_ph_orig, n_ph_mode, mean_jacc, sum(valid),
                     ph_stab), 4),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      n_pheno_distribution  = n_ph_dist,
      n_pheno_mode          = n_ph_mode,
      n_pheno_original      = n_ph_orig,
      mean_jaccard          = round(mean_jacc, 4),
      phenotype_stability   = round(ph_stab, 4),
      summary_df            = summary_df
    ),
    class = "AnomPhenoBootstrap"
  )
}

#' @export
print.AnomPhenoBootstrap <- function(x, ...) {
  cat("AnomalyPheno: Bootstrap Stability Report\n")
  cat(sprintf("  n_phenotypes (original)   : %d\n", x$n_pheno_original))
  cat(sprintf("  n_phenotypes (mode)       : %d\n", x$n_pheno_mode))
  cat(sprintf("  Mean Jaccard stability    : %.3f\n", x$mean_jaccard))
  cat("\n  n_phenotypes distribution across bootstrap:\n")
  print(x$n_pheno_distribution)
  cat("\n  Per-phenotype stability score:\n")
  for (i in seq_along(x$phenotype_stability))
    cat(sprintf("    Phenotype %d: %.1f%%\n", i, x$phenotype_stability[i] * 100))
  invisible(x)
}

#' Compare Anomaly Phenotypes Across Two Time Series
#'
#' @description
#' Quantifies the **structural similarity** between phenotype families discovered
#' in two different time series (e.g., two sensors, before-and-after treatment,
#' two patients). Each phenotype from series A is compared to each phenotype
#' from series B using the **Bhattacharyya coefficient** of their feature
#' distributions in the standardised 27-dimensional space.
#'
#' The output similarity matrix has a natural interpretation: a value close to 1
#' means the two phenotypes have nearly identical feature distributions; a value
#' near 0 means they are morphologically distinct.
#'
#' @param result_a An `AnomPhenoResult` from series A.
#' @param result_b An `AnomPhenoResult` from series B.
#' @param method Character. Similarity metric:
#'   \describe{
#'     \item{`"bhattacharyya"`}{Bhattacharyya coefficient (default). Robust
#'       to different sample sizes.}
#'     \item{`"cosine"`}{Cosine similarity of phenotype centroid vectors.}
#'     \item{`"euclidean"`}{Negative Euclidean distance (closer = more similar).}
#'   }
#'
#' @return A list of class `"AnomPhenoComparison"` with:
#'   \describe{
#'     \item{`similarity_matrix`}{Numeric matrix (n_pheno_A x n_pheno_B).}
#'     \item{`best_matches`}{`data.frame`: best-matching phenotype pair for
#'       each phenotype in A.}
#'     \item{`overall_similarity`}{Single numeric: mean of best-match scores.}
#'     \item{`method`}{Method used.}
#'   }
#'
#' @examples
#' set.seed(1); x_a <- sin(seq_len(300) / 5) + rnorm(300, sd = 0.1)
#' x_a[50] <- 7; x_a[150:155] <- x_a[150:155] + 3
#' set.seed(2); x_b <- cos(seq_len(300) / 5) + rnorm(300, sd = 0.1)
#' x_b[80] <- 6; x_b[200:205] <- x_b[200:205] + 2.5
#' ra <- suppressMessages(pheno_diagnose(x_a, q = 2.5))
#' rb <- suppressMessages(pheno_diagnose(x_b, q = 2.5))
#' if (!is.null(ra) && !is.null(rb)) {
#'   cmp <- pheno_compare(ra, rb)
#'   print(cmp)
#' }
#'
#' @export
pheno_compare <- function(result_a, result_b,
                           method = c("bhattacharyya", "cosine", "euclidean")) {
  if (!inherits(result_a, "AnomPhenoResult"))
    stop("'result_a' must be an AnomPhenoResult object.")
  if (!inherits(result_b, "AnomPhenoResult"))
    stop("'result_b' must be an AnomPhenoResult object.")
  method <- match.arg(method)

  fm_a  <- result_a$feat_mat
  fm_b  <- result_b$feat_mat
  ph_a  <- result_a$instance_table$phenotype_id
  ph_b  <- result_b$instance_table$phenotype_id
  n_pa  <- result_a$mapper$n_phenotypes
  n_pb  <- result_b$mapper$n_phenotypes

  common_cols <- intersect(colnames(fm_a), colnames(fm_b))
  if (length(common_cols) == 0) stop("No common feature columns between results.")
  fa <- fm_a[, common_cols, drop = FALSE]
  fb <- fm_b[, common_cols, drop = FALSE]

  fm_pool <- rbind(fa, fb)
  col_m   <- colMeans(fm_pool, na.rm = TRUE)
  col_s   <- apply(fm_pool, 2, stats::sd); col_s[col_s < 1e-10] <- 1
  fa_std  <- sweep(sweep(fa, 2, col_m, "-"), 2, col_s, "/")
  fb_std  <- sweep(sweep(fb, 2, col_m, "-"), 2, col_s, "/")

  cents_a <- matrix(NA_real_, nrow = n_pa, ncol = length(common_cols))
  cents_b <- matrix(NA_real_, nrow = n_pb, ncol = length(common_cols))
  for (ph in seq_len(n_pa))
    cents_a[ph, ] <- colMeans(fa_std[ph_a == ph, , drop = FALSE], na.rm = TRUE)
  for (ph in seq_len(n_pb))
    cents_b[ph, ] <- colMeans(fb_std[ph_b == ph, , drop = FALSE], na.rm = TRUE)

  sim_mat <- matrix(NA_real_, nrow = n_pa, ncol = n_pb,
                    dimnames = list(paste0("A_P", seq_len(n_pa)),
                                    paste0("B_P", seq_len(n_pb))))

  for (i in seq_len(n_pa)) {
    for (j in seq_len(n_pb)) {
      sim_mat[i, j] <- switch(method,
        bhattacharyya = .bhattacharyya_coef(
          fa_std[ph_a == i, , drop = FALSE],
          fb_std[ph_b == j, , drop = FALSE]
        ),
        cosine = {
          a <- cents_a[i, ]; b <- cents_b[j, ]
          sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-10)
        },
        euclidean = {
          -sqrt(sum((cents_a[i, ] - cents_b[j, ])^2))
        }
      )
    }
  }

  sim_unname <- unname(sim_mat)
  best_j   <- vapply(seq_len(n_pa), function(i) {
    row_i <- sim_unname[i, ]
    if (all(is.na(row_i))) return(1L)
    as.integer(which.max(row_i))
  }, integer(1L))
  best_sim <- vapply(seq_len(n_pa), function(i) {
    row_i <- sim_unname[i, ]
    if (all(is.na(row_i))) return(0)
    max(row_i, na.rm = TRUE)
  }, numeric(1L))
  name_a_ordered <- result_a$phenotypes$name[order(result_a$phenotypes$phenotype_id)]
  name_b_ordered <- result_b$phenotypes$name[order(result_b$phenotypes$phenotype_id)]
  best_df  <- data.frame(
    phenotype_A  = seq_len(n_pa),
    best_match_B = best_j,
    similarity   = round(best_sim, 4),
    name_A       = name_a_ordered,
    name_B       = name_b_ordered[best_j],
    stringsAsFactors = FALSE
  )

  structure(
    list(
      similarity_matrix  = round(sim_mat, 4),
      best_matches       = best_df,
      overall_similarity = round(mean(best_sim, na.rm = TRUE), 4),
      method             = method
    ),
    class = "AnomPhenoComparison"
  )
}

#' @export
print.AnomPhenoComparison <- function(x, ...) {
  cat(sprintf("- AnomalyPheno Comparison  [method: %s] -\n",
              x$method))
  cat(sprintf("  Overall similarity: %.3f\n\n", x$overall_similarity))
  cat("  Similarity matrix (rows = Series A, cols = Series B):\n")
  print(round(x$similarity_matrix, 3))
  cat("\n  Best matches:\n")
  print(x$best_matches[, c("phenotype_A", "name_A", "best_match_B",
                            "name_B", "similarity")], row.names = FALSE)
  invisible(x)
}

#' @export
plot.AnomPhenoComparison <- function(x, ...) {
  sm <- x$similarity_matrix
  df <- expand.grid(
    A = rownames(sm),
    B = colnames(sm),
    stringsAsFactors = FALSE
  )
  df$similarity <- as.vector(sm)
  df$A <- factor(df$A, levels = rownames(sm))
  df$B <- factor(df$B, levels = colnames(sm))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$B, y = .data$A,
                                    fill = .data$similarity)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$similarity)),
                       size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_gradientn(
      colours = c("#1a237e", "#1565c0", "#42a5f5", "#e3f2fd", "#ffffff",
                  "#fff8e1", "#ffb300", "#e65100"),
      limits = c(-1, 1), name = "Similarity"
    ) +
    ggplot2::labs(
      title    = "Cross-Series Phenotype Similarity",
      subtitle = sprintf("Method: %s | Overall: %.3f",
                         x$method, x$overall_similarity),
      x = "Series B Phenotypes", y = "Series A Phenotypes"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold"),
      axis.text    = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

.adjusted_rand_index <- function(true_vec, pred_vec) {
  n    <- length(true_vec)
  tbl  <- table(true_vec, pred_vec)
  a    <- sum(choose(tbl, 2))
  b    <- sum(choose(rowSums(tbl), 2)) - a
  c_v  <- sum(choose(colSums(tbl), 2)) - a
  d    <- choose(n, 2) - a - b - c_v
  denom <- (a + b) * (a + c_v) / choose(n, 2)
  if (choose(n, 2) == 0) return(0)
  numerator <- a - denom
  denominator <- (choose(n, 2) - denom +
                    choose(n, 2) - denom) / 2 - denom + choose(n, 2)
  if (denominator == 0) return(1)
  expected <- (a + b) * (a + c_v) / choose(n, 2)
  max_val  <- ((a + b) + (a + c_v)) / 2
  if (max_val - expected == 0) return(1)
  (a - expected) / (max_val - expected)
}

.bhattacharyya_coef <- function(X, Y) {
  if (nrow(X) == 0 || nrow(Y) == 0) return(0)
  mu_x <- colMeans(X, na.rm = TRUE)
  mu_y <- colMeans(Y, na.rm = TRUE)

  if (nrow(X) == 1L || nrow(Y) == 1L) {
    denom <- sqrt(sum(mu_x^2)) * sqrt(sum(mu_y^2))
    if (denom < 1e-10) return(0)
    sim <- sum(mu_x * mu_y) / denom
    return(pmax(0, pmin(1, (sim + 1) / 2)))
  }

  var_x <- apply(X, 2, stats::var); var_x[is.na(var_x) | var_x < 1e-10] <- 1e-10
  var_y <- apply(Y, 2, stats::var); var_y[is.na(var_y) | var_y < 1e-10] <- 1e-10
  sig   <- (var_x + var_y) / 2
  d_mah <- (1/8) * sum((mu_x - mu_y)^2 / sig)
  d_cov <- 0.5 * sum(log(sig) - 0.5 * (log(var_x) + log(var_y)))
  bc    <- exp(-(d_mah + d_cov))
  pmax(0, pmin(1, bc))
}

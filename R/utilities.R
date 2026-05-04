
#' Simulate Time Series with Controlled Anomaly Phenotypes
#'
#' @description
#' Generates synthetic univariate time series with **known, controlled anomaly
#' phenotypes** injected at non-overlapping positions. The function is designed
#' for three primary use cases:
#'
#' \enumerate{
#'   \item **Benchmarking** -- generate series where ground truth is known, then
#'     evaluate how accurately [pheno_diagnose()] recovers the structure.
#'   \item **Power analysis** -- vary amplitude, duration, or noise to find the
#'     detection boundary for each phenotype type.
#'   \item **Reproducible research** -- specify `seed` to ensure exact
#'     replicability across platforms.
#' }
#'
#' ## Phenotype Types
#' \describe{
#'   \item{`"spike"`}{Gaussian-shaped sharp peak (high excess kurtosis, short).}
#'   \item{`"depression"`}{Gaussian-shaped trough (negative spike).}
#'   \item{`"plateau"`}{Sustained level shift (low oscillation, long duration).}
#'   \item{`"oscillation"`}{Sinusoidal burst (high HF energy ratio).}
#'   \item{`"drift"`}{Monotonic linear trend (high trend slope, high Hurst).}
#'   \item{`"sawtooth"`}{Repeating ramp pattern (high turning-point rate).}
#' }
#'
#' @param n Integer. Total series length. Default: 600.
#' @param phenotype_defs List of phenotype definition lists. Each element may
#'   contain: `type` (character), `amplitude` (numeric), `duration` (integer),
#'   `n_instances` (integer), `frequency` (numeric, for `"oscillation"` type).
#'   If `NULL`, uses a canonical 3-phenotype benchmark definition.
#' @param noise_sd Numeric. Standard deviation of Gaussian baseline noise.
#'   Default: 0.1.
#' @param baseline_freq Numeric. Period (in steps) of the sinusoidal baseline.
#'   Default: 30.
#' @param min_gap Integer. Minimum spacing between injected anomaly windows.
#'   Default: 20.
#' @param seed Integer. Random seed for reproducibility. Default: 42.
#'
#' @return A list with:
#'   \describe{
#'     \item{`x`}{Numeric vector: the simulated time series.}
#'     \item{`true_labels`}{Integer vector (length `n`): phenotype ID at each
#'       time step (0 = normal).}
#'     \item{`anomaly_idx`}{Integer vector: start indices of injected anomalies.}
#'     \item{`true_phenotype_id`}{Integer vector: phenotype ID per anomaly
#'       instance (parallel to `anomaly_idx`).}
#'     \item{`phenotype_defs`}{The phenotype definitions used.}
#'     \item{`seed`}{The seed used.}
#'   }
#'
#' @examples
#' sim <- pheno_simulate(n = 600, seed = 42)
#' result <- suppressMessages(pheno_diagnose(sim$x, q = 2.5))
#' if (!is.null(result)) {
#'   ev <- pheno_evaluate(result, sim$true_labels)
#'   print(ev)
#' }
#'
#' @export
pheno_simulate <- function(n             = 600L,
                            phenotype_defs = NULL,
                            noise_sd      = 0.1,
                            baseline_freq = 30,
                            min_gap       = 20L,
                            seed          = 42L) {
  n       <- as.integer(n)
  min_gap <- as.integer(min_gap)
  set.seed(seed)

  if (is.null(phenotype_defs)) {
    phenotype_defs <- list(
      list(type = "spike",       amplitude = 6.0, duration = 5L,  n_instances = 3L),
      list(type = "plateau",     amplitude = 3.0, duration = 15L, n_instances = 2L),
      list(type = "oscillation", amplitude = 2.5, duration = 10L, frequency = 0.4,
           n_instances = 2L)
    )
  }

  x <- sin(2 * pi * seq_len(n) / baseline_freq) + rnorm(n, sd = noise_sd)

  true_labels       <- rep(0L, n)
  anomaly_idx       <- integer(0)
  true_phenotype_id <- integer(0)
  occupied          <- integer(0)

  for (ph in seq_along(phenotype_defs)) {
    def    <- phenotype_defs[[ph]]
    type   <- .null_coalesce(def$type,        "spike")
    amp    <- .null_coalesce(def$amplitude,    4.0)
    dur    <- as.integer(.null_coalesce(def$duration, 10L))
    n_inst <- as.integer(.null_coalesce(def$n_instances, 2L))
    freq   <- .null_coalesce(def$frequency,    0.4)

    placed   <- 0L
    attempts <- 0L

    while (placed < n_inst && attempts < 500L) {
      attempts <- attempts + 1L
      st <- sample(seq(min_gap, n - dur - min_gap), 1L)
      en <- min(n, st + dur - 1L)

      buffer <- seq(max(1L, st - min_gap), min(n, en + min_gap))
      if (any(buffer %in% occupied)) next

      t_seg <- seq_len(en - st + 1L)
      slen  <- length(t_seg)

      x[st:en] <- x[st:en] + switch(type,
        spike = {
          pk <- ceiling(slen / 2)
          amp * exp(-((t_seg - pk)^2) / max(2, (slen / 4)^2))
        },
        depression = {
          pk <- ceiling(slen / 2)
          -amp * exp(-((t_seg - pk)^2) / max(2, (slen / 4)^2))
        },
        plateau = {
          amp + rnorm(slen, sd = noise_sd * 0.3)
        },
        oscillation = {
          amp * sin(2 * pi * freq * t_seg) + rnorm(slen, sd = noise_sd * 0.5)
        },
        drift = {
          amp * (t_seg / slen) + rnorm(slen, sd = noise_sd * 0.5)
        },
        sawtooth = {
          period <- max(3, ceiling(slen / 4))
          amp * ((t_seg - 1) %% period) / period
        },
        stop(sprintf("Unknown phenotype type: '%s'", type))
      )

      true_labels[st:en] <- ph
      anomaly_idx        <- c(anomaly_idx, st)
      true_phenotype_id  <- c(true_phenotype_id, ph)
      occupied           <- c(occupied, st:en)
      placed             <- placed + 1L
    }

    if (placed < n_inst)
      warning(sprintf(
        "Phenotype %d ('%s'): only %d of %d instances placed (try smaller n_instances or larger n).",
        ph, type, placed, n_inst))
  }

  ord <- order(anomaly_idx)

  list(
    x                = x,
    true_labels      = true_labels,
    anomaly_idx      = anomaly_idx[ord],
    true_phenotype_id = true_phenotype_id[ord],
    phenotype_defs   = phenotype_defs,
    seed             = seed
  )
}

#' Compute a Continuous Anomaly Severity Score Series
#'
#' @description
#' Returns a numeric vector of length `length(result$x)` where each anomaly
#' segment is assigned a continuous score and non-anomaly regions are zero.
#' Three scoring modes:
#'
#' \describe{
#'   \item{`"zscore"`}{Raw global z-score of each segment mean. Reflects how
#'     far the anomaly is from the overall series distribution.}
#'   \item{`"confidence"`}{Phenotype confidence score (0-1). Reflects how
#'     unambiguously the instance belongs to its assigned phenotype.}
#'   \item{`"combined"`}{Product of absolute z-score and confidence: a holistic
#'     score that is large only when the anomaly is both severe and clearly
#'     phenotyped. Recommended for monitoring dashboards.}
#' }
#'
#' When anomaly segments overlap (rare with default `gap` settings), the
#' maximum score across overlapping assignments is used.
#'
#' @param result An `AnomPhenoResult` from [pheno_diagnose()].
#' @param type Character. Scoring mode: `"zscore"` (default), `"confidence"`,
#'   or `"combined"`.
#' @param normalise Logical. If `TRUE`, scale scores to \[0, 1\] by dividing
#'   by the maximum observed score. Default: `FALSE`.
#'
#' @return Named numeric vector of length `length(result$x)`.
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result)) {
#'   scores <- pheno_score(result, type = "combined", normalise = TRUE)
#'   plot(scores, type = "l", ylab = "Anomaly score", xlab = "Time index")
#' }
#'
#' @export
pheno_score <- function(result, type = c("combined", "zscore", "confidence"),
                         normalise = FALSE) {
  if (!inherits(result, "AnomPhenoResult"))
    stop("'result' must be an AnomPhenoResult object.")
  type     <- match.arg(type)
  x        <- result$x
  n        <- length(x)
  inst_tab <- result$instance_table
  seg_len  <- result$seg_len

  scores <- numeric(n)

  for (i in seq_len(nrow(inst_tab))) {
    st       <- inst_tab$anomaly_idx[i]
    en       <- min(n, st + seg_len - 1L)
    gz       <- abs(inst_tab$global_zscore[i])
    conf     <- inst_tab$phenotype_confidence[i]
    score_i  <- switch(type,
      zscore     = gz,
      confidence = conf,
      combined   = gz * conf
    )
    scores[st:en] <- pmax(scores[st:en], score_i)
  }

  if (normalise) {
    mx <- max(scores)
    if (mx > 0) scores <- scores / mx
  }

  scores
}

#' Coerce an AnomPhenoResult to a Tidy Data Frame
#'
#' @description
#' Converts an `AnomPhenoResult` object into a **tidy data frame** with one
#' row per time step, making results directly compatible with
#' `dplyr`, `ggplot2`, or CSV export workflows.
#'
#' @param x An `AnomPhenoResult` from [pheno_diagnose()].
#' @param score_type Character. Score type passed to `pheno_score()`. One of
#'   `"combined"` (default), `"zscore"`, `"confidence"`.
#' @param ... Ignored (for S3 consistency).
#'
#' @return A `data.frame` with columns:
#'   \describe{
#'     \item{`time`}{Integer time index.}
#'     \item{`value`}{Original series value.}
#'     \item{`is_anomaly`}{Logical.}
#'     \item{`anomaly_score`}{Continuous score (0 = normal).}
#'     \item{`phenotype_id`}{Integer, `NA` for normal points.}
#'     \item{`phenotype_name`}{Character, `NA` for normal points.}
#'     \item{`phenotype_confidence`}{Numeric, `NA` for normal points.}
#'     \item{`severity`}{Character, `NA` for normal points.}
#'   }
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
#' if (!is.null(result)) {
#'   df <- as.data.frame(result)
#'   anomaly_rows <- df[df$is_anomaly, ]
#'   print(head(anomaly_rows))
#' }
#'
#' @export
as.data.frame.AnomPhenoResult <- function(x, score_type = "combined", ...) {
  result   <- x
  n        <- length(result$x)
  inst_tab <- result$instance_table
  ph_df    <- result$phenotypes
  seg_len  <- result$seg_len
  scores   <- pheno_score(result, type = score_type)

  out <- data.frame(
    time                 = seq_len(n),
    value                = result$x,
    is_anomaly           = FALSE,
    anomaly_score        = scores,
    phenotype_id         = NA_integer_,
    phenotype_name       = NA_character_,
    phenotype_confidence = NA_real_,
    severity             = NA_character_,
    stringsAsFactors     = FALSE
  )

  for (i in seq_len(nrow(inst_tab))) {
    inst <- inst_tab[i, ]
    st   <- inst$anomaly_idx
    en   <- min(n, st + seg_len - 1L)
    ph   <- inst$phenotype_id
    ph_row <- ph_df[ph_df$phenotype_id == ph, ]

    out$is_anomaly[st:en]           <- TRUE
    out$phenotype_id[st:en]         <- ph
    out$phenotype_name[st:en]       <- ph_row$name
    out$phenotype_confidence[st:en] <- inst$phenotype_confidence
    out$severity[st:en]             <- ph_row$severity
  }

  out
}

.null_coalesce <- function(x, default) if (is.null(x)) default else x

#' Benchmark Anomaly Dataset
#'
#' @description
#' A synthetic benchmark time series with **six injected anomalies of three
#' distinct types**, designed for illustrating and validating the AnomalyPheno
#' phenotyping pipeline. The series has n = 600 observations with a sinusoidal
#' seasonal pattern (period = 30), Gaussian noise, and three phenotype classes:
#'
#' \describe{
#'   \item{**Phenotype A** -- Sharp Transient Spike}{
#'     A single-observation amplitude excursion, 6-8 SD above baseline.
#'     Occurs at indices 50, 300, and 500. Morphologically: very short duration,
#'     high kurtosis, symmetric onset/recovery slope.}
#'   \item{**Phenotype B** -- Sustained Elevation Drift}{
#'     A 10-12 step window of elevated values (mean shift +3 SD).
#'     Occurs at indices 120-130 and 420-432. Morphologically: long duration,
#'     low kurtosis, gradual onset and recovery.}
#'   \item{**Phenotype C** -- Oscillatory Burst}{
#'     A 4-step alternating-sign pattern, +/-3 SD.
#'     Occurs at index 220. Morphologically: high-frequency energy, asymmetric,
#'     biphasic shape.}
#' }
#'
#' @format A list with three elements:
#' \describe{
#'   \item{`x`}{Numeric vector of length 600. The time series.}
#'   \item{`true_labels`}{Integer vector of length 600. 0 = normal,
#'     1 = Phenotype A (spike), 2 = Phenotype B (drift), 3 = Phenotype C (burst).}
#'   \item{`anomaly_idx`}{Integer vector. True start indices of each anomaly
#'     interval.}
#' }
#'
#' @examples
#' data(anomaly_benchmark)
#' plot(anomaly_benchmark$x, type = "l",
#'      main = "AnomalyPheno Benchmark Dataset")
#'
#' result <- pheno_diagnose(anomaly_benchmark$x, q = 2.5)
#' plot(result)
#'
"anomaly_benchmark"

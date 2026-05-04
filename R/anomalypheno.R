
#' Diagnose Anomaly Phenotypes in a Time Series
#'
#' @description
#' The single entry-point for the **AnomalyPheno** pipeline. Orchestrates four
#' layers automatically:
#'
#' 1. **Detection** -- flags anomaly intervals using the Spectral Residual (SR)
#'    or Local Outlier Factor (LOF) engine.
#' 2. **Feature extraction** -- computes a **27-dimensional** signature per
#'    anomaly covering morphological, spectral, context, severity, and
#'    complexity domains.
#' 3. **Multi-Resolution Mapper** -- runs topological phenotyping at three
#'    resolution scales and combines assignments by majority vote, producing
#'    a **stability score** for each instance.
#' 4. **Characterisation** -- assigns auto-names, canonical templates, severity
#'    classes, phenotype confidence scores, and a feature importance ranking.
#'
#' @param x Numeric vector. The input time series.
#' @param q Numeric. Detection threshold. SR mode: q-sigma multiplier
#'   (typical: 2.5-3.5). LOF mode: LOF score cutoff (typical: 1.5-3.0).
#'   Default: 3.0.
#' @param engine Character. Detection engine: `"sr"` (Spectral Residual,
#'   default) or `"lof"` (Local Outlier Factor sliding window).
#' @param seg_len Integer. Length of each anomaly segment window. Default: 10.
#' @param gap Integer. Minimum index gap between distinct anomaly intervals.
#'   Default: 3.
#' @param n_intervals Integer or `NULL`. Base Mapper cover intervals. If `NULL`
#'   (default), auto-selected via Sturges rule.
#' @param overlap Numeric. Mapper interval overlap fraction. Default: 0.5.
#'
#' @return An `AnomPhenoResult` object, or `NULL` if no anomalies are detected.
#'   Print, summary, and plot methods are available.
#'
#' @seealso [detect_anomaly_intervals()], [extract_anomaly_features()],
#'   [mapper_phenotype()], [characterise_phenotypes()], [pheno_predict()],
#'   [plot_phenotype_timeline()], [plot_phenotype_gallery()],
#'   [plot_phenotype_profile()], [plot_phenotype_landscape()]
#'
#' @examples
#' data(anomaly_benchmark)
#' result <- pheno_diagnose(anomaly_benchmark$x, q = 2.5)
#' print(result)
#' plot(result)                          # canonical template gallery
#' plot(result, type = "timeline")      # phenotype-coloured series
#' plot(result, type = "profile")       # feature discriminability chart
#' plot(result, type = "landscape")     # Mapper topological graph
#'
#' # Classify a new anomaly
#' pred <- pheno_predict(result, new_segment = c(0.1, 7.0, 0.2))
#'
#' @export
pheno_diagnose <- function(x, q = 3.0, engine = c("sr", "lof"),
                            seg_len = 10L, gap = 3L,
                            n_intervals = NULL, overlap = 0.5) {
  engine  <- match.arg(engine)
  seg_len <- as.integer(seg_len)
  gap     <- as.integer(gap)

  message("AnomalyPheno [1/4] Detecting anomaly intervals (engine: ",
          toupper(engine), ")...")
  idx <- detect_anomaly_intervals(x, q = q, gap = gap, engine = engine)
  if (length(idx) == 0L) {
    message("AnomalyPheno: No anomaly intervals detected. ",
            "Try lowering 'q' or switching engine.")
    return(NULL)
  }
  message(sprintf("  -> %d anomaly interval(s) detected.", length(idx)))

  message("AnomalyPheno [2/4] Extracting 27-dimensional feature vectors...")
  feat_mat <- extract_anomaly_features(x, idx, seg_len = seg_len)
  if (is.null(feat_mat)) {
    message("AnomalyPheno: Feature extraction returned NULL.")
    return(NULL)
  }

  if (nrow(feat_mat) < 2L) {
    message("AnomalyPheno: Too few anomaly instances for phenotyping (need >= 2). ",
            "Try lowering 'q'.")
    return(NULL)
  }

  message("AnomalyPheno [3/4] Running Multi-Resolution Mapper phenotyping...")
  mp <- mapper_phenotype(feat_mat, n_intervals = n_intervals, overlap = overlap)
  message(sprintf("  -> %d phenotype(s) discovered (Mapper nodes: %d).",
                  mp$n_phenotypes, length(mp$nodes)))

  message("AnomalyPheno [4/4] Characterising phenotypes...")
  result <- characterise_phenotypes(x, idx, feat_mat, mp, seg_len = seg_len)
  message("  -> Done.")
  return(result)
}

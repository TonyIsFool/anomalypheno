library(testthat)
library(AnomalyPheno)

# -- Helpers ------------------------------------------------------------------

make_test_series <- function(seed = 42) {
  set.seed(seed)
  x <- sin(2 * pi * seq_len(300) / 30) + rnorm(300, sd = 0.1)
  x[50]      <- 6
  x[150:160] <- x[150:160] + 3
  x[250:253] <- c(2, -2, 2, -2)
  x
}

make_rich_series <- function(seed = 7) {
  set.seed(seed)
  x <- sin(2 * pi * seq_len(500) / 30) + rnorm(500, sd = 0.15)
  x[50] <- 8; x[51] <- 7.5
  x[100:110] <- x[100:110] + 4
  x[200:202] <- c(4, -4, 4)
  x[300] <- -7
  x[400:405] <- x[400:405] - 3
  x[450] <- 6
  x
}

# -- Layer 1: Detection -------------------------------------------------------

test_that("SR engine returns integer vector", {
  idx <- detect_anomaly_intervals(make_test_series(), q = 2.5)
  expect_type(idx, "integer")
  expect_gte(length(idx), 1L)
})

test_that("SR engine finds spike near index 50", {
  idx <- detect_anomaly_intervals(make_test_series(), q = 3.0)
  expect_true(any(abs(idx - 50L) <= 5L))
})

test_that("LOF engine returns integer vector", {
  idx <- detect_anomaly_intervals(make_test_series(), q = 2.0, engine = "lof")
  expect_type(idx, "integer")
})

test_that("Flat series at extreme q returns zero anomalies", {
  set.seed(99)
  expect_equal(length(detect_anomaly_intervals(rnorm(100), q = 99)), 0L)
})

# -- Layer 2: Feature Extraction (27-dim) -------------------------------------

test_that("Feature matrix has exactly 27 columns", {
  x   <- make_test_series()
  idx <- detect_anomaly_intervals(x, q = 3.0)
  fm  <- extract_anomaly_features(x, idx, seg_len = 10L)
  expect_equal(ncol(fm), 27L)
  expect_equal(nrow(fm), length(idx))
})

test_that("Feature values are finite", {
  x   <- make_test_series()
  idx <- detect_anomaly_intervals(x, q = 3.0)
  fm  <- extract_anomaly_features(x, idx, seg_len = 10L)
  expect_true(all(is.finite(fm)))
})

test_that("peak_abs is non-negative", {
  x   <- make_test_series()
  idx <- detect_anomaly_intervals(x, q = 3.0)
  fm  <- extract_anomaly_features(x, idx, seg_len = 10L)
  expect_true(all(fm[, "peak_abs"] >= 0))
})

# -- Layer 3: Mapper ----------------------------------------------------------

test_that("mapper_phenotype returns AnomalyMapper", {
  inp <- { x <- make_rich_series(); idx <- detect_anomaly_intervals(x, q = 2.0)
            fm <- extract_anomaly_features(x, idx, seg_len = 10L)
            list(x=x,idx=idx,fm=fm) }
  if (is.null(inp$fm) || nrow(inp$fm) < 2) skip("Too few instances")
  mp <- mapper_phenotype(inp$fm)
  expect_s3_class(mp, "AnomalyMapper")
})

test_that("Mapper phenotype_id has length == n_instances", {
  x   <- make_rich_series()
  idx <- detect_anomaly_intervals(x, q = 2.0)
  fm  <- extract_anomaly_features(x, idx, seg_len = 10L)
  if (is.null(fm) || nrow(fm) < 2) skip("Too few instances")
  mp <- mapper_phenotype(fm)
  expect_equal(length(mp$phenotype_id), nrow(fm))
})

# -- Layer 4: pheno_diagnose (full pipeline) -----------------------------------

test_that("pheno_diagnose returns AnomPhenoResult", {
  result <- suppressMessages(pheno_diagnose(make_rich_series(), q = 2.0))
  expect_s3_class(result, "AnomPhenoResult")
})

test_that("AnomPhenoResult has required fields", {
  result <- suppressMessages(pheno_diagnose(make_test_series(), q = 2.5))
  if (is.null(result)) skip("No anomalies detected")
  expect_true(all(c("phenotypes","instance_table","canonical_templates",
                     "feature_importance","mapper","feat_mat") %in% names(result)))
})

test_that("Instance confidence is in [0, 1]", {
  result <- suppressMessages(pheno_diagnose(make_test_series(), q = 2.5))
  if (is.null(result)) skip("No anomalies")
  conf <- result$instance_table$phenotype_confidence
  expect_true(all(conf >= 0 & conf <= 1))
})

# -- Benchmark dataset --------------------------------------------------------

test_that("anomaly_benchmark loads correctly", {
  data(anomaly_benchmark)
  expect_true(is.list(anomaly_benchmark))
  expect_true("x" %in% names(anomaly_benchmark))
  expect_gte(length(anomaly_benchmark$x), 500L)
})

test_that("pheno_diagnose runs on benchmark data", {
  data(anomaly_benchmark)
  result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
  expect_s3_class(result, "AnomPhenoResult")
  expect_gte(result$mapper$n_phenotypes, 1L)
})

# -- Evaluation ---------------------------------------------------------------

test_that("pheno_evaluate returns AnomPhenoEval", {
  data(anomaly_benchmark)
  result <- suppressMessages(pheno_diagnose(anomaly_benchmark$x, q = 2.5))
  if (is.null(result)) skip("No result")
  gt <- anomaly_benchmark$ground_truth_labels
  if (is.null(gt)) skip("No ground truth")
  ev <- pheno_evaluate(result, gt)
  expect_s3_class(ev, "AnomPhenoEval")
  expect_true(ev$purity >= 0 && ev$purity <= 1)
  expect_true(ev$ari >= -1 && ev$ari <= 1)
})

test_that("pheno_cluster_quality returns valid metrics", {
  result <- suppressMessages(pheno_diagnose(make_rich_series(), q = 2.0))
  if (is.null(result) || result$mapper$n_phenotypes < 2) skip("Too few phenotypes")
  qlt <- pheno_cluster_quality(result)
  expect_true(is.list(qlt))
  expect_true(all(c("silhouette_mean","calinski_harabasz","davies_bouldin",
                     "n_phenotypes","n_instances") %in% names(qlt)))
})

# -- Advanced: prune, rename, track -------------------------------------------

test_that("pheno_prune returns AnomPhenoResult", {
  result <- suppressMessages(pheno_diagnose(make_rich_series(), q = 2.0))
  if (is.null(result)) skip("No result")
  pruned <- suppressMessages(pheno_prune(result, min_instances = 2L))
  expect_s3_class(pruned, "AnomPhenoResult")
  expect_lte(pruned$mapper$n_phenotypes, result$mapper$n_phenotypes)
})

test_that("pheno_rename returns AnomPhenoResult with non-empty names", {
  result <- suppressMessages(pheno_diagnose(make_rich_series(), q = 2.0))
  if (is.null(result)) skip("No result")
  renamed <- pheno_rename(result)
  expect_s3_class(renamed, "AnomPhenoResult")
  expect_true(all(nchar(renamed$phenotypes$name) > 0))
})

test_that("pheno_track returns AnomPhenoTrack", {
  result <- suppressMessages(pheno_diagnose(make_rich_series(), q = 2.0))
  if (is.null(result)) skip("No result")
  trk <- pheno_track(result, window_size = 80L, step = 20L)
  expect_s3_class(trk, "AnomPhenoTrack")
  expect_equal(ncol(trk$track_matrix), result$mapper$n_phenotypes)
})

# -- Inference: auto-tune, transition, regime ---------------------------------

test_that("pheno_auto_tune returns AnomPhenoTune with best_q in grid", {
  g    <- c(2.0, 2.5, 3.0)
  tune <- suppressMessages(pheno_auto_tune(make_test_series(), q_grid = g, verbose = FALSE))
  expect_s3_class(tune, "AnomPhenoTune")
  expect_true(tune$best_q %in% g)
})

test_that("pheno_transition matrix is row-stochastic", {
  result <- suppressMessages(pheno_diagnose(make_test_series(), q = 2.0))
  if (is.null(result) || nrow(result$instance_table) < 2) skip("Too few instances")
  tr <- pheno_transition(result)
  rs <- rowSums(tr$transition_matrix)
  expect_true(all(abs(rs - 1) < 1e-9 | rs == 0))
})

test_that("pheno_regime returns AnomPhenoRegime", {
  result <- suppressMessages(pheno_diagnose(make_test_series(), q = 2.0))
  if (is.null(result) || nrow(result$instance_table) < 3) skip("Too few instances")
  reg <- pheno_regime(result, h = 2)
  expect_s3_class(reg, "AnomPhenoRegime")
  expect_gte(length(reg$cusum_pos), 1L)
})

# -- Probabilistic: GMM, FPCA, permtest ---------------------------------------

test_that("pheno_gmm posterior rows sum to 1", {
  result <- suppressMessages(pheno_diagnose(make_test_series(), q = 2.0))
  if (is.null(result) || nrow(result$instance_table) < 4) skip("Too few instances")
  gmm <- suppressMessages(pheno_gmm(result, k_max = 3L, n_init = 2L, seed = 1L, verbose = FALSE))
  expect_s3_class(gmm, "AnomPhenoGMM")
  expect_true(all(abs(rowSums(gmm$posterior) - 1) < 1e-6))
})

test_that("pheno_fpca variance explained <= 1", {
  result <- suppressMessages(pheno_diagnose(make_test_series(), q = 2.0))
  if (is.null(result) || nrow(result$instance_table) < 3) skip("Too few instances")
  fp <- pheno_fpca(result, n_comp = 2L)
  expect_s3_class(fp, "AnomPhenoFPCA")
  expect_lte(sum(fp$var_explained), 1 + 1e-9)
})

test_that("pheno_permtest p_value in [0, 1]", {
  result <- suppressMessages(pheno_diagnose(make_test_series(), q = 2.0))
  if (is.null(result) || result$mapper$n_phenotypes < 2 ||
      nrow(result$instance_table) < 4) skip("Too few")
  pt <- suppressMessages(pheno_permtest(result, B = 49L, seed = 1L, verbose = FALSE))
  expect_s3_class(pt, "AnomPhenoPermTest")
  expect_gte(pt$p_value, 0); expect_lte(pt$p_value, 1)
})

# -- Simulation & utilities ---------------------------------------------------

test_that("pheno_simulate returns valid series with ground truth", {
  sim <- pheno_simulate(n = 300L, seed = 1L)
  expect_true(is.list(sim))
  expect_true(all(c("x", "anomaly_idx", "true_phenotype_id") %in% names(sim)))
  expect_equal(length(sim$anomaly_idx), length(sim$true_phenotype_id))
})

test_that("pheno_predict classifies new segment", {
  result <- suppressMessages(pheno_diagnose(make_rich_series(), q = 2.0))
  if (is.null(result)) skip("No result")
  pred <- pheno_predict(result, new_segment = c(0, 3, 5, 3, 0))
  expect_true(is.list(pred))
  expect_named(pred, c("phenotype_id", "phenotype_name", "confidence"),
               ignore.order = TRUE)
  expect_gte(pred$confidence, 0); expect_lte(pred$confidence, 1)
})

test_that("as.data.frame.AnomPhenoResult returns data.frame", {
  result <- suppressMessages(pheno_diagnose(make_rich_series(), q = 2.0))
  if (is.null(result)) skip("No result")
  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
  expect_true("phenotype_id" %in% colnames(df))
})

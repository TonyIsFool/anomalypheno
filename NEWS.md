# News for AnomalyPheno

## AnomalyPheno 0.1.0

### New functionality
* `pheno_diagnose()`: unified pipeline API with SR and LOF detection engines
* `detect_anomaly_intervals()`: SR and LOF detection, gap-merge logic
* `extract_anomaly_features()`: 27-dimensional feature space (22 original +
  5 complexity features: approximate entropy, turning-point rate, trend slope,
  Hurst exponent approximation, zero-crossing rate)
* `mapper_phenotype()`: Multi-Resolution Mapper (three scales, majority-vote
  aggregation, per-instance stability scores via Sturges auto-selection)
* `characterise_phenotypes()`: phenotype confidence scores, feature importance
  F-statistic table, Hurst/ApEn in phenotype summary
* `pheno_predict()`: prospective classification of new anomaly segments
* `pheno_evaluate()`: external validation — purity, Adjusted Rand Index,
  per-phenotype precision/recall/F1, confusion matrix
* `pheno_bootstrap()`: bootstrap Jaccard stability analysis (B resamples,
  Gaussian noise perturbation, per-phenotype stability scores)
* `pheno_compare()`: cross-series phenotype comparison via Bhattacharyya
  coefficient, cosine similarity, or negative Euclidean distance; includes
  heatmap visualisation via `plot.AnomPhenoComparison()`
* Four plot functions: `plot_phenotype_timeline()`, `plot_phenotype_gallery()`,
  `plot_phenotype_profile()`, `plot_phenotype_landscape()`
* `plot.AnomPhenoResult()` S3 dispatch with `type` argument
* Built-in dataset `anomaly_benchmark`: 600-point series with 3 known phenotype
  classes (Transient Spike, Sustained Drift, Oscillatory Burst)

### Architecture
* 5-file structure: `anomalypheno.R`, `core_detection.R`,
  `core_classification.R`, `core_evaluation.R`, `visualise.R`
* Pure-R implementation (no C++ / external TDA library dependencies)
* Full `testthat` test suite (58+ tests)

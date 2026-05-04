# AnomalyPheno

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)]()
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**AnomalyPheno** is an R package for unsupervised anomaly phenotyping in univariate time series. It extends conventional anomaly detection by automatically classifying detected anomalies into distinct, interpretable types (phenotypes) using topological data analysis.

The pipeline operates in two stages:

1. **Detection** — Identify anomalous intervals via Spectral Residual or Local Outlier Factor engines.
2. **Phenotyping** — Cluster anomaly instances into phenotypes using the Multi-Resolution Mapper algorithm, with automated characterisation and confidence scoring.

---

## The Core Analogy

| Medical Diagnostics | AnomalyPheno |
|---|---|
| ECG monitor flags abnormal heart rhythm | SR engine flags anomaly interval |
| Cardiologist extracts waveform features | 27-dim feature extraction |
| Diagnosis: "Type II AV Block" | Auto-named phenotype: "Oscillatory Burst" |
| Treatment protocol by diagnosis type | Downstream action by phenotype |

---

## Installation

```r
# From GitHub
devtools::install_github("TonyIsFool/anomalypheno")

# (CRAN submission pending)
```

---

## Quick Start

```r
library(AnomalyPheno)

data(anomaly_benchmark)

result <- pheno_diagnose(anomaly_benchmark$x, q = 2.5)
print(result)

# Visualisations
plot(result, type = "gallery")
plot_phenotype_timeline(result)
plot_phenotype_profile(result)
plot_phenotype_landscape(result)

# Classify a new anomaly segment
pred <- pheno_predict(result, new_segment = c(0.1, 7.5, 0.2))
cat("Phenotype:", pred$phenotype_name, "| Confidence:", pred$confidence)
```

---

## Comparison with Existing Packages

| Package | Detects | Cannot do |
|---|---|---|
| `anomalize` | Point anomalies (value spikes) | Classify anomaly types |
| `tsmp` | Shape anomalies (Matrix Profile) | Explain what the shape means |
| `ICSOutlier` | Multivariate outliers | Handle temporal structure |
| **AnomalyPheno** | All of the above + **phenotype classification** | — |

### Why Mapper?

Unlike k-means or hierarchical clustering, the Mapper algorithm preserves the **topological structure** of the feature space, revealing branching phenotype relationships that flat clustering methods miss.

---

## Visualisations

### 1. Phenotype-Coloured Timeline
```r
plot_phenotype_timeline(result)
```

### 2. Canonical Template Gallery
```r
plot(result, type = "gallery")
```

### 3. Feature Profile Chart
```r
plot_phenotype_profile(result)
```

### 4. Mapper Landscape Graph
```r
plot(result, type = "landscape")
```

---

## References

- Ren, H., et al. (2019). Time-series anomaly detection service at Microsoft. *KDD '19*. <https://doi.org/10.1145/3292500.3330680>
- Singh, G., Mémoli, F., & Carlsson, G. (2007). Topological methods for the analysis of high dimensional data sets and 3D object recognition. *Eurographics Symposium on Point-Based Graphics*.
- Carlsson, G. (2009). Topology and data. *Bulletin of the American Mathematical Society*, 46(2), 255–308.

---

## Citation

```bibtex
@misc{anomalypheno2024,
  title  = {{AnomalyPheno}: Unsupervised Anomaly Phenotyping for Time Series
             via Topological Data Analysis},
  author = {Lu, Tony},
  year   = {2024},
  note   = {R package version 0.1.0}
}
```

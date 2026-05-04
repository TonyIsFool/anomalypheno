# AnomalyPheno <img src="man/figures/logo.png" align="right" height="120"/>

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)]()
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![JCGS Submission](https://img.shields.io/badge/Target-JCGS-orange)]()
<!-- badges: end -->

> **"Don't just detect anomalies. Diagnose them."**

**AnomalyPheno** is the first R package that treats time series anomaly
detection as a *two-stage medical diagnostic problem*:

1. **Detection** — *Is there something wrong?*  
2. **Phenotyping** — *What type of abnormality is it?*

Every existing package stops at stage 1. AnomalyPheno does both, using the
**Mapper** topological data analysis algorithm (Singh et al., 2007) to discover
anomaly *phenotypes* — interpretable, reproducible categories of anomalous
behaviour — directly from the data, with no human labels required.

---

## The Core Analogy

| Medical Diagnostics | AnomalyPheno |
|---|---|
| ECG monitor flags abnormal heart rhythm | SR engine flags anomaly interval |
| Cardiologist extracts waveform features | 22-dim feature extraction |
| Diagnosis: "Type II AV Block" | Auto-named phenotype: "Oscillatory Burst" |
| Treatment protocol by diagnosis type | Downstream action by phenotype |

---

## Installation

```r
# From local source
devtools::install("path/to/Package4")

# (CRAN submission pending)
```

---

## Quick Start

```r
library(AnomalyPheno)

data(anomaly_benchmark)          # built-in 3-phenotype benchmark

result <- pheno_diagnose(anomaly_benchmark$x, q = 2.5)
print(result)
#> ── AnomalyPheno Diagnostic Report ────────────────────────────────
#>   Total anomaly instances : 6
#>   Phenotypes discovered   : 3
#>   Mapper nodes            : 12
#>
#>   Phenotype Summary:
#>  phenotype_id                       name count severity mean_amplitude mean_duration
#>             1    Transient Elevation Spike     3   Severe          6.43            10
#>             2  Sustained Elevation Plateau     2 Moderate          3.25            10
#>             3    Brief Oscillatory Burst     1   Severe          3.50            10

# Four complementary visualisations
plot(result, type = "gallery")      # canonical templates
plot_phenotype_timeline(result)     # time series coloured by phenotype
plot_phenotype_profile(result)      # feature radar per phenotype
plot_phenotype_landscape(result)    # Mapper topological graph

# Prospective classification of new anomaly
pred <- pheno_predict(result, new_segment = c(0.1, 7.5, 0.2))
cat("Phenotype:", pred$phenotype_name, "| Confidence:", pred$confidence)
```

---

## What Makes AnomalyPheno Novel?

### Problem with existing approaches

| Package | Detects | Cannot do |
|---|---|---|
| `anomalize` | Point anomalies (value spikes) | Classify anomaly types |
| `tsmp` | Shape anomalies (Matrix Profile) | Explain what the shape means |
| `ICSOutlier` | Multivariate outliers | Handle temporal structure |
| **AnomalyPheno** | All of the above + **phenotype classification** | — |

### The Mapper Advantage

Unlike k-means or hierarchical clustering, the Mapper algorithm preserves
the **topological structure** of the feature space:

- **K-means** assumes spherical clusters → misses anomaly "families" that share boundaries
- **Mapper** builds a graph where edges connect overlapping clusters → reveals branching phenotype relationships invisible to flat clustering

---

## The Four Visualisations

### 1. Phenotype-Coloured Timeline
The signature figure for any paper or report — shows when and what type:

```r
plot_phenotype_timeline(result)
```

### 2. Canonical Template Gallery  
The "face" of each phenotype — the medoid segment:

```r
plot(result, type = "gallery")
```

### 3. Feature Profile Chart
What makes each phenotype distinctive across 8 top features:

```r
plot_phenotype_profile(result)
```

### 4. Mapper Landscape Graph
The topological structure of the anomaly space:

```r
plot(result, type = "landscape")
```

---

## Algorithm Pipeline

```
Input: numeric time series x
          │
          ▼
[1] Spectral Residual Detection (Ren et al., 2019)
    → anomaly_idx: integer vector of interval starts
          │
          ▼
[2] 22-Dimensional Feature Extraction
    ┌─────────────────────────────────────────────────┐
    │ Morphological (8): peak_abs, onset_slope, AUC,  │
    │   duration, asymmetry, sharpness, ...           │
    │ Spectral     (5): dom_freq, spec_entropy, HF%   │
    │ Context      (5): pre/post mean & var, local_z  │
    │ Severity     (4): global_z, pct_rank, IQR_score │
    └─────────────────────────────────────────────────┘
          │
          ▼
[3] Mapper Topological Clustering (Singh et al., 2007)
    - Filter: PC1 of feature matrix
    - Cover: N overlapping intervals with δ overlap
    - Cluster: single-linkage per interval
    - Graph: connect nodes with shared members
    - Phenotypes: connected components
          │
          ▼
[4] Phenotype Characterisation
    - Medoid (canonical template)
    - Auto-name from dominant features
    - Severity: Mild / Moderate / Severe
    - Prospective classification via pheno_predict()
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
  author = {L, Tony},
  year   = {2024},
  note   = {R package version 0.1.0}
}
```

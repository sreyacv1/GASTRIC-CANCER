# Tables

**Table 1. Study cohorts and prognostic performance.** Discovery and external-validation cohorts for the transcriptomic prognostic signature. Two distinct effect measures are reported and must not be conflated: the median-split high-versus-low hazard ratio (unadjusted), and the age/stage-adjusted hazard ratio per standard deviation of the continuous risk score — the latter is the estimate carried into the meta-analysis.

| Cohort | Role | n | Events | Harrell C | HR high-vs-low, median split (95% CI) | Age/stage-adjusted HR per SD (95% CI) |
|---|---|---|---|---|---|---|
| TCGA-STAD | Discovery (training) | 383 | 156 | 0.61 (nested-CV) | — | — |
| GSE62254 (ACRG) | External validation | 300 | 152 | 0.61 | 1.90 (1.37–2.62) | 1.30 (1.10–1.54) |
| GSE15459 | External validation | 191 | 95 | 0.58 | 1.68 (1.11–2.54) | 1.20 (0.97–1.48) |
| GSE84437 | External validation | 431 | 207 | 0.53 | 1.11 (0.84–1.46) | 1.11 (0.97–1.27) |
| **Pooled (Hartung–Knapp, per SD)** | Meta-analysis | 922 | 454 | — | — | **1.19 (0.96–1.47), p=0.073, not significant** |

*CV, cross-validation; HR, hazard ratio; SD, standard deviation of the risk score; CI, confidence interval. The TCGA C-index is the leakage-free nested cross-validation estimate (apparent C 0.72). GSE84437 adjustment used age + T-stage (full TNM stage was unavailable in that series). The pooled estimate is a REML random-effects meta-analysis of the age/stage-adjusted per-SD log-hazard ratios with the Hartung–Knapp adjustment (I²=19.2%); its confidence interval includes 1, so the pooled effect is not statistically significant.*

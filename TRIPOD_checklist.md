# Supplementary Table S3 — TRIPOD Checklist

Prognostic prediction model: a 25-gene stromal/fibroblast signature and clinical
nomogram for overall survival in gastric cancer. Model **development** in TCGA-STAD
with internal optimism-corrected (bootstrap) validation and **external validation**
in independent GEO cohorts (ACRG/GSE62254, GSE15459, GSE84437). This is a Type 3
study (development + external validation) under the TRIPOD taxonomy.

| # | TRIPOD item | Addressed | Location |
|---|---|---|---|
| 1 | Title identifies study as developing/validating a prediction model | Yes | Title; Abstract |
| 2 | Abstract: objectives, data, methods, results, conclusions | Yes | Abstract |
| 3a | Background and rationale | Yes | §1 Introduction |
| 3b | Objectives | Yes | §1; end of Introduction |
| 4a | Source of data (development & validation) | Yes | §2.1–2.2; Data availability |
| 4b | Dates / study period | Yes (public cohorts, accessions dated) | §2.1; refs |
| 5a | Study setting / eligibility | Yes | §2.1 (cohorts, inclusion) |
| 5b | Details of treatments received | N/A (retrospective transcriptomic cohorts) | §2.1 |
| 6a | Outcome definition (overall survival) | Yes | §2.7 |
| 6b | Outcome assessment blinding | N/A (registry survival) | §2.7 |
| 7a | Predictors: definition and measurement | Yes | §2.5–2.7 (hub genes, Stage, Age, Grade, TMB) |
| 7b | Predictor assessment blinding | N/A | §2.7 |
| 8 | Sample size | Yes (EPV reported) | §2.7; §3.4 |
| 9 | Missing data handling | Yes (complete-case, no imputation; stated) | §2.7 |
| 10a | Predictor handling in analysis | Yes | §2.6–2.7 (LASSO; backward-AIC) |
| 10b | Model-building procedure | Yes (LASSO signature; backward-AIC clinical) | §2.6–2.7 |
| 10c | Model-updating (recalibration) at validation | Yes (out-of-sample, no re-fit) | §3.7 |
| 10d | Predictive performance measures | Yes (Harrell/Uno C, time-AUC, calibration, DCA) | §2.7; §3.4, §3.7 |
| 10e | Model comparison (added value) | Yes (ΔC, LRT, IDI, NRI vs staging) | §3.7 |
| 11 | Risk-group definition | Yes (median-split; continuous per-SD also reported) | §3.4 |
| 12 | Development vs validation differences | Yes (platform/cohort heterogeneity discussed) | §3.4; §4.1 |
| 13a | Participant flow / numbers | Yes | §2.1; §3.4 (n, events per cohort) |
| 13b | Participant characteristics | Yes | §2.1; §3.7 |
| 13c | Missing-data numbers | Yes (383→183 complete-case attrition stated) | §3.6 |
| 14a | Model specification (final) | Yes (25-gene signature; nomogram covariates) | §2.6–2.7; §3.5 |
| 14b | Model performance | Yes | §3.4, §3.7 |
| 15a | Full model / how to obtain predictions | Yes (coefficients in `results/`; nomogram) | §3.5; Code availability |
| 15b | Model presentation | Yes (nomogram; Fig 2) | §3.5 |
| 16 | Validation results (discrimination, calibration) | Yes (external C 0.53–0.61 across three cohorts, incl. GSE84437 null C=0.53; HK meta; DCA null) | §3.4, §3.7 |
| 17 | Comparison with prior models | Yes (Li et al. 2023 CAF signature; established biology) | §4 Discussion |
| 18 | Limitations | Yes (dedicated section incl. shared-data circularity) | §4.1 |
| 19a | Interpretation vs objectives | Yes | §4 Discussion |
| 19b | Overall interpretation | Yes | §4; Clinical implications callout |
| 20 | Implications for practice | Yes ("not a clinical tool"; AJCC standard of care) | Clinical implications callout |
| 21 | Supplementary information | Yes | Supplementary Tables S1–S3; `results/` |
| 22 | Funding | Yes | Declarations |

**Notes.** The signature and nomogram are reported as research instruments, not
validated clinical tools (item 20). The most important honestly-reported deviation
from a positive-model narrative is that external added value over AJCC stage was
decision-analytically negligible (ΔC +0.005; no decision-curve benefit; items 10e, 16),
and that model/signature agreement is partly internal-consistency on shared TCGA data
rather than independent replication (item 18).

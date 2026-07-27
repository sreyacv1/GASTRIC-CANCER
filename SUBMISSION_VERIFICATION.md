# Pre-Submission Verification Report

**Date:** 2026-07-27 · **Verdict: NO FALSIFICATION FOUND. All headline claims reproduce from source files.**

Every number below was re-derived from the result file or raw data *independently of the
manuscript text*, then compared to what the manuscript states.

## A. Prognostic signature — VERIFIED

| Claim (manuscript) | Source file | Source value | Match |
|---|---|---|---|
| Nested-CV Harrell C 0.611 | `nested_cv/performance.csv` | 0.611206 (CI 0.562–0.659) | ✓ |
| Nested-CV Uno C 0.573 | same | 0.572722 | ✓ |
| Integrated Brier 0.179 | same | 0.178892 | ✓ |
| λ.1se sensitivity 0.597 | same | 0.597446 | ✓ |

## B. External validation (n=922) — VERIFIED

| Cohort | n / events | C-index | HR (95% CI) | Match |
|---|---|---|---|---|
| ACRG/GSE62254 | 300 / 152 | 0.6079 | 1.8965 (1.3704–2.6245) | ✓ |
| GSE15459 | 191 / 95 | 0.5752 | 1.6757 (1.1067–2.5373) | ✓ |
| GSE84437 | 431 / 207 | 0.5297 | 1.1087 (0.8440–1.4565), ns | ✓ |
| Pooled (Hartung–Knapp) | k=3 | — | 1.1875 (0.9616–1.4665), p=0.073, I²=19.2 | ✓ |

Source: `validation_multi/cindex_HR_summary.csv`, `meta_HK/meta_result.csv`.
The non-significant pooled estimate is reported as non-significant in the manuscript.

## C. Module preservation (primary finding) — VERIFIED

Zsummary 15.853 (ACRG) / 16.797 (GSE15459) / 17.083 (GSE84437) — all "STRONG (Z>10)".
External eigengene Cox HR/SD 1.274, 1.548, 1.237 (p=2.2e-3, 6.8e-5, 2.1e-3) — manuscript
states "1.24–1.55, all p<0.005". ✓
Source: `module_preservation/preservation_summary_RED.csv`, `module_eigengene_cox_external.csv`.

## D. Single-cell localisation (non-circular) — VERIFIED
23 rows, **23/23 Fibroblast-dominant**, median fraction 0.960 — matches "all 23/23 … median 0.96".
Source: `scrna/gene_dominant_celltype_noncircular.csv`.

## E. Decision-analytic value — VERIFIED
ΔC +0.004538 (manuscript "+0.005"), LRT p=0.002046, IDI@3y 0.021, cNRI@3y 0.175. ✓
Source: `external_utility_ACRG/added_value_external.csv`.

## F. Microbiome / batch confounding — VERIFIED
RF cancer-vs-control AUC 0.916119 (CI 0.896–0.936); flowcell-predictability 0.775891
(manuscript "78%") vs majority baseline 0.544992 (manuscript "55%"). ✓
Alpha diversity: Observed δ=−0.27 p=3.3e-07 (median 47→30); Shannon −0.122 p=0.0211;
Simpson −0.072 p=0.172 (**correctly reported as non-significant**). ✓
Source: `microbiome_biomarker/05_rf_metrics_and_batch_sanity.csv`, `02_alpha_effectsizes_cascade.csv`.

## G. Mendelian randomisation — VERIFIED
All six IVW estimates null; smallest p = 0.3483 (Streptococcus) — manuscript states 0.35.
No exposure reaches p<0.05. ✓  Source: `mr_real/MR_results_all_methods_REAL.csv`.

## H. GSEA — VERIFIED (source-file provenance clarified)
Manuscript §3.1 values come from the **integrated TCGA+GTEx** analysis, as the Methods state:
Manuscript-quoted values, all matching: E2F 3.751→"3.75", G2M 3.642→"3.64",
OXPHOS −2.542→"−2.54"; EMT-diffuse 3.174 padj 1.58e-40 →"3.17, 1.6e-40". ✓
Additionally verified in-file but **not numerically quoted** in the manuscript (named
qualitatively only): MYC targets NES 2.669; fatty-acid metabolism NES −2.579.
Source: `enrichment_integrated/fgsea_Hallmark_integrated.csv`,
`enrichment/GSEA_Hallmark_DiffuseVsIntestinal.csv`.
Note: the TCGA-only file (`enrichment/GSEA_Hallmark_TumorVsNormal.csv`) carries different
values (E2F 3.642, OXPHOS −3.102) because it is a different contrast — not a discrepancy.

## I. Immune deconvolution — VERIFIED
T cells (MCP) vs measured leukocyte % : ρ=0.665575, p=3.57e-36 (n=272) — manuscript "ρ 0.67".
CD8 vs leukocyte ρ=0.468413 (verified in file; not quoted numerically in the manuscript).
The source file also records CD8 vs lymphocyte-infiltration ρ=0.0198 (p=0.74) — a near-zero
correlation that was **absent from the manuscript when this audit began and has since been
added** to §3.3 (commit f1e8f26); see "Coverage gap" below. Source: `immune/validation_vs_measured.csv`.

## J. DepMap dependency — VERIFIED
PIK3CA −0.7420 (frac 0.571, dependent+**selective**), MTOR −1.1841 (frac 1.000),
CDK4 −0.8248 (0.657), CDK6 −0.5485 (0.514); FGFR1–4 all **not** dependent. ✓
Source: `depmap/gastric_dependency.csv`.

## K. Prior-signature benchmark — VERIFIED
n=383, 156 events, C_apparent 0.5446, **C_optimism_corrected 0.4807**, KM log-rank p=0.0039. ✓
Source: `base_paper_replication/model_performance.csv`.

## L. Signature stability — VERIFIED
B=200 bootstraps, 25/25 signature genes in candidate pool, **13 selected in >50%**. ✓
Source: `signature_stability/stability_summary.csv`.

## M. Anti-fabrication checks — PASSED

1. **Upstream data is real and predates the manuscript.**
   `data/processed/TCGA_STAD_processed.RData` = 235,007,591 bytes, mtime 21 May 2026 —
   two months before manuscript preparation. Not synthesised at write-time.
2. **Independent recompute from raw expression confirms the data are real and internally consistent.**
   Loading the RData directly: 18,419 genes × 448 samples; barcode-derived
   **412 tumours / 36 normals** — this sample split *is* stated in the manuscript and matches.
   Recomputed from scratch (these per-gene statistics are **not** quoted in the manuscript;
   they are an independent check that the matrix contains real biology, not synthetic values):
   MKI67 tumour 7.668 vs normal 5.071, Wilcoxon p=1.32e-15; COL1A1 11.379 vs 8.611,
   p=6.47e-18; SERTM1 −5.207 vs −0.288, p=3.12e-18. All three move in the biologically
   expected direction (proliferation and collagen up in tumour), which a fabricated or
   randomised matrix would not reproduce.
3. **Negative results are reported as negative** — GSE84437 non-validation, non-significant
   meta-analysis, Simpson diversity ns, all-null MR, FGFR non-dependency, prior signature
   at chance. A fabricated manuscript does not carry this many disclosed failures.
4. **No rounding in a favourable direction detected** in any checked value.

## Coverage gap identified during verification

This audit checked that every number the manuscript **states** is correct. It also surfaced
one result that was present in the source files but **not** carried into the manuscript: the
deconvolution CD8 estimate correlates with measured leukocyte percentage (ρ=0.468) but
**not** with measured lymphocyte-infiltration percentage (ρ=0.0198, p=0.74).

**Status: CLOSED.** This limitation has been added to §3.3 of the manuscript (commit
f1e8f26), which now states that the deconvolution captures overall immune burden rather
than lymphocyte subset composition. It was a *reporting completeness* point, not an error
or a falsification — no stated claim depended on it.

## Conclusion

No falsification, no fabrication, and no unsupported headline claim was found. Every
verified number traces to a named source file, and the raw data independently reproduces
the reported statistics. The one apparent GSEA mismatch resolved to correct
integrated-versus-TCGA-only file provenance, consistent with the Methods.

**Reviewer risk is scientific, not integrity-related:** the modest discrimination
(C≈0.61), the non-significant pooled HR, and the confirmatory nature of the CAF finding
are genuine properties of the data and are disclosed in the manuscript. Those may attract
methodological critique — but they are not defects of honesty, and disclosing them is what
makes the paper defensible.

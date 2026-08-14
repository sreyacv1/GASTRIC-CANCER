# Tables and Supplementary Tables

*Gastric-cancer multi-omics study. Main Tables 1-2 summarize results reported in the manuscript text. Each supplementary table below reproduces the full contents of its source file, with the source path given so every value can be traced to the script that produced it.*

## Table 1 — External validation of the 25-gene signature across four cohorts

| Cohort | Platform | N (events) | C-index | Log-rank p | HR high-vs-low (95% CI) |
|---|---|---|---|---|---|
| TCGA-STAD (training) | RNA-seq | 383 (156) | 0.72 | 2.0×10⁻¹² | — |
| GSE62254 / ACRG | Affymetrix | 300 (152) | 0.61 | 8.6×10⁻⁵ | 1.90 (1.37–2.62) |
| GSE15459 | Affymetrix | 191 (95) | 0.58 | 0.014 | 1.68 (1.11–2.54) |
| GSE84437 | Illumina | 431 (207) | 0.53 | 0.46 | 1.11 (0.84–1.46) |

C-index is Harrell's concordance. The TCGA row is the training cohort and
its C-index is apparent, not leakage-corrected; the honest nested-CV estimate
is C = 0.611 (95% CI 0.562-0.659). HR is high-versus-low risk at the median split.

## Table 2 — Two-sample Mendelian-randomization estimates for six microbial exposures

| Exposure | Axis | nSNP | IVW OR (95% CI) | p |
|---|---|---|---|---|
| Anti-*H. pylori* IgG seropositivity | *H. pylori* | 17 | 0.96 (0.71–1.30) | 0.79 |
| *Streptococcus* (genus) | oral-origin | 15 | 1.10 (0.90–1.34) | 0.35 |
| *Fusobacterium* A | oral-origin | 23 | 1.04 (0.79–1.36) | 0.79 |
| *Prevotella 9* | oral-origin | 15 | 0.98 (0.84–1.14) | 0.76 |
| *Veillonella* | oral-origin | 8 | 1.04 (0.84–1.29) | 0.69 |
| *Lactobacillus* | protective-cand. | 10 | 0.96 (0.84–1.09) | 0.51 |

Inverse-variance-weighted (IVW) odds ratios per unit exposure, ancestry-matched
European outcome (ebi-a-GCST90018849; 1,029 cases). nSNP is the number of instruments retained after harmonization and clumping.
No exposure yielded three or more genome-wide-significant SNPs, so all six use a
locus-wide-suggestive threshold of p<1x10-5 (recorded per exposure in
Supplementary Table S1); this can bias toward the null through weak-instrument
effects, though the per-SNP F-statistics (min F 19.3) indicate the retained
instruments are not weak.
No exposure reaches significance; per-exposure instrument strength and
MR-PRESSO global tests are in Supplementary Tables S1 and S2.

## Supplementary Table S1 — Per-exposure Mendelian-randomization instruments

Instrument-selection threshold, single-nucleotide-polymorphism counts before clumping and after harmonization, per-SNP F statistics, inverse-variance-weighted odds ratio with confidence interval and P value, and the MR-PRESSO global test P value, for all six microbial exposures.

| Exposure | P threshold | SNPs pre-clump | SNPs used | Mean F | Min F | IVW OR | 95% CI | IVW P | MR-PRESSO global P |
|---|---|---|---|---|---|---|---|---|---|
| H. pylori IgG seropositivity | 5e-8 -> 1e-5 (locus-wide) | 21 | 17 | 21.6 | 19.7 | 0.96 | 0.71-1.30 | 0.79 | 0.052 |
| Streptococcus (genus) | 5e-8 -> 1e-5 (locus-wide) | 17 | 15 | 22.5 | 19.3 | 1.1 | 0.90-1.34 | 0.35 | 0.159 |
| Fusobacterium | 5e-8 -> 1e-5 (locus-wide) | 25 | 23 | 22.8 | 19.6 | 1.04 | 0.79-1.36 | 0.79 | 0.611 |
| Prevotella | 5e-8 -> 1e-5 (locus-wide) | 19 | 15 | 20.9 | 19.4 | 0.98 | 0.84-1.14 | 0.76 | 0.605 |
| Veillonella | 5e-8 -> 1e-5 (locus-wide) | 11 | 8 | 21.2 | 19.9 | 1.04 | 0.84-1.29 | 0.69 | 0.676 |
| Lactobacillus | 5e-8 -> 1e-5 (locus-wide) | 12 | 10 | 22.1 | 20.3 | 0.96 | 0.84-1.09 | 0.51 | 0.542 |

All instruments exceed the conventional F > 10 weak-instrument threshold (minimum F 19.3). No exposure reaches nominal significance.

Source: `results/mr_real/MR_per_exposure_instruments_REAL.csv` (6 exposures).

## Supplementary Table S2 — MR-PRESSO global tests

Residual sum of squares and global-test P value per exposure. The global test detects horizontal pleiotropy across the instrument set; no exposure reaches P < 0.05, so no outlier correction was applied.

| Exposure | RSS observed | Global-test P |
|---|---|---|
| H. pylori IgG seropositivity | 32.601 | 0.052 |
| Streptococcus (genus) | 23.073 | 0.159 |
| Lactobacillus | 10.286 | 0.542 |
| Prevotella | 13.969 | 0.605 |
| Fusobacterium | 21.781 | 0.611 |
| Veillonella | 6.394 | 0.676 |

Source: `results/mr_real/MR_PRESSO_global_REAL.csv` (6 exposures).

## Supplementary Table S3 — TRIPOD checklist

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

Source: `TRIPOD_checklist.md`.

## Supplementary Table S4 — Alpha-diversity effect sizes

Median values and Cliff's delta for the non-ulcer-control versus cancer-adjacent-mucosa (GCN) contrast. The non-ulcer, ulcer and GCN groups share sequencing flowcells, so contrasts among them are not confounded by batch (unlike any contrast involving tumor tissue, whose libraries sit on separate flowcells). Negative delta indicates lower diversity in GCN.

| Metric | Median (non-ulcer control) | Median (cancer-adjacent, GCN) | Cliff's delta | Wilcoxon P |
|---|---|---|---|---|
| Observed | 47 | 30 | -0.27 | 3.3e-07 |
| Shannon | 2.176 | 1.923 | -0.122 | 0.0211 |
| Simpson | 0.781 | 0.756 | -0.072 | 0.172 |

Richness (Observed) and Shannon diversity decline significantly, whereas the evenness-weighted Simpson index does not (P = 0.172) — the decline is driven by loss of rare taxa rather than by a shift in community evenness.

Source: `results/microbiome_biomarker/02_alpha_effectsizes_cascade.csv`.

## Supplementary Table S5 — Non-circular single-cell localization

Dominant cell type and the fraction of total expression in that type, for the 23 hub genes that were **not** used to annotate the fibroblast cluster. Excluding the annotation markers removes the circularity in assigning the module to fibroblasts.

| Gene | Dominant cell type | Fraction in dominant type | Panel |
|---|---|---|---|
| ANTXR1 | Fibroblast | 0.863 | stromal_hub |
| COL6A3 | Fibroblast | 0.977 | stromal_hub |
| SPARC | Fibroblast | 0.653 | stromal_hub |
| POSTN | Fibroblast | 0.973 | stromal_hub |
| MXRA8 | Fibroblast | 0.98 | stromal_hub |
| FN1 | Fibroblast | 0.914 | stromal_hub |
| COL8A1 | Fibroblast | 0.998 | stromal_hub |
| VCAN | Fibroblast | 0.972 | stromal_hub |
| GFPT2 | Fibroblast | 0.527 | stromal_hub |
| THBS2 | Fibroblast | 0.82 | stromal_hub |
| AEBP1 | Fibroblast | 0.924 | stromal_hub |
| BGN | Fibroblast | 0.972 | stromal_hub |
| SULF1 | Fibroblast | 0.773 | stromal_hub |
| COL5A1 | Fibroblast | 0.937 | stromal_hub |
| BICC1 | Fibroblast | 0.97 | stromal_hub |
| FBN1 | Fibroblast | 0.932 | stromal_hub |
| ISLR | Fibroblast | 0.902 | stromal_hub |
| CDH11 | Fibroblast | 0.96 | stromal_hub |
| SFRP4 | Fibroblast | 0.987 | stromal_hub |
| SCARF2 | Fibroblast | 0.736 | stromal_hub |
| PRRX1 | Fibroblast | 0.979 | stromal_hub |
| FNDC1 | Fibroblast | 0.983 | stromal_hub |
| ITGA11 | Fibroblast | 0.984 | stromal_hub |

All 23 genes remain fibroblast-dominant.

Source: `results/scrna/gene_dominant_celltype_noncircular.csv`.

## Supplementary Table S6 — GSE84437 pT-stage-stratified validation

Discrimination of the 25-gene signature within pT-stage strata of GSE84437, the cohort in which the signature failed. Stratification does not recover discrimination: the C-index remains below 0.5 in every stratum.

| Subset | n | Events | Harrell C | C standard error | HR per SD | 95% CI | P | C below 0.5 |
|---|---|---|---|---|---|---|---|---|
| All (T1-T4) | 431 | 207 | 0.465 | 0.021 | 1.14 | 1.00-1.29 | 0.0512 | TRUE |
| Early (T1-T3) | 141 | 47 | 0.444 | 0.046 | 1.21 | 0.96-1.54 | 0.108 | TRUE |
| Late (T4) | 290 | 160 | 0.483 | 0.024 | 1.09 | 0.93-1.27 | 0.29 | TRUE |
| T2-T3 only | 130 | 45 | 0.442 | 0.047 | 1.23 | 0.96-1.56 | 0.0978 | TRUE |

Reported as a tested negative result.

Source: `results/validation_multi/GSE84437_Tstage_stratified.csv`.

## Supplementary Table S7 — WGCNA soft-thresholding power sensitivity

Network construction repeated across four candidate soft-thresholding powers. Supports Figure 6(c).

| Soft-thresholding power | Scale-free topology R2 | Mean connectivity | Modules recovered | CAF module | Hubs in CAF module (of 12) | HR per SD | HR lower | HR upper | Cox P | Hub-8 co-membership |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 0.877 | 90.1 | 11 | red | 12 | 1.307 | 1.116 | 1.532 | 0.00093 | 1 |
| 6 | 0.878 | 16.9 | 9 | green | 11 | 1.308 | 1.113 | 1.537 | 0.00114 | 0.875 |
| 9 | 0.879 | 5.4 | 7 | yellow | 11 | 1.293 | 1.095 | 1.526 | 0.00239 | 0.875 |
| 12 | 0.873 | 2.2 | 6 | yellow | 12 | 1.283 | 1.096 | 1.501 | 0.0019 | 1 |

The fibroblast/CAF module is recovered at every power, and its eigengene hazard ratio is stable (1.283-1.308 per standard deviation, Cox P 0.0009-0.0024), so the prognostic result does not depend on the power chosen.

Source: `results/wgcna_real/power_robustness_summary.csv`.

---

All other result tables and figures are provided in the project repository under `results/`, with each mapped to its generating script in `PIPELINE.md`.

## Supplementary Table S8 — Data-acquisition and reproducibility checklist

*Gastric-cancer multi-omics study. Every dataset used, its accession, access route,
version/build, the raw→processed entry point, and the pipeline script that consumes it.
Datasets are listed in the order the pipeline uses them. This lets a reader reacquire
every input from its primary repository and re-enter the pipeline at the correct stage.*

### A. Host transcriptome (differential expression, WGCNA, signature, immune, drug-repurposing)

- [ ] **TCGA-STAD** — primary tumor + solid-tissue normal RNA-seq.
  - Accession/portal: GDC / UCSC Xena (project TCGA-STAD).
  - Build: STAR gene counts, GENCODE v36; 448 samples (412 tumor, 36 normal).
  - Clinical: Lauren class, AJCC stage, grade, TCGA molecular subtype, driver mutations,
    TMB, measured leukocyte/lymphocyte %, overall survival.
  - Raw→processed entry point: count matrix → `results/rdata/tcga_processed.RData`.
  - Consuming scripts: DEG (`analysis/…DEG…`), WGCNA (`analysis/…wgcna…`),
    immune deconvolution (`08_immune_deconvolution.R`), enrichment,
    microbe-response GSEA (`31_microbe_response_enrichment.R`).
- [ ] **GTEx v10** — normal stomach reference (DE denominator).
  - Portal: gtexportal.org. Combined with TCGA normals for the integrated tumor-vs-normal contrast.

### B. External survival / DE validation cohorts (transcriptome)

- [ ] **GSE62254 (ACRG)** — 300 patients; Affymetrix. Primary external validation of the signature.
  - Access: GEO. Local: `data/geo/GSE62254.rda`. Consumer: `12_multicohort_validation.R`, `32_nested_cv_signature.R`.
- [ ] **GSE15459** — 191 patients. Access: GEO. Local: `data/geo/GSE15459_es.rds` + `GSE15459_outcome.xls`.
- [ ] **GSE84437** — 431 patients. Access: GEO. Local: `data/geo/GSE84437_es.rds`. (Reported negative — retained transparently.)
- [ ] **GSE27342, GSE63089** — independent DE-concordance cohorts. Access: GEO.
  Consumer: DEG cross-cohort concordance (`results/tables/DEG_GEO_*.csv`).

### C. Single-cell transcriptome (module localization)

- [ ] **GSE134520** — premalignant-to-early gastric-cancer scRNA-seq.
  - Access: GEO. Processing: Seurat v4 (QC 200–6,000 genes, MT<20% → 43,992 cells;
    LogNormalize; 2,000 HVGs; PCA 30; SNN res 0.5; UMAP). 8 cell types by canonical markers.
  - Consumer: scRNA localization → `results/scrna/`.

### D. Tissue microbiome (16S)

- [ ] **PRJDB20660 (Japan, discovery)** — 944 libraries, 16S V3–V4.
  - Access: DDBJ, **raw FASTQ**. Processing: DADA2 from raw → ASV table
    (897 samples × 314 genera after host/off-target removal).
  - Read-tracking validation: per-sample counts vs published Supplementary Table 3, Pearson r=0.98
    (`results/microbiome_biomarker/00_readtracking_concordance.csv`).
  - Consumer: `23_dada2_16S.R` → `24_microbiome_real.R`; reference DB `silva_nr99_v138.1`.
- [ ] **PRJNA641258 (Italy, validation)** — paired tumor/non-tumor, V3–V4. Access: NCBI/ENA. (Ravegnini 2020, PMC7766162.)
- [ ] **PRJNA413125 (Portugal, validation)** — gastritis vs carcinoma, V5–V6. Access: NCBI/ENA. (Ferreira 2018.)
  - **Batch caveat:** the discovery tumor/normal split is confounded with sequencing flowcell (§3.8);
    validation cohorts are the batch-clean test.

### E. GWAS summary statistics (Mendelian randomization)

*All accessions independently verified against GWAS-Catalog records; raw `gwasinfo()`/lookup output
archived at `results/mr_real/gwas_catalog_verification.json`.*

- [ ] **Exposures (IEU OpenGWAS):**
  - `ebi-a-GCST90006910` — anti-*H. pylori* IgG seropositivity (Butler-Laporte 2020; European, PMID 33204752).
  - `ebi-a-GCST90017070` — genus *Streptococcus* (id.1853; MiBioGen/Kurilshikov 2021; European, PMID 33462485).
  - `ebi-a-GCST90017045` — genus *Prevotella 9* (id.11183; MiBioGen; European, PMID 33462485).
  - `ebi-a-GCST90017088` — genus *Veillonella* (id.2198; MiBioGen; European, PMID 33462485).
  - `ebi-a-GCST90017030` — genus *Lactobacillus* (id.1837; MiBioGen; European, PMID 33462485).
  - `ebi-a-GCST90032406` — *Fusobacterium* A abundance in stool (Qin 2022; European/Finland, PMID 35115689).
- [ ] **Outcomes:**
  - `ebi-a-GCST90018849` — gastric cancer, European, 1,029 cases (Sakaue 2021, PMID 34594039). **Primary.**
  - `ebi-a-GCST90018629` — gastric cancer, East-Asian, 7,921 cases (Sakaue 2021). Cross-ancestry power sensitivity only.
  - Access token: OpenGWAS JWT via `OPENGWAS_JWT` env var. Consumer: `11_real_mr.R`.
  - Harmonised per-SNP instrument tables archived (`results/mr_real/`) so exact SNP sets survive upstream reprocessing.

### F. Environment & reproducibility

- [ ] R 4.3.3; pinned package versions in `package_versions.csv`; full `sessionInfo()` with the code.
- [ ] Complete-case analysis throughout; no imputation.
- [ ] Bundled `r_env/` has a build-time hardcoded path — after relocation use the `R_HOME`-override
      invocation or a system R (documented in `PIPELINE.md`).
- [ ] Every `results/` output maps to its generating script in `PIPELINE.md`;
      figure→source mapping in `FIGURE_SOURCES.md`.

### G. Analyses NOT performed (stated to prevent over-reading)

- [ ] **No genuinely paired per-patient tissue-microbiome + host-transcriptome cohort exists** in this
      study (16S and transcriptome are from different, unpaired cohorts). A per-genus host-gene MaAsLin2
      association was run in a superseded pipeline on a bridged, non-verifiable sample set; its outputs are
      not reproducible from retained files, so it is excluded and no host-gene↔genus claim is made. The
      microbe-response result (§3.9) is a pathway-level GSEA on the TCGA tumor-vs-normal ranking only —
      hypothesis-generating context, not a host–microbiome axis.

Source: `DATA_ACQUISITION_CHECKLIST.md`.

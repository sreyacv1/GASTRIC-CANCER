# Tables

| Exposure | p_threshold | nSNP_prclump | nSNP_used | mean_F | min_F | IVW_OR | IVW_CI | IVW_p | PRESSO_global_p |
|---|---|---|---|---|---|---|---|---|---|
| H. pylori IgG seropositivity | 5e-8 -> 1e-5 (locus-wide) | 21 | 17 | 21.6 | 19.7 | 0.96 | 0.71-1.30 | 0.79 | 0.052 |
| Streptococcus (genus) | 5e-8 -> 1e-5 (locus-wide) | 17 | 15 | 22.5 | 19.3 | 1.1 | 0.90-1.34 | 0.35 | 0.159 |
| Fusobacterium | 5e-8 -> 1e-5 (locus-wide) | 25 | 23 | 22.8 | 19.6 | 1.04 | 0.79-1.36 | 0.79 | 0.611 |
| Prevotella | 5e-8 -> 1e-5 (locus-wide) | 19 | 15 | 20.9 | 19.4 | 0.98 | 0.84-1.14 | 0.76 | 0.605 |
| Veillonella | 5e-8 -> 1e-5 (locus-wide) | 11 | 8 | 21.2 | 19.9 | 1.04 | 0.84-1.29 | 0.69 | 0.676 |
| Lactobacillus | 5e-8 -> 1e-5 (locus-wide) | 12 | 10 | 22.1 | 20.3 | 0.96 | 0.84-1.09 | 0.51 | 0.542 |

Table 1: Per-exposure Mendelian-randomization instrument summary. Threshold, pre-clump and harmonized SNP counts with mean and minimum F-statistics for each exposure. Source: `results/mr_real/MR_per_exposure_instruments_REAL.csv`.

| exposure | RSSobs | PRESSO_global_p |
|---|---|---|
| H. pylori IgG seropositivity | 32.601 | 0.052 |
| Streptococcus (genus) | 23.073 | 0.159 |
| Lactobacillus | 10.286 | 0.542 |
| Prevotella | 13.969 | 0.605 |
| Fusobacterium | 21.781 | 0.611 |
| Veillonella | 6.394 | 0.676 |

Table 2: MR-PRESSO global tests. Observed residual sum of squares and global-test P value per exposure. Source: `results/mr_real/MR_PRESSO_global_REAL.csv`.

| metric | median_Nonul | median_GCN | cliffs_delta | wilcox_p |
|---|---|---|---|---|
| Observed | 47 | 30 | -0.27 | 3.3e-07 |
| Shannon | 2.176 | 1.923 | -0.122 | 0.0211 |
| Simpson | 0.781 | 0.756 | -0.072 | 0.172 |

Table 3: Alpha-diversity effect sizes along the gastritis-to-carcinoma cascade. Median differences and Cliff's delta for Observed richness, Shannon and Simpson indices. Source: `results/microbiome_biomarker/02_alpha_effectsizes_cascade.csv`.


# Supplementary Tables

- **Supplementary Table S1 — Per-exposure MR instruments.** Threshold, pre-clump and harmonized-used SNP counts, mean/min per-SNP F, IVW OR/CI/p, and MR-PRESSO global p for all six exposures. File: `results/mr_real/MR_per_exposure_instruments_REAL.csv`.
- **Supplementary Table S2 — MR-PRESSO global tests.** RSSobs and global-test p per exposure. File: `results/mr_real/MR_PRESSO_global_REAL.csv`.
- **Supplementary Table S3 — TRIPOD checklist.** Item-by-item reporting map for the prognostic model. File: `TRIPOD_checklist.md`.
- **Supplementary Table S4 — Alpha-diversity effect sizes.** Median differences and Cliff's δ (Observed / Shannon / Simpson) across the non-confounded gastritis-to-cancer cascade. File: `results/microbiome_biomarker/02_alpha_effectsizes_cascade.csv`.
- **Supplementary Table S5 — Non-circular single-cell localization.** Dominant cell type and fraction for the 23 hub genes not used to annotate the fibroblast cluster. File: `results/scrna/gene_dominant_celltype_noncircular.csv`.
- **Supplementary Table S6 — GSE84437 pT-stage-stratified validation.** Harrell C-index and per-SD hazard ratios for the 25-gene signature within pT-stage strata (all / early pT1–T3 / pT4 / pT2–T3), showing discrimination is not recovered by stratification (C<0.5 throughout). File: `results/validation_multi/GSE84437_Tstage_stratified.csv`.
- **Supplementary Table S7 — Data-acquisition & reproducibility checklist.** Every dataset used, its accession, access route, version/build, raw→processed entry point, and consuming pipeline script, ordered by pipeline stage; includes analyses explicitly not performed. File: `DATA_ACQUISITION_CHECKLIST.md`.
- **Supplementary Table S8 — WGCNA soft-thresholding power sensitivity.** For each candidate soft-thresholding power (3, 6, 9, 12): scale-free topology fit, mean connectivity, number of modules recovered, the module carrying the fibroblast/CAF hub genes, its eigengene hazard ratio per standard deviation with 95% confidence interval and Cox P value, and the fraction of the eight hub genes remaining co-clustered. Supports Figure 6(c). File: `results/wgcna_real/power_robustness_summary.csv`.
All other result tables and figures are provided in the project repository under `results/`, with each mapped to its generating script in `PIPELINE.md`.

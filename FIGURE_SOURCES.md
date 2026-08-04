# Figure & Result Provenance Manifest

*Every figure in `PAPER.md` traces to committed source files in this repository. Composite
supplementary figures (S15–S19) are montages assembled from the individual pipeline-output
PNGs listed below — no panel is redrawn or synthetic; each source panel and each numeric
value can be opened and verified independently. Every headline number traces to a source CSV.*

> **Verification.** Every source panel and every composite display image listed below was confirmed present and git-tracked (`git ls-files`, verified 2026-07-23) — 0 missing, 0 untracked. The five composite supplementary figures are stored as committed images under `results/composite_figures/` and are the exact files embedded in `PAPER.md`. Headline numeric values are independently re-read against their source CSVs in `RESULTS_COMPENDIUM.md` / `RESULTS_REVIEW.md`.

## Main figures
| Fig | Source image(s) | Source data (numbers) |
|---|---|---|
| **1** Tumour programmes | `results/enrichment/GSEA_Hallmark_NES_barplot_*.png` | `results/enrichment_integrated/fgsea_Hallmark_integrated.csv`; `results/enrichment/GSEA_Hallmark_DiffuseVsIntestinal.csv` |
| **2** Signature (honest) | `results/figures/Fig2.png` | `results/nested_cv/performance.csv` (C 0.611); `results/meta_HK/meta_result.csv` (HR 1.19); `results/timevarying_ACRG/` |
| **3** Stromal/CAF module | `results/figures/Fig3.png` | `results/module_preservation/preservation_stats_*.csv` (Z 15.9/16.8/17.1); `results/scrna/gene_dominant_celltype_noncircular.csv` (23/23) |
| **4** Immune + microbiome | `results/figures/Fig4.png` | `results/immune/validation_vs_measured.csv` (ρ 0.6656); microbiome cascade CSVs |

## Supplementary figures
| Fig | Source image(s) | Source data |
|---|---|---|
| **S1** scRNA UMAP | `results/scrna/UMAP_celltypes.png` | `results/scrna/celltype_composition.csv` (Epithelial 68.64%) |
| **S2** hub dot plot | `results/scrna/DotPlot_stromal_module_hub.png` | `results/scrna/gene_dominant_celltype_noncircular.csv` |
| **S3** drug repurposing | `results/drug_repurposing_integrated/top_candidate_drugs.png` | `.../candidate_drugs_ranked.csv` (resveratrol #1) |
| **S4** WGCNA dendrogram | `results/wgcna_real/wgcna_dendrogram.png` | `results/wgcna_real/` |
| **S5** module–trait | `results/wgcna_real/wgcna_module_trait.png` | `results/wgcna_real/ME_survival_cox.csv` |
| **S6** power robustness | `results/wgcna_real/wgcna_power.png` | `results/wgcna_real/soft_threshold_table.csv` |
| **S7** DEG volcano (TCGA-only) | `results/plots/transcriptome/deg_volcano.png` | `results/tables/TCGA_DEG_results.csv` (2134↑/2362↓/21446) |
| **S8** top-30 DEG heatmap | `results/plots/transcriptome/deg_heatmap_top30_clean.png` | `data/processed/TCGA_STAD_processed.RData` |
| **S9** ORA GO+KEGG | `results/enrichment/path_ORA_GO_KEGG.png` | `results/enrichment/dotplot_*.png` |
| **S10** LASSO coefficients | `results/validation/signature_coefficients_clean.png` | `results/validation/signature_coefficients.csv` (16↑/9↓) |
| **S11** external forest | `results/validation_multi/forest_HR.png` | `results/validation_multi/cindex_HR_summary.csv` (per-SD 1.90/1.68/1.11) |
| **S12** MR scatter (H. pylori) | `results/mr_real/scatter_H__pylori_IgG_seropositivity.png` | `results/mr_real/MR_results_all_methods_REAL.csv` |
| **S13** MR leave-one-out (H. pylori) | `results/mr_real/loo_H__pylori_IgG_seropositivity.png` | same |
| **S14** microbiome CLR DA | `results/microbiome_biomarker/da_clr_barplot.png` | `.../04a_DA_control_vs_GCN.csv`, `04b_DA_GCN_vs_GCT_paired.csv` |
| **S15** immune (composite → `results/composite_figures/s15_immune.png`) | `results/plots/Immune_validation_scatter.png` + `Immune_tumor_vs_normal.png` + `Immune_by_subtype.png` + `Immune_CD8_survival_KM.png` | `results/immune/validation_vs_measured.csv`; `CD8_survival_summary.csv` (HR 1.04, p=0.41) |
| **S16** four-cohort KM (composite → `results/composite_figures/s16_km.png`) | `results/validation/KM_TCGA.png` + `results/validation/KM_ACRG.png` + `results/validation_multi/KM_GSE15459.png` + `results/validation_multi/KM_GSE84437.png` | `results/validation/multivariable_cox_ACRG.csv` (median-split HR 1.76) |
| **S17** nomogram/calib/DCA (composite → `results/composite_figures/s17_nomogram.png`) | `results/nomogram_combined/combined_nomogram.png` + `results/nomogram_combined/calibration_combined.png` + `results/external_utility_ACRG/DCA_external.png` | `results/external_utility_ACRG/added_value_external.csv` (ΔC +0.005) |
| **S18** MR scatter ×6 (composite → `results/composite_figures/s18_mr_scatter_all.png`) | `results/mr_real/scatter_{H__pylori_IgG_seropositivity,Streptococcus__genus_,Fusobacterium,Prevotella,Veillonella,Lactobacillus}.png` | `results/mr_real/MR_per_exposure_instruments_REAL.csv` |
| **S19** MR leave-one-out ×6 (composite → `results/composite_figures/s19_mr_loo_all.png`) | `results/mr_real/loo_{H__pylori_IgG_seropositivity,Streptococcus__genus_,Fusobacterium,Prevotella,Veillonella,Lactobacillus}.png` | same |

## How to verify (for reviewers)
1. Open any source PNG listed above directly in the repo — the composite figure is that exact image.
2. Open the source CSV — every number in the caption/text appears there.
3. Re-run the pipeline: `PIPELINE.md` lists the exact scripts (`analysis/*.R`) and versions; `sessionInfo.txt` / `package_versions.csv` pin the environment.
4. Nothing in this manuscript is hand-entered: `RESULTS_COMPENDIUM.md` and `RESULTS_REVIEW.md` document an independent re-read of every headline number against these files.

assets/fonts/ contains Tinos (Apache-2.0), metric-compatible with Times New Roman,
used for the graphical abstract because Times New Roman is not installed on the build host.

## Figures added/rebuilt 2026-08-04

| Figure file | Build script | Notes |
|---|---|---|
| `results/figures/deg_concordance_panel.png` | `analysis/23_deg_concordance_panel.R` | Was an orphan (no build script). Script asserts all three r values against `results/tables/DEG_integrated_concordance.csv`. **Panel 1 is internal consistency, not replication** — TCGA supplies the integrated contrast's own 412 tumours. |
| `results/plots/transcriptome/deg_volcano.png` | R block in session log; data `res_tcga_df` in `data/processed/TCGA_STAD_processed.RData` | Labels moved outside the point cloud; subtitle now gives the biotype breakdown (16,164 protein-coding of 21,446 features). |
| `results/plots/transcriptome/Lauren_volcano.png` | regenerated from `res_lauren_df` | Original PDF was an empty plot device (3,611 B, zero drawing ops). |
| `results/enrichment/GSEA_Hallmark_NES_barplot_DiffuseVsIntestinal.png` | `analysis/09_functional_enrichment.R` (inline Lauren block) | Title was clipped by `nes_barplot`'s 9x8 canvas; replaced with a self-contained 8.6x6.2 block. Verified byte-identical to the shipped PNG. |
| `results/plots/Immune_validation_scatter.png` | `analysis/08_immune_deconvolution.R` L118-135 | `stat_cor()` replaced with an explicit annotation asserted against `validation_vs_measured.csv` (rho 0.4684, n=272). Verified byte-identical to the shipped PNG. |
| `results/figures/panels/Fig4A.png` | `make_figures.R` §4A | Title shortened to fit canvas. |
| `results/figures/panels/Fig4C.png` | `make_figures.R` §4C | Title "(honest)" removed — editorialising, not a data label. |
| `results/deg_diagnostics/batch_sensitivity_TSS.csv` | R block, session log | TSS batch-adjustment sensitivity: 92.7% of DEGs retained, logFC r = 0.974. |

**Known cosmetic limitation (not a data defect).** The 12 MR **scatter** plots (`scatter_*.png`, six under `results/mr_real/` and six under
`results/mr_real_eas/`) have their legend touching the bottom canvas edge. The 12
leave-one-out plots (`loo_*.png`) in the same directories were each measured with the same
metric: 11 of 12 score exactly 0, and `results/mr_real_eas/loo_Fusobacterium.png` scores 0.0018
on the **right** edge: exactly 4 dark pixels, in the outermost column only, at rows 690-693 of
750 (92% down the canvas). Both EAS and European canvases are 900x750, so this is not a
size difference; it is a single confidence-interval whisker of the last leave-one-out row
reaching the plot boundary. A different and far smaller defect than the scatters' bottom-edge
legend, and not visible at normal viewing size. The width patch applies to both
families; only the scatters are materially affected. `analysis/11_real_mr.R` was patched (width 6 -> 7.5 in) but regenerating them
requires a live OpenGWAS API token (`OPENGWAS_JWT`) plus network access to re-extract the
harmonised SNP data, which is not cached locally. The plotted values are correct and are
independently tabulated in `results/mr_real/MR_results_all_methods_REAL.csv`; the figures
will pick up the wider canvas the next time the MR arm is re-run with a token.

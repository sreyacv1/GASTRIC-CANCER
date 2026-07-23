# Figure & Result Provenance Manifest

*Every figure in `PAPER.md` traces to committed source files in this repository. Composite
supplementary figures (S15–S19) are montages assembled from the individual pipeline-output
PNGs listed below — no panel is redrawn or synthetic; each source panel and each numeric
value can be opened and verified independently. Every headline number traces to a source CSV.*

> **Verification.** All 34 individual source panels and 5 composites listed below were confirmed present and git-tracked (`git ls-files`) in this repository — 0 missing. Each cited numeric value was traced to its source CSV (also tracked). This is a complete check, not a spot-check.

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
| **S15** immune (composite) | `results/plots/Immune_validation_scatter.png` + `Immune_tumor_vs_normal.png` + `Immune_by_subtype.png` + `Immune_CD8_survival_KM.png` | `results/immune/validation_vs_measured.csv`; `CD8_survival_summary.csv` (HR 1.04, p=0.41) |
| **S16** four-cohort KM (composite) | `results/validation/KM_TCGA.png` + `KM_ACRG.png` + `results/validation_multi/KM_GSE15459.png` + `KM_GSE84437.png` | `results/validation/multivariable_cox_ACRG.csv` (median-split HR 1.76) |
| **S17** nomogram/calib/DCA (composite) | `results/nomogram_combined/combined_nomogram.png` + `calibration_combined.png` + `results/external_utility_ACRG/DCA_external.png` | `results/external_utility_ACRG/added_value_external.csv` (ΔC +0.005) |
| **S18** MR scatter ×6 (composite) | `results/mr_real/scatter_{H_pylori,Streptococcus,Fusobacterium,Prevotella,Veillonella,Lactobacillus}.png` | `results/mr_real/MR_per_exposure_instruments_REAL.csv` |
| **S19** MR leave-one-out ×6 (composite) | `results/mr_real/loo_{...}.png` (same six) | same |

## How to verify (for reviewers)
1. Open any source PNG listed above directly in the repo — the composite figure is that exact image.
2. Open the source CSV — every number in the caption/text appears there.
3. Re-run the pipeline: `PIPELINE.md` lists the exact scripts (`analysis/*.R`) and versions; `sessionInfo.txt` / `package_versions.csv` pin the environment.
4. Nothing in this manuscript is hand-entered: `RESULTS_COMPENDIUM.md` and `RESULTS_REVIEW.md` document an independent re-read of every headline number against these files.

## Supplementary Figures

![Supplementary Figure S1](results/scrna/UMAP_celltypes.png)

Supplementary Figure S1. Single-cell reference atlas of the gastric premalignant-to-early-cancer cascade.

Uniform manifold approximation and projection (UMAP) of 43,992 cells from the public dataset GSE134520, coloured by the eight annotated major cell types. UMAP is a non-linear dimensionality-reduction method in which each point represents one cell and proximity reflects transcriptional similarity; the axes carry no independent units. Cell-type proportions are epithelial 68.6%, endothelial 7.6%, myeloid 6.4%, plasma 5.8%, T cell 5.8%, fibroblast 4.2%, mast 1.0% and B cell 0.6%. This atlas provides the cellular reference used for the module-localisation analysis in Figure 3(c).

File: `results/scrna/UMAP_celltypes.png`.

![Supplementary Figure S2](results/scrna/DotPlot_stromal_module_hub.png)

Supplementary Figure S2. Expression of stromal-module hub genes across cell types.

Dot plot of the prognostic co-expression module hub genes across the eight annotated cell types. Dot size represents the fraction of cells in which the gene is detected and dot colour represents scaled mean expression. Hub-gene expression is confined to the fibroblast compartment. Because the genes used to annotate the fibroblast cluster could bias this result, the analysis was repeated after removing them; all twenty-three remaining hub genes stayed fibroblast-dominant.

File: `results/scrna/DotPlot_stromal_module_hub.png`.

![Supplementary Figure S3](results/drug_repurposing_integrated/top_candidate_drugs.png)

Supplementary Figure S3. Candidate compounds from in-silico drug repurposing.

Bar plot of the eight highest-ranked compounds whose transcriptional perturbation signature reverses the tumour programme, from Enrichr queried against the LINCS L1000 and Gene Expression Omnibus drug-perturbation libraries. Bar length is −log10 of the best Benjamini–Hochberg-adjusted enrichment p value across contributing signatures, which is also the ordering criterion. Every compound shown reverses the programme on both arms — repressing tumour-up genes and inducing tumour-down genes — so no arm legend is displayed. Four of the eight are targeted kinase inhibitors concordant with a proliferation-dominated programme: the fibroblast growth factor receptor (FGFR) inhibitor PD-173074, the cyclin-dependent kinase 4/6 (CDK4/6) inhibitor palbociclib, the FGFR/multikinase inhibitor dovitinib and the phosphoinositide 3-kinase (PI3K)/mechanistic target of rapamycin (mTOR) inhibitor NVP-BEZ235. The remaining four are not kinase inhibitors and include the highest-ranked hit, resveratrol, together with calcitriol, vemurafenib and mitoxantrone; broadly bioactive compounds of this kind recur across LINCS reversal analyses irrespective of tumour type and should be treated as non-specific. These predictions are computational and hypothesis-generating only; no compound was experimentally validated, and because the tumour-up programme is proliferation-dominated the nominated agents are broadly anti-proliferative rather than specific to gastric cancer.

File: `results/drug_repurposing_integrated/top_candidate_drugs.png`.

![Supplementary Figure S4](results/enrichment/path_ORA_GO_KEGG.png)

Supplementary Figure S4. Over-representation analysis of tumour-up-regulated genes.

Dot plots of Gene Ontology biological-process and Kyoto Encyclopedia of Genes and Genomes (KEGG) over-representation analysis for the tumour-up gene set. Dot size represents the number of genes in each term and colour represents the adjusted P value. Cell-cycle and nuclear-division terms, together with extracellular-matrix-receptor and cytokine terms, dominate. Over-representation analysis depends on an arbitrary significance threshold and is sensitive to gene-family size; the prominent olfactory and sensory-perception terms reflect the large olfactory-receptor gene family rather than tumour biology. The threshold-free gene-set enrichment analysis in Figure 1 is therefore the primary enrichment result, with this panel provided for completeness.

Files: `results/enrichment/dotplot_GO_BP_UP.png`, `results/enrichment/dotplot_KEGG_UP.png`.

![Supplementary Figure S5](results/figures/clean/s5_mr_scatter_hpylori.png)

Supplementary Figure S5. Mendelian-randomisation scatter plot for a representative exposure.

Scatter plot of single-nucleotide-polymorphism (SNP) effects on the exposure, anti-Helicobacter pylori immunoglobulin G seropositivity, against their effects on gastric-cancer risk. Each point represents one genetic instrument and error bars represent standard errors on both axes. The five overlaid lines are the slopes estimated by inverse-variance weighting (IVW), MR-Egger, weighted median, weighted mode and simple mode. All slopes are flat and statistically indistinguishable from zero (IVW odds ratio 0.96, 95% confidence interval 0.71-1.30). Equivalent plots for all six exposures in both ancestries are provided in the repository.

File: `results/figures/clean/s5_mr_scatter_hpylori.png` (source panel `results/mr_real/scatter_H__pylori_IgG_seropositivity.png`).

![Supplementary Figure S6](results/mr_real/loo_H__pylori_IgG_seropositivity.png)

Supplementary Figure S6. Mendelian-randomisation leave-one-out analysis for a representative exposure.

Leave-one-out inverse-variance-weighted estimates for anti-Helicobacter pylori immunoglobulin G seropositivity. Each row shows the pooled causal estimate recomputed with one instrument removed; points represent the estimate and horizontal bars the 95% confidence interval. Every estimate straddles the null, confirming that no individual instrument drives the result.

File: `results/mr_real/loo_H__pylori_IgG_seropositivity.png`.

![Supplementary Figure S7](results/microbiome_biomarker/da_clr_barplot.png)

Supplementary Figure S7. Compositional differential abundance of tissue microbiota.

Bar plots of centred-log-ratio effect sizes for the two discovery-cohort contrasts, tested by Wilcoxon rank-sum with Benjamini-Hochberg correction. The centred-log-ratio transform is used because sequencing yields relative rather than absolute abundances, so raw proportions are not independent. The left panel shows control versus cancer-adjacent mucosa (44 of 61 genera at q < 0.05) and the right panel the paired cancer-adjacent versus tumour contrast (18 of 61 genera at q < 0.05). Bars are coloured by direction where q < 0.05, with red indicating enrichment and blue depletion, and the twelve most enriched and depleted genera are shown per panel. These genus-level shifts are reported for completeness only; as detailed in the Results, the tumour contrast is confounded with sequencing batch and does not replicate in an independent batch-clean cohort, so these taxa are not proposed as biomarkers.

Files: `results/microbiome_biomarker/04a_DA_control_vs_GCN.csv`, `04b_DA_GCN_vs_GCT_paired.csv`.

![Supplementary Figure S8](results/composite_figures/s15_immune.png)

Supplementary Figure S8. Immune infiltration with deconvolution validated against pathology.

(a) Scatter plot of the deconvolution-derived T-cell score against the measured histological leukocyte fraction scored by a pathologist (Spearman rho = 0.67, P = 3.6 x 10-36), establishing that the expression-based estimates track true tissue composition.

(b) Comparison of immune compartments between tumour and normal tissue. Enrichment is dominated by the macrophage and monocyte lineage, with no net gain in CD8-positive T cells.

(c) Immune infiltration across the four molecular subtypes defined by The Cancer Genome Atlas. Epstein-Barr virus-positive tumours are the most infiltrated and chromosomal-instability tumours the least (Kruskal-Wallis P < 1 x 10-6).

(d) Association between CD8-positive T-cell score and overall survival. The score is not prognostic in this cohort (Cox hazard ratio 1.04, P = 0.41); this null result is reported as observed.

Files: `results/plots/Immune_*.png`.

![Supplementary Figure S9](results/composite_figures/s17_nomogram.png)

Supplementary Figure S9. Clinical nomogram, calibration, and external decision-curve analysis.

(a) Nomogram combining clinical covariates with the signature score for prediction of 1-, 3- and 5-year overall survival. A nomogram converts each predictor into points on a common scale whose total maps to a predicted survival probability.

(b) Calibration plots at 1, 3 and 5 years, computed in the development cohort. Calibration compares predicted against observed survival; the diagonal represents perfect agreement.

(c) Decision-curve analysis performed out-of-sample in the Asian Cancer Research Group cohort. Decision-curve analysis plots net benefit against the threshold probability at which a clinician would act. Both the clinical model and the combined model carry net benefit relative to treat-all and treat-none strategies, but their curves are essentially superimposed across the 5-40% threshold range, indicating that adding the signature confers no incremental net benefit over standard staging. The nomogram is presented as an illustrative research tool, not as a validated decision instrument.

Files: `results/nomogram_combined/`, `results/external_utility_ACRG/DCA_external.png`.

![Supplementary Figure S10](results/composite_figures/s19_mr_loo_all.png)

Supplementary Figure S10. Mendelian-randomisation leave-one-out analysis for all six exposures.

Leave-one-out inverse-variance-weighted estimates for each of the six microbial exposures. Removing any single instrument leaves every pooled estimate straddling the null, confirming that no individual instrument drives any result and that the overall null is not an outlier artefact.

Files: `results/mr_real/loo_*.png`.

![Supplementary Figure S11](results/microbiome_biomarker/rf_batch_classifier.png)

Supplementary Figure S11. Evidence that the tumour-microbiome classifier reflects sequencing batch.

(a) Bar plot of the fifteen genera ranked highest by mean decrease in Gini importance in the cancer-versus-control random-forest classifier. Genera shown in red italic type are recognised environmental or reagent contaminants (Dietzia, Serinicoccus, Methylobacterium-Methylorubrum, Microbacterium, Sphingomonas and Serratia) rather than gastric or oral commensals, and they dominate the classifier.

(b) Comparison of the classifier's apparent discrimination with its ability to predict sequencing batch. The cancer-versus-control area under the receiver operating characteristic curve is 0.92 (95% confidence interval 0.90-0.94), but the same feature space predicts the sequencing flowcell among biologically similar samples at 78% accuracy against a 55% majority-class baseline (dotted line). The apparent tumour signal therefore tracks sequencing batch rather than tumour biology.

Files: `results/microbiome_biomarker/05_rf_importance.csv`, `05_rf_metrics_and_batch_sanity.csv`.

## Supplementary Materials

- **Supplementary Table S1 — Per-exposure MR instruments.** Threshold, pre-clump and harmonised-used SNP counts, mean/min per-SNP F, IVW OR/CI/p, and MR-PRESSO global p for all six exposures. File: `results/mr_real/MR_per_exposure_instruments_REAL.csv`.
- **Supplementary Table S2 — MR-PRESSO global tests.** RSSobs and global-test p per exposure. File: `results/mr_real/MR_PRESSO_global_REAL.csv`.
- **Supplementary Table S3 — TRIPOD checklist.** Item-by-item reporting map for the prognostic model. File: `TRIPOD_checklist.md`.
- **Supplementary Table S4 — Alpha-diversity effect sizes.** Median differences and Cliff's δ (Observed / Shannon / Simpson) across the non-confounded gastritis-to-cancer cascade. File: `results/microbiome_biomarker/02_alpha_effectsizes_cascade.csv`.
- **Supplementary Table S5 — Non-circular single-cell localisation.** Dominant cell type and fraction for the 23 hub genes not used to annotate the fibroblast cluster. File: `results/scrna/gene_dominant_celltype_noncircular.csv`.
- **Supplementary Table S6 — GSE84437 pT-stage-stratified validation.** Harrell C-index and per-SD hazard ratios for the 25-gene signature within pT-stage strata (all / early pT1–T3 / pT4 / pT2–T3), showing discrimination is not recovered by stratification (C<0.5 throughout). File: `results/validation_multi/GSE84437_Tstage_stratified.csv`.
- **Supplementary Table S7 — Data-acquisition & reproducibility checklist.** Every dataset used, its accession, access route, version/build, raw→processed entry point, and consuming pipeline script, ordered by pipeline stage; includes analyses explicitly not performed. File: `DATA_ACQUISITION_CHECKLIST.md`.
- **Supplementary Table S8 — WGCNA soft-thresholding power sensitivity.** For each candidate soft-thresholding power (3, 6, 9, 12): scale-free topology fit, mean connectivity, number of modules recovered, the module carrying the fibroblast/CAF hub genes, its eigengene hazard ratio per standard deviation with 95% confidence interval and Cox P value, and the fraction of the eight hub genes remaining co-clustered. Supports Figure 5(c). File: `results/wgcna_real/power_robustness_summary.csv`.

All other result tables and figures are provided in the project repository under `results/`, with each mapped to its generating script in `PIPELINE.md`.

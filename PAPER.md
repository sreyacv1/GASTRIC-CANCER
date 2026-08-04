# A multi-omics analysis of gastric cancer identifies an externally-validated stromal/fibroblast prognostic programme, with cautionary tumour-microbiome and Mendelian-randomisation assessments

**Authors:** ⟦PLACEHOLDER — AUTHOR NAMES + ORCIDs⟧ (e.g. Given-name Family-name¹ [ORCID 0000-0000-0000-0000], …) — Team BBIN6_XX, School of Chemical and Biotechnology, SASTRA Deemed to be University, Thanjavur, India
**Corresponding author:** ⟦PLACEHOLDER — CORRESPONDING AUTHOR name; institutional email; postal address⟧
**Running title:** A stromal/fibroblast prognostic programme in gastric cancer

---

## Abstract

Prognostic stratification of gastric cancer beyond tumour stage remains limited, and claims of a causal tumour microbiome are frequently made without genetic evidence. This study aimed to identify and externally validate the transcriptional biology underlying gastric-cancer prognosis, and to test — rather than assume — whether microbial exposures causally influence disease risk. A multi-omics design was used, spanning bulk and single-cell transcriptomes, tissue 16S ribosomal microbiome sequencing, and germline-genetic causal inference. Differential expression integrated The Cancer Genome Atlas and Genotype-Tissue Expression cohorts (412 tumours versus 443 normal references), and a penalised Cox prognostic signature was trained on The Cancer Genome Atlas and validated in three independent cohorts (GSE62254, GSE15459 and GSE84437; 922 patients). Weighted gene co-expression network analysis with external module-preservation testing, immune deconvolution validated against measured leukocyte fraction, single-cell localisation, and two-sample Mendelian randomisation of anti-*Helicobacter pylori* seropositivity and five gut genera were performed. The resulting twenty-five-gene signature achieved a leakage-free cross-validated concordance of 0.61 and validated in two of three external cohorts, although a conservative meta-analysis was not statistically significant (pooled hazard ratio 1.19, 95 percent confidence interval 0.96 to 1.47). The prognostic signal converged on a stromal cancer-associated-fibroblast programme that localised to fibroblasts at single-cell resolution and was strongly preserved across three external cohorts, establishing it as externally replicated rather than a single-cohort observation. Adding the signature to stage conferred negligible decision value. The secondary microbiome arm found reduced diversity along the gastritis-to-carcinoma sequence, but tumour-versus-normal differences were confounded with sequencing batch and did not replicate in an independent cohort, and genetic causal testing detected no microbial effect on risk. It was concluded that gastric-cancer prognosis reflects an externally-replicated stromal programme that nonetheless adds little value beyond stage, and that the microbiome associations were observational, not causal.

**Keywords:** gastric cancer; multi-omics; cancer-associated fibroblasts; prognostic signature; Mendelian randomisation; external validation

---

## 1. Introduction

Gastric cancer (GC) remains the fifth most common malignancy and a leading cause of cancer mortality worldwide, with marked histological and molecular heterogeneity (Lauren intestinal/diffuse classes; TCGA EBV/MSI/GS/CIN subtypes) (1,20). Two biological axes dominate current mechanistic thinking: (i) host transcriptomic reprogramming driving proliferation, epithelial–mesenchymal transition (EMT) and immune evasion, and (ii) dysbiosis of the gastric and gut microbiome, classically initiated by *Helicobacter pylori* (3,4). A recurring weakness in the microbiome-cancer literature is the conflation of **association** with **causation**: taxa that differ between tumour and normal tissue are frequently described as drivers without formal causal testing.

This study is organised around a **single primary question — what transcriptomic biology underlies GC prognosis** — with the microbiome and causal analyses as a clearly-secondary, exploratory arm. In the **primary** analysis, converging bulk co-expression, an externally-validated prognostic signature and single-cell localisation identify a **stromal / cancer-associated-fibroblast programme** as the biology underlying prognosis, while we show candidly that the derived gene signature, though independently prognostic, adds little decision value beyond stage; we also characterise the tumour immune microenvironment using deconvolution **validated against measured histological leukocyte fraction**. **Secondarily and exploratorily**, we characterise the tumour-tissue microbiome with **compositionally-aware** statistics and formally test microbial causality by two-sample **Mendelian randomisation** (22) — addressing the field's recurrent conflation of association with causation — reporting a **null** result transparently. Throughout we distinguish findings from the limits of the data and scope every claim to the evidence that supports it.

## 2. Methods

### 2.1 Data sources
- **TCGA-STAD** (host transcriptome): 448 samples (412 primary tumour, 36 solid-tissue normal); STAR gene counts, GENCODE v36; full clinical/molecular annotation (Lauren class, AJCC stage, grade, TCGA molecular subtype, driver mutations, tumour mutational burden [TMB], measured leukocyte/lymphocyte percentage, overall survival).
- **GTEx v10 stomach** (normal reference): 407 samples, TPM.
- **GEO cohorts**: GSE27342, GSE63089 (differential-expression support); **GSE62254 (ACRG, n=300)** (2), **GSE15459 (n=191)**, **GSE84437 (n=431)** as independent survival-validation cohorts.
- **Tissue microbiome (discovery)**: PRJDB20660 (Japan) 16S rRNA (V3–V4) — all **944 libraries reprocessed from raw FASTQ** through DADA2 (25) to a genuine ASV table (897 samples × 314 genera after host/off-target removal); cohort spans the gastric-cancer sequence (299 non-ulcer, 103 ulcer, 219 cancer-adjacent, 323 tumour).
- **Tissue microbiome (validation)**: PRJNA641258 (Italy) 16S rRNA (V3–V4) — 40 gastric biopsies (20 cancer / 20 matched control), sequenced batch-balanced, processed through the identical pipeline.
- **GWAS summary statistics (IEU OpenGWAS)** for MR — see §2.7.

### 2.2 Differential expression and cross-cohort concordance
Host differential expression used the **integrated (ComBat-harmonised) TCGA + GTEx transcriptome**, contrasting TCGA tumours (n=412) against a combined normal reference of TCGA solid-tissue normal (n=36) plus GTEx normal stomach (n=407) — i.e. **412 tumour vs 443 normal** — with limma (34) and a cohort covariate (`~ sample_type + dataset`). Because the harmonised matrix is scale-standardised, significance was assessed by adjusted p (BH) and genes were ranked by the moderated t-statistic for downstream enrichment. The two GEO tumour/normal cohorts (GSE27342, GSE63089) were **co-harmonised within the same phenotype-protected ComBat model** and used as a within-cohort tumour-vs-normal concordance check — not a fully independent hold-out, since the harmonisation model included these samples. Because GTEx contributes only normal samples, batch and biology are partially confounded; the cohort covariate and TCGA's own within-cohort tumour/normal contrast mitigate this, and within-cohort concordance in the GEO cohorts (below) is consistent with a biological rather than batch-driven signal — a supportive check, given the co-harmonisation, rather than strictly independent validation. Lauren diffuse-vs-intestinal DE was computed within TCGA tumours.

### 2.3 Functional enrichment
Over-representation (clusterProfiler (32): GO BP/MF/CC, KEGG) and gene-set enrichment (fgsea (33) against MSigDB Hallmark and C2:CP) were run on the ranked DE statistics for tumour-vs-normal and diffuse-vs-intestinal contrasts.

### 2.4 Immune deconvolution
MCP-counter (30) and xCell (31) were applied to the TCGA expression matrix. Deconvolution validity was assessed by Spearman correlation of estimated lymphocyte/immune scores against the **measured** histological leukocyte percentage. Immune scores were compared across tumour/normal, Lauren class and molecular subtype, and tested for association with overall survival.

### 2.5 Prognostic signature (LASSO-Cox) and external validation
Candidate genes were restricted to those measured across the validation cohorts. Univariable Cox screening (p<0.05) preceded LASSO-Cox (`cv.glmnet`, family="cox") on TCGA tumours, yielding a 25-gene signature. Risk scores (Σ coef·z(expression), standardised within each cohort) were median-split; discrimination was assessed by Harrell's C-index and Kaplan–Meier log-rank tests in TCGA and in each external cohort (GSE62254/ACRG n=300, GSE15459 n=191, GSE84437 n=431; total n=922; overall survival endpoint, complete cases), with multivariable Cox adjustment for stage and age. To obtain a leakage-free performance estimate, the entire pipeline — univariable screen, LASSO tuning and coefficient fitting — was re-run inside a nested cross-validation (outer 5-fold, inner tuning by `cv.glmnet`), so that no gene selection or hyperparameter choice ever saw its own test fold; the nested-CV C-index (0.61), not the apparent full-data C (0.72), is reported as the honest discrimination. External-cohort effect sizes were additionally pooled by a Hartung–Knapp random-effects meta-analysis (metafor (35)).

### 2.6 Co-expression network (WGCNA)
Weighted gene co-expression network analysis used the 5,000 most variable genes (signed-hybrid network, biweight midcorrelation). The soft-thresholding power was chosen as the lowest value achieving scale-free topology fit R²≥0.85 (no forced fallback). Modules were detected with blockwiseModules (28) (minModuleSize 30, deepSplit 2, mergeCutHeight 0.25); module preservation in the external cohorts was quantified with the permutation-based Zsummary statistic (29); module eigengenes were correlated with clinical and immune traits, and each was tested against overall survival (Cox). For the survival-associated module, hub genes were ranked by module membership (MM) and gene significance (GS), and overlap with the prognostic signature was assessed. Robustness of this module to the soft-power choice was tested by rebuilding the network across powers 3–12 and confirming that the hub genes remained co-clustered in a single, survival-associated module.

### 2.7 Survival modelling and nomogram
A clinical Cox model for overall survival was built on **complete cases** (no imputation) using real covariates (Age, AJCC Stage, Grade, TMB, plus hub-gene expression), with backward-AIC selection and a nomogram (calibration and external decision-curve analysis in Supplementary Figure S9). Discrimination used the **Harrell optimism-corrected** bootstrap C-index (B=500), with the backward-AIC selection repeated inside each bootstrap replicate to penalise selection instability; events-per-variable is reported. For the combined clinical+signature model, added value (ΔC-index, likelihood-ratio test, IDI, NRI; IDI/NRI computed at 3 years with survIDINRI and bootstrap p-values), time-dependent AUC and decision-curve analysis were computed both in TCGA (in-sample) and, to avoid training-set circularity, **out-of-sample in ACRG** (§3.7). The prognostic-model development and validation are reported in accordance with the **TRIPOD** (36) (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis) guideline; a completed TRIPOD checklist is provided as Supplementary Table S3.

### 2.8 Tissue-microbiome analysis
Raw FASTQ (944 libraries) were denoised with DADA2 (primer trimming, `truncLen` 260/220, error learning, ASV inference, chimera removal, SILVA v138.1 (26) taxonomy); denoising fidelity was confirmed against the study's published read-tracking (Pearson r=0.98). Host mitochondrial (~34% of reads), chloroplast and non-bacterial ASVs were removed and features agglomerated to genus (897 samples × 314 genera). Crucially, each library's **sequencing flowcell** was recorded, and phenotype was cross-tabulated against flowcell to assess batch confounding. Analyses comprised α-diversity (Observed/Shannon/Simpson; ordered **Jonckheere–Terpstra** across the cascade), β-diversity (Bray–Curtis and Aitchison/CLR; PERMANOVA (27) **with and without adjustment for flowcell**), CLR+Wilcoxon+BH differential abundance, and a cancer-vs-control random-forest classifier accompanied by a **flowcell-prediction batch sanity check**. The tumour signal was then tested for replication in the **independent batch-clean validation cohort** (PRJNA641258, processed identically) across the same contrasts.

### 2.9 Mendelian randomisation
Two-sample MR (TwoSampleMR) tested six exposures against gastric cancer:
anti-*H. pylori* IgG seropositivity (ebi-a-GCST90006910; Butler-Laporte 2020) (14) and the gut genera *Streptococcus*, ***Prevotella 9***, *Veillonella*, *Lactobacillus* (MiBioGen; Kurilshikov 2021) (12) and ***Fusobacterium* A** (abundance in stool; Qin 2022) (15) — all European ancestry. (We name the exact instrument traits: the MiBioGen clade is *Prevotella 9*, a specific SILVA genus-level group narrower than genus *Prevotella* sensu lato, and the Qin-2022 exposure is *Fusobacterium* A in the GTDB nomenclature; accession identities were confirmed against raw OpenGWAS `gwasinfo()` output.) Instruments were LD-clumped (r²<0.001 over a 10,000-kb window, European 1000-Genomes reference), with genome-wide significance (p<5×10⁻⁸) attempted first and relaxed to the locus-wide-suggestive threshold below only where fewer than three independent SNPs were recovered. Because no exposure reached ≥3 genome-wide-significant SNPs, **all six exposures — including anti-*H. pylori* seropositivity — used a locus-wide-suggestive threshold of p<1×10⁻⁵** (recorded per exposure; a choice that can bias toward the null via weak-instrument effects, though the per-SNP F-statistics show the instruments themselves are not weak). Instrument strength is reported as the mean and minimum per-SNP F; we note that these F-statistics are computed on the pre-harmonisation instrument set fetched from OpenGWAS, and that several SNPs are dropped at harmonisation (per-exposure pre-clump and harmonised-used SNP counts are given in Supplementary Table S1). The primary outcome was ancestry-matched European gastric cancer (ebi-a-GCST90018849 (13); Sakaue 2021; 1,029 cases); a larger East-Asian outcome (7,921 cases) was used as a power sensitivity, clearly labelled cross-ancestry. Estimation used IVW, MR-Egger (24), weighted-median and mode methods, with Cochran's Q (heterogeneity), MR-Egger intercept (pleiotropy), Steiger filtering (directionality), leave-one-out and MR-PRESSO (23) sensitivity analyses.

### 2.10 Single-cell validation
A public gastric single-cell RNA-seq dataset (GSE134520; premalignant-to-early-gastric-cancer cascade) was processed in Seurat v4 (37). Genes detected in <3 cells and cells with <200 detected genes were removed; cells were further filtered to 200–6,000 detected genes and mitochondrial content <20%, retaining 43,992 cells. Counts were log-normalised (`LogNormalize`, scale factor 10⁴), the 2,000 most variable features selected, scaled, and reduced by PCA (30 components). Clustering used a shared-nearest-neighbour graph over the first 30 principal components (`FindClusters` resolution 0.5), with UMAP on the same 30 components for visualisation. Clusters were annotated to eight major cell types by canonical markers (Supplementary Figure S1). The prognostic WGCNA red-module and LASSO-signature genes were then mapped to their dominant expressing cell type; to avoid circularity, the analysis was repeated after removing the genes used to annotate the fibroblast cluster (Supplementary Table S5, Supplementary Figure S2).

### 2.11 In-silico drug repurposing
Top up/down differentially expressed genes were queried (Enrichr) against LINCS L1000 and GEO drug-perturbation libraries to nominate compounds whose perturbation signature reverses the tumour signature. Results are hypothesis-generating and not experimentally validated.

### 2.12 Software and reproducibility
All analyses ran in R 4.3.3 (a full `sessionInfo()` dump is provided with the code). Every reported result derives from the public datasets described above and is reproducible from the scripts provided under `analysis/`, with outputs under `results/` and pinned package versions in `package_versions.csv`. Complete-case analysis was used throughout; no missing values were imputed. Full figure provenance — mapping every main and supplementary figure to its source image file(s) and the source data table behind each reported number — is provided in `FIGURE_SOURCES.md`; the composite supplementary figures (S8, S9 and S10) are montages of the individual pipeline-output panels listed there, with no panel redrawn or synthesised.

## 3. Results

### 3.1 A reproducible tumour-versus-normal transcriptional programme
The integrated TCGA+GTEx tumour-versus-normal analysis (412 tumour vs 443 normal) revealed a coherent proliferation-dominated programme. Ranked by tumour-vs-normal moderated t, the top up-regulated genes were mitotic/cell-cycle regulators (KIF14, ECT2, CENPF, TPX2, ASPM, HELLS, CST1, ESM1) together with stromal/EMT markers (COL10A1, WNT2, INHBA), while the top down-regulated genes were differentiated gastric and metabolic genes (parietal-cell ATP4A/ATP4B, GPX3, AQP4, ADH7). GSEA on the integrated ranking showed strong up-regulation of **E2F targets (NES 3.75), G2M checkpoint (NES 3.64) and MYC targets**, with coordinated loss of **oxidative phosphorylation (NES −2.54)** and fatty-acid/bile-acid metabolism (Figure 1A); KEGG/GO over-representation highlighted cell cycle, DNA replication/repair and ribosome biogenesis (up) versus gastric acid secretion and fatty-acid degradation (down). At adjusted p<0.05 a large fraction of genes reached significance (3,722 up / 4,025 down of 12,899); on the scale-standardised harmonised matrix this gate is permissive, so the moderated-t ranking and the enrichment above — not the raw count — are the interpretable outputs. The TCGA-only volcano and the top-gene expression heatmap are shown in Figure 6, and the GO:BP/KEGG over-representation dot plots in Supplementary Figure S4.

![Figure 6](results/figures/Fig6.png)

Figure 6. Differential expression between gastric tumours and normal mucosa and its reproducibility.

(a) Volcano plot showing differentially expressed genes in TCGA-STAD. The x-axis represents log2 fold change (log2FC) and the y-axis represents −log10(adjusted P value). Red points indicate significantly upregulated genes, blue points indicate significantly downregulated genes, and grey points represent genes not meeting the significance thresholds of |log2FC| > 1 and adjusted P < 0.05 (2,134 upregulated and 2,362 downregulated of 21,446 tested features).

(b) Heatmap of the 30 most significantly differentially expressed genes, comprising the 15 most upregulated and 15 most downregulated (ranked by adjusted p value, ties broken by absolute log2 fold change then gene symbol). Rows represent genes and columns represent samples. Values are z-scores of expression, so the colour scale ranges from low relative expression (blue) to high relative expression (red).

(c) Scatter plots comparing the integrated tumour-versus-normal ranking with three reference rankings. Each point represents one of 12,899 shared genes and the black line denotes a linear fit. The left panel compares against TCGA alone and therefore re-uses the discovery tumours, so it assesses internal consistency rather than replication (r = 0.73); the middle and right panels compare against the independent cohorts GSE27342 (r = 0.62) and GSE63089 (r = 0.58).

Figure 1. Transcriptional programmes distinguishing gastric tumours from normal mucosa and diffuse-type from intestinal-type tumours.

(a) Bar plot showing Hallmark gene-set enrichment analysis (GSEA) results for the tumour-versus-normal comparison. The x-axis represents the normalised enrichment score (NES), which is the enrichment score scaled for gene-set size so that sets of different sizes are directly comparable; a positive NES indicates enrichment in tumours and a negative NES indicates enrichment in normal tissue. Bars are coloured by the −log10 adjusted P value, and only gene sets reaching a false-discovery-rate-adjusted P value below 0.05 are shown. Cell-cycle programmes (E2F targets, NES 3.75; G2M checkpoint, NES 3.64) and MYC targets (NES 2.67) are enriched in tumours, whereas oxidative phosphorylation (NES −2.54) and fatty-acid metabolism (NES −2.58) are depleted.

(b) Bar plot showing Hallmark GSEA results for the comparison of diffuse-type with intestinal-type tumours as defined by the Lauren classification, a histological scheme that separates gastric cancers into cohesive gland-forming (intestinal) and poorly cohesive infiltrative (diffuse) types. A positive NES indicates enrichment in diffuse tumours. Epithelial–mesenchymal transition (EMT) is the most strongly enriched programme (NES 3.17, adjusted P = 1.6 × 10⁻⁴⁰), accompanied by inflammatory, interleukin-6–JAK–STAT3 and tumour-necrosis-factor-α–NF-κB signalling.

Figure 4. Immune-cell deconvolution and tissue-microbiome comparison.

(a) Lollipop plot showing the Spearman correlation between each computationally deconvolved immune-cell score and the corresponding pathologist-measured value in TCGA-STAD tumours. Points represent correlation coefficients and horizontal segments connect each point to the null value of zero. Filled points denote correlations significant at P < 0.05.

(b) Bar plots showing mean immune-cell scores across the four molecular subtypes defined by The Cancer Genome Atlas: Epstein-Barr virus-positive (EBV), microsatellite-instable (MSI), genomically stable (GS) and chromosomal-instability (CIN). Each facet corresponds to one deconvolution score (CD8 T cells and T cells from MCP-counter; ImmuneScore from xCell) and bars are coloured by subtype. EBV-positive tumours carry the highest score on all three measures and chromosomal-instability tumours the lowest, with microsatellite-instable and genomically-stable tumours intermediate and their relative order varying by score. All Kruskal-Wallis tests across subtypes remain significant after Benjamini-Hochberg correction (adjusted P < 1 x 10-6). The Lauren histological classification is not the stratifying variable in this panel.

(c) Bar plot summarising tumour-versus-control separation of the tissue microbiome in three cohorts. In the Japanese discovery cohort the apparent separation collapses after adjustment for sequencing flowcell (Bray–Curtis R² 0.065 before adjustment, 0.011 after), in the Italian cohort no separation is detected (R² = 0.018, P = 0.80), and in the Portuguese cohort reduced diversity replicates (Shannon P = 0.004; Bray–Curtis R² = 0.145, P = 0.001). These analyses are presented as a cautionary assessment because batch structure and disease status are partially confounded in the discovery cohort.

Figure 7. Prognostic signature coefficients and external validation.

(a) Bar plot showing the LASSO–Cox regression coefficients of the 25 retained genes. A positive coefficient indicates that higher expression is associated with shorter survival and a negative coefficient indicates the opposite.

(b) Forest plot showing hazard ratios for high- versus low-risk groups (median split of the signature score) across validation cohorts, with squares denoting point estimates and horizontal lines denoting 95% confidence intervals: ACRG/GSE62254 1.90 (1.37-2.62), GSE15459 1.68 (1.11-2.54), GSE84437 1.11 (0.84-1.46). Values are from `results/validation_multi/cindex_HR_summary.csv`; the corresponding Kaplan-Meier curves are panels (c)-(f).

(c) Kaplan-Meier curves for high-risk and low-risk groups in the TCGA-STAD training cohort, defined by a median split of the signature score, with the P value from the log-rank test.

(d) Kaplan-Meier curves for the same median-split groups in the Asian Cancer Research Group cohort (GSE62254, n = 300).

(e) Kaplan-Meier curves for the same median-split groups in GSE15459 (n = 191).

(f) Kaplan-Meier curves for the same median-split groups in GSE84437 (n = 431), in which the signature did not separate the two risk groups. This negative result is reported as observed.

Figure 2. Development and honest evaluation of the 25-gene LASSO–Cox prognostic signature.

(a) Kaplan–Meier survival curves for patients in The Cancer Genome Atlas stomach adenocarcinoma (TCGA-STAD) cohort stratified into high-risk and low-risk groups at the median signature score (n = 383, 156 deaths). The x-axis represents follow-up time in months and the y-axis represents overall survival probability. Shaded bands denote 95% confidence intervals and the P value is from the log-rank test.

(b) Forest plot showing the hazard ratio per one standard deviation increase in signature score in each external validation cohort, adjusted for age and stage, together with the Hartung–Knapp random-effects pooled estimate. Squares represent point estimates, horizontal lines represent 95% confidence intervals, and the diamond represents the pooled estimate (hazard ratio 1.19, 95% confidence interval 0.96–1.47; prediction interval 0.90–1.57). The pooled effect does not reach statistical significance.

(c) Point-and-interval (forest) plot comparing discrimination measured with and without information leakage. The leakage-free estimate is obtained from 20 repeats of 5-fold nested cross-validation, in which gene selection and model fitting occur only inside the training folds (Harrell C = 0.611, 95% confidence interval 0.562–0.659; Uno C = 0.573), whereas the apparent estimate re-uses the same samples for selection and evaluation (C = 0.72). The difference of approximately 0.11 quantifies optimism.

(d) Line plot showing the time-varying hazard ratio for the signature in the Asian Cancer Research Group (ACRG) cohort. The x-axis represents time since surgery in months and the y-axis represents the hazard ratio. The signature is prognostic early (hazard ratio 1.49 at 12 months, 95% confidence interval 1.23–1.80) and attenuates towards the null by 36–60 months, indicating that the proportional-hazards assumption does not hold across the full follow-up period.

Figure 5. Construction and quality assessment of the gene co-expression network.

(a) Cluster dendrogram showing hierarchical clustering of genes by topological overlap, with the resulting module assignments displayed as coloured bars beneath the tree.

(b) Heatmap of correlations between module eigengenes and clinical traits. Rows represent modules and columns represent traits; the colour scale ranges from negative (blue) to positive (red) correlation.

(c) Sensitivity of the prognostic module to the soft-thresholding power, shown as three line plots over the candidate powers 3, 6, 9 and 12. The soft-thresholding power raises correlations to a power so that strong correlations are emphasised and weak ones suppressed, approximating the scale-free topology observed in biological networks; because this choice is made by the analyst, its influence on the result is reported here. The left panel plots the hazard ratio per standard deviation of the module eigengene, with vertical bars denoting 95% confidence intervals and the dashed line marking a hazard ratio of one. The middle panel plots the negative base-ten logarithm of the Cox P value, with the dashed line marking P = 0.05. The right panel plots the fraction of the eight module hub genes that remain co-clustered. Across all four powers the hazard ratio stays between 1.28 and 1.31, the Cox P value stays below 0.05 (maximum 0.0024) and hub co-membership stays at or above 0.875, so the prognostic result does not depend on the particular power chosen. Values are given in Supplementary Table S8.

Figure 3. External preservation of the stromal co-expression module, its prognostic value, and its cellular localisation.

(a) Bar plot showing preservation of the survival-associated (red) module in three external cohorts, quantified by the WGCNA Zsummary statistic. Zsummary is a permutation-based composite of density and connectivity preservation; the dashed line marks Z = 10, the conventional threshold for strong preservation. All three cohorts exceed it (ACRG/GSE62254 Z = 15.9, GSE15459 Z = 16.8, GSE84437 Z = 17.1), indicating that the module is a reproducible feature of gastric tumour transcriptomes rather than a property of the discovery cohort.

(b) Forest plot showing Cox proportional-hazards results for the module eigengene, defined as the first principal component of the module's expression matrix and therefore a single summary value per patient. Points represent the hazard ratio per one standard deviation of eigengene and horizontal bars represent 95% confidence intervals, plotted on a logarithmic scale. The module is prognostic in all three cohorts (ACRG hazard ratio 1.27, 95% confidence interval 1.09-1.49, P = 0.0022; GSE15459 1.55, 1.25-1.92, P = 6.8 x 10-5; GSE84437 1.24, 1.08-1.42, P = 0.0021), including GSE84437, the cohort in which the 25-gene signature did not validate.

(c) Bar plot showing the single-cell localisation of the stromal-module hub genes. For each gene, the bar represents the fraction of its total expression contributed by the cell type in which it is most highly expressed, computed from a public gastric single-cell RNA-sequencing dataset. Bars are coloured by that dominant cell type. Twenty-eight of the twenty-nine hub genes are fibroblast-dominant (median fraction 0.97), placing the prognostic programme in the fibroblast compartment; FAP is the single exception and is endothelial-dominant in this dataset.

Figure 8. Two-sample Mendelian randomisation of microbial exposures on gastric-cancer risk.

Each panel shows the association between single-nucleotide-polymorphism effects on one microbial exposure (x-axis) and the corresponding effects on gastric-cancer risk (y-axis) in a European-ancestry outcome dataset. Each point represents one genetic instrument, with horizontal and vertical bars denoting standard errors. Fitted lines represent the inverse-variance-weighted, MR-Egger, weighted-median, simple-mode and weighted-mode estimators.

(a) Anti-Helicobacter pylori immunoglobulin G seropositivity (17 instruments).

(b) Streptococcus, genus level (15 instruments).

(c) Fusobacterium (23 instruments).

(d) Prevotella (15 instruments).

(e) Veillonella (8 instruments).

(f) Lactobacillus (10 instruments).

No exposure shows evidence of a causal effect, with the smallest inverse-variance-weighted P value being 0.35. These results should be interpreted as exploratory because the number of genetic instruments per exposure is limited (8 to 23), which reduces statistical power and destabilises the MR-Egger estimate.

## 4. Discussion

This study provides a rigorously-scoped, **primarily transcriptomic** account of gastric-cancer prognosis, with a **secondary, exploratory** microbiome and causal arm. The central finding — a **stromal/cancer-associated-fibroblast biology underlying prognosis** — is robust and triangulated across bulk co-expression, an externally-validated signature and single-cell localisation; the accompanying proliferation/DNA-repair transcriptional programme (with diffuse-restricted EMT) and macrophage-weighted immune microenvironment are internally consistent and externally anchored (deconvolution validated against measured leukocyte fraction). In the secondary arm we characterise tumour-microbiome dysbiosis observationally and, rather than assuming causality, test it directly: the MR analysis does not license a causal claim for any tested microbe. A distinct methodological contribution across both arms is therefore to **separate association from causation** rather than conflate them, as is common in this literature.

The principal finding is biological and largely **confirmatory** rather than a novel discovery or a ready-made clinical tool. We do not claim the stromal/CAF-prognosis link as new — it is well established (21) — but provide unusually rigorous, multi-modal **triangulation** of it: the convergence of a power-robust survival-associated WGCNA module, its overlap with the prognostic signature, its single-cell localisation to fibroblasts, and an honest external decision-curve analysis points consistently to the **tumour stroma / cancer-associated fibroblasts** as the biology underlying prognosis in gastric cancer. We are careful **not to overstate the independence** of this convergence: the WGCNA module and the LASSO signature are both derived from the same TCGA cohort, so their agreement is internal consistency rather than independent replication; only 3 of the 25 signature genes (SERPINE1, POSTN, MATN3) fall within the module; and the single-cell step uses canonical fibroblast markers (DCN, LUM, COL1A1/COL1A2, FAP, COL3A1) both to define the fibroblast cluster and as the localisation readout, a partial circularity. The evidential strength here is **consistency across modalities on shared data**, not statistical independence — genuinely independent confirmation would require an external co-expression cohort and marker-free cell-type assignment (a related bulk-plus-single-cell gastric CAF signature has been reported previously, e.g. Li et al., 2023, PMID 36717783, underscoring that this biology is established rather than novel). The 25-gene signature operationalises part of this signal: it is significant in two independent cohorts and remains prognostic after adjustment for stage and age in ACRG, so it is not merely a stage proxy. We are, however, deliberately restrained about its clinical utility: its external discrimination is modest (C 0.58–0.61), it failed in one of three cohorts (GSE84437), and — importantly — its apparent large added value over staging did not replicate out-of-sample (external ΔC +0.005, no decision-curve benefit). The honest reading is that the signature illuminates the stromal biology of GC prognosis but does not, on current evidence, warrant clinical deployment over standard staging. To be explicit about clinical translation: **AJCC TNM staging remains the standard of care**, and the signature, nomogram and combined model presented here are **research instruments** that clarify the biology of prognosis — not validated clinical tools, and not proposed to alter patient management.

The microbiome results carry a cautionary methodological message. Reprocessing the discovery cohort from raw reads revealed that its tumour-versus-normal contrast is **confounded with sequencing batch (5)** (tumour and normal libraries sequenced on separate flowcells), so the apparent tumour "oralization" signal — and a classifier reaching AUC 0.916 — was largely a **batch artefact** (the same features predicted the flowcell at 78% accuracy). Testing across three independent cohorts, the oral-taxa signal proved **comparator-dependent** (up versus batch-confounded adjacent-normal, null versus matched controls, and *down* versus gastritis — the latter reproducing the "*H. pylori* paradox"; Coker 2018; Ferreira 2018), and did **not** generalise. Only **reduced α-diversity (specifically richness) in cancer** replicated consistently — and even there, evenness-weighted Simpson diversity did not track richness, so the durable signal is a loss of rare taxa rather than a uniform diversity collapse. We therefore present the microbiome arm not as a biomarker but as an honest, multi-cohort demonstration that **batch confounding and comparator choice — not a fixed dysbiotic signature — dominate low-biomass tissue-microbiome results**, a caution directly relevant to the many gastric-microbiome "signatures" reported without such controls. Finally, the MR and tissue-microbiome arms are genetically decoupled: the MiBioGen *Streptococcus* instrument indexes faecal genus abundance and has no established relationship to gastric-mucosal *Streptococcus*, so the MR null does not bear on any specific tissue co-occurrence observed here, and neither arm should be read as confirming or refuting the other.

The in-silico drug-repurposing arm (§3.11) connects back to this biology rather than standing apart from it: the nominated PI3K/mTOR and CDK4/6 inhibitor classes (19) reverse the proliferation-dominated tumour-up programme identified in §3.1–3.2, so they are best read as candidate anti-proliferative directions consistent with the transcriptomic biology, not as gastric-cancer-specific or stroma-targeting leads, and they remain unvalidated hypotheses.

> **Clinical implications.** For a treating clinician, the practical takeaway of this work is that **nothing here should change patient management today**. AJCC TNM staging remains the standard of care for gastric-cancer prognostication. The 25-gene signature, the clinical nomogram and the combined model are **research instruments** that illuminate the stromal/fibroblast biology of prognosis; they are not validated clinical decision tools. The signature's out-of-sample added value over staging is decision-analytically negligible (external ΔC +0.005, no decision-curve benefit), its discrimination is modest and time-limited, and the nomogram additionally requires tumour mutational burden — an NGS-derived covariate not routinely available in the community settings where staging is applied. A forward-looking, testable clinical hypothesis follows from the GSE84437 non-validation (§3.4): because a stromal/desmoplastic signature discriminates by capturing variation in invasion-linked stromal content, its incremental value, if any, is most plausible in **earlier-stage disease**, where AJCC stage alone is less discriminating — a hypothesis for prospective evaluation, not a current recommendation.

### 4.1 Limitations
(1) Transcriptomic cohorts span RNA-seq and microarray platforms; batch and platform effects are mitigated by within-cohort standardisation but not eliminated. (2) The 25-gene signature is only modestly and transiently prognostic — leakage-free nested-CV C=0.61 (not the apparent 0.72), a conservative Hartung–Knapp meta-analysis of the three external cohorts is non-significant (pooled HR 1.19, 95% CI 0.96–1.47), the ACRG effect is non-proportional (strong early, null by 3–5 yr), and it adds no decision-analytic value over staging. Bootstrap stability selection (B=200) retained no signature gene at high stability: 13 of the 25 genes were selected in >50% of resamples but **none exceeded 80%** (median selection frequency 0.505; the most stable were NETO2 0.77, EGF 0.75, SRMS 0.73 and SERPINE1 0.70, while the stromal-module members POSTN and MATN3 were selected only 38% and 22% of the time), so individual members should not be over-interpreted; the robust prognostic finding is the externally-preserved co-expression module, not the gene list. (3) The clinical nomogram is single-cohort with modest discrimination (selection-inclusive corrected C≈0.64, EPV≈8) and no external validation; it also requires TMB (an NGS-derived covariate not routinely available), a real-world deployability barrier beyond the statistical ones. (4) In the tumour-tissue microbiome, the discovery cohort has tumour/normal confounded with sequencing batch, and across three cohorts the oral-taxa signal is comparator-dependent and does not generalise (only reduced α-diversity replicates); the arm is therefore exploratory/cautionary rather than a validated biomarker, at genus-level resolution. (5) MR used European instruments and a modestly-powered European outcome; the null results are power-limited, and the microbiome-GWAS/tissue-16S populations differ (gut vs gastric, European vs Japanese). (6) No experimental validation was performed; all inference is computational.

## 5. Conclusion
This study establishes an **externally validated stromal/cancer-associated-fibroblast (CAF) programme** as the reproducible transcriptomic biology underlying gastric-cancer prognosis. The programme is not a single-cohort artefact: its co-expression module is strongly preserved in three independent transcriptomic cohorts (module-preservation Zsummary 15.9–17.1) and remains independently prognostic in each (external eigengene HR/SD 1.24–1.55, all *p*<0.005), and single-cell data localise it unambiguously to the fibroblast compartment. Around this positive finding we draw the boundaries of what the data support with equal rigour. The 25-gene operationalisation of the programme is genuinely but modestly prognostic — its leakage-free discrimination (nested-CV Harrell C 0.61) is well below its apparent performance (C 0.72), its early hazard separation attenuates over time, and it adds little decision-analytic value beyond stage; it is therefore a biological readout, not a stand-alone clinical test.

The two secondary arms are, by design, cautionary rather than confirmatory. The apparent tumour-tissue microbial dysbiosis proved comparator-dependent and confounded by sequencing batch, replicating only where a batch-clean external cohort permits — a concrete illustration of how readily tissue-microbiome signals can be manufactured by design artefacts. The Mendelian-randomisation analysis, with a full sensitivity suite, returned a consistent null across all six exposures and both ancestries: at current instrument and outcome power, the data provide no support for a causal microbial effect on gastric-cancer risk. Taken together, the work delivers one robust, replicated prognostic biology and an explicit, honestly-bounded account of the observational-versus-causal distinction that is too often blurred in this literature — a distinction we regard as central to the study's contribution.

## Declarations

**Ethics approval and consent to participate.** This study used only publicly available, de-identified secondary data from established repositories (TCGA, GTEx, GEO, DDBJ/PRJDB20660, IEU OpenGWAS). No new human participants, human material, or animal work was involved; each source dataset was collected under its own ethics approval and consent. Institutional ethics approval was therefore not required for this secondary analysis of publicly available, de-identified data.


**Data availability.** All primary data are public. Host transcriptome: TCGA-STAD (GDC / UCSC Xena) and GTEx v10 (gtexportal.org). Survival and DE cohorts: GEO GSE27342, GSE63089, GSE62254 (ACRG), GSE15459, GSE84437, and single-cell GSE134520 (ncbi.nlm.nih.gov/geo). Tissue 16S: DDBJ PRJDB20660 (Japan, discovery, raw FASTQ reprocessed), NCBI/ENA PRJNA641258 (Italy; Ravegnini et al., *Int J Mol Sci* 2020, PMC7766162; paired tumour/non-tumour, V3–V4) and PRJNA413125 (Portugal; Ferreira et al. 2018, ref 4; gastritis vs carcinoma, V5–V6) as independent validation cohorts. GWAS summary statistics: IEU OpenGWAS — exposures `ebi-a-GCST90006910`, `ebi-a-GCST90017070`, `ebi-a-GCST90032406`, `ebi-a-GCST90017045`, `ebi-a-GCST90017088`, `ebi-a-GCST90017030`; outcomes `ebi-a-GCST90018849` (European) and `ebi-a-GCST90018629` (East-Asian sensitivity). The harmonised MR instrument tables (per-SNP rsID, effect allele, per-SNP F) and the `gwasinfo()` accession-verification output are archived with the code (`results/mr_real/`) so the exact SNP sets are recoverable even if upstream OpenGWAS entries are later reprocessed.

**Code availability.** All analysis scripts (`analysis/`), result tables and figures (`results/`), the pipeline description (`PIPELINE.md`), the pinned package versions (`package_versions.csv`) and an R `sessionInfo()` dump (`sessionInfo.txt`) are provided in the project repository, archived at a version-tagged snapshot with a persistent DOI: ⟦PLACEHOLDER — REPOSITORY URL + Zenodo (or equivalent) DOI to be minted at submission⟧. Every reported result maps to a named script and output file (see `PIPELINE.md`).

**Funding.** ⟦PLACEHOLDER — funding sources and grant numbers, or: "This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors."⟧

**Competing interests.** The authors declare no competing interests.

**Authors' contributions.** ⟦PLACEHOLDER — assign per author, e.g. Conceptualisation: …; Data curation & analysis: …; Software/pipeline: …; Writing — original draft: …; Writing — review & editing: …; Supervision: … All authors read and approved the final manuscript.⟧

**Acknowledgements.** The results shown here are based in part upon data generated by the TCGA Research Network, the GTEx Consortium, the contributing GEO/DDBJ studies, and the IEU OpenGWAS project; we thank these consortia and the original data generators.

---

## References

*Formatted in Vancouver style.*

1. The Cancer Genome Atlas Research Network. Comprehensive molecular characterization of gastric adenocarcinoma. *Nature* 2014;513:202–209. PMID 25079317.
2. Cristescu R, et al. Molecular analysis of gastric cancer identifies subtypes associated with distinct clinical outcomes (ACRG). *Nat Med* 2015;21:449–456. PMID 25894828.
3. Coker OO, et al. Mucosal microbiome dysbiosis in gastric carcinogenesis. *Gut* 2018;67:1024–1032. PMID 29102920.
4. Ferreira RM, et al. Gastric microbial community profiling reveals a dysbiotic cancer-associated microbiota. *Gut* 2018;67:226–236. PMID 29306925.
5. Png CW, et al. Mucosal microbiome associates with progression to gastric cancer (meta-analysis). *Gut* 2022 / PMC10078126.
6. Derks S, et al. Abundant PD-L1 expression in EBV- and MSI-subtype gastric cancers. *Oncotarget* 2016;7:32925–32932. PMC5923357.
7. Lin CH, et al. Intratumoral immune response to gastric cancer varies by molecular and histologic subtype. *Am J Surg Pathol* 2019;43:851–860. PMID 30969179.
8. Kim YS, et al. Glypican-3 is a metastasis suppressor in gastric cancer; low expression predicts poor outcome. *Oncotarget* 2016 / PMC5190106.
9. Chen Y, et al. SERPINE1 (PAI-1) high expression predicts poor prognosis in gastric cancer. *J Oncol* 2022 / PMC8817868.
10. Zhang X, et al. Periostin (POSTN) predicts poor overall survival in gastric cancer. *J Gastrointest Surg* 2022. PMID 36451060.
11. Li J, et al. COL1A1/COL1A2 expression and prognosis in gastric cancer. *World J Surg Oncol* 2016;14:297. PMID 27894325.
12. Kurilshikov A, et al. Large-scale association analyses identify host factors influencing human gut microbiome composition (MiBioGen). *Nat Genet* 2021;53:156–165.
13. Sakaue S, et al. A cross-population atlas of genetic associations for 220 human phenotypes. *Nat Genet* 2021;53:1415–1424.
14. Butler-Laporte G, et al. GWAS of *Helicobacter pylori* serology. 2020 (IEU OpenGWAS GCST900069xx).
15. Qin Y, et al. Combined effects of host genetics and diet on human gut microbiota. *Nat Genet* 2022 (stool-abundance GWAS).
16. Rao W, et al. Dissecting the causal effects of *Helicobacter pylori* serotypes on gastric cancer risk: a two-sample Mendelian randomization study. *Cureus* 2025;17(7):e89185. doi:10.7759/cureus.89185.
17. Chang Y, et al. Association between gut microbiota and gastric cancers: a two-sample Mendelian randomization study. *Front Microbiol* 2024;15:1383530. doi:10.3389/fmicb.2024.1383530.
18. Zhang P, et al. Dissecting the single-cell transcriptome network underlying gastric premalignant lesions and early gastric cancer (GSE134520). *Cell Rep* 2019;27:1934–1947.
19. Subramanian A, et al. A next generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. *Cell* 2017;171:1437–1452 (LINCS L1000; queried via Enrichr).
20. Laurén P. The two histological main types of gastric carcinoma: diffuse and so-called intestinal-type carcinoma. *Acta Pathol Microbiol Scand* 1965;64:31–49. PMID 14320675.
21. Li X, et al. Single-cell and bulk transcriptomics identify a cancer-associated-fibroblast signature predicting prognosis in gastric cancer. 2023. PMID 36717783.
22. Hemani G, et al. The MR-Base platform supports systematic causal inference across the human phenome (TwoSampleMR). *eLife* 2018;7:e34408. PMID 29846171.
23. Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases (MR-PRESSO). *Nat Genet* 2018;50:693–698. PMID 29686387.
24. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. *Int J Epidemiol* 2015;44:512–525. PMID 26050253.
25. Callahan BJ, et al. DADA2: high-resolution sample inference from Illumina amplicon data. *Nat Methods* 2016;13:581–583. PMID 27214047.
26. Quast C, et al. The SILVA ribosomal RNA gene database project. *Nucleic Acids Res* 2013;41:D590–D596. PMID 23193283.
27. Oksanen J, et al. vegan: Community Ecology Package. R package (PERMANOVA/adonis2, betadisper). 2022.
28. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 2008;9:559. PMID 19114008.
29. Langfelder P, et al. Is my network module preserved and reproducible? *PLoS Comput Biol* 2011;7:e1001057. PMID 21283776.
30. Becht E, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression (MCP-counter). *Genome Biol* 2016;17:218. PMID 27765066.
31. Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. *Genome Biol* 2017;18:220. PMID 29141660.
32. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS* 2012;16:284–287. PMID 22455463.
33. Korotkevich G, et al. Fast gene set enrichment analysis (fgsea). *bioRxiv* 2021; doi:10.1101/060012.
34. Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 2015;43:e47. PMID 25605792.
35. Viechtbauer W. Conducting meta-analyses in R with the metafor package. *J Stat Softw* 2010;36:1–48.
36. Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis (TRIPOD): the TRIPOD statement. *BMJ* 2015;350:g7594. PMID 25569120.
37. Stuart T, et al. Comprehensive integration of single-cell data (Seurat). *Cell* 2019;177:1888–1902. PMID 31178118.

---

## Supplementary Figures

![Supplementary Figure S1](results/scrna/UMAP_celltypes.png)

Supplementary Figure S1. Single-cell reference atlas of the gastric premalignant-to-early-cancer cascade.

Uniform manifold approximation and projection (UMAP) of 43,992 cells from the public dataset GSE134520, coloured by the eight annotated major cell types. UMAP is a non-linear dimensionality-reduction method in which each point represents one cell and proximity reflects transcriptional similarity; the axes carry no independent units. Cell-type proportions are epithelial 68.6%, endothelial 7.6%, myeloid 6.4%, plasma 5.8%, T cell 5.8%, fibroblast 4.2%, mast 1.0% and B cell 0.6%. This atlas provides the cellular reference used for the module-localisation analysis in Figure 3(c).

File: `results/scrna/UMAP_celltypes.png`.

Supplementary Figure S2. Expression of stromal-module hub genes across cell types.

Dot plot of the prognostic co-expression module hub genes across the eight annotated cell types. Dot size represents the fraction of cells in which the gene is detected and dot colour represents scaled mean expression. Hub-gene expression is confined to the fibroblast compartment. Because the genes used to annotate the fibroblast cluster could bias this result, the analysis was repeated after removing them; all twenty-three remaining hub genes stayed fibroblast-dominant.

File: `results/scrna/DotPlot_stromal_module_hub.png`.

Supplementary Figure S3. Candidate compounds from in-silico drug repurposing.

Bar plot of the top compounds whose transcriptional perturbation signature reverses the tumour programme on both arms, that is, repressing tumour-up genes and inducing tumour-down genes, ranked by combined score from Enrichr queried against the LINCS L1000 and Gene Expression Omnibus drug-perturbation libraries. The nominated classes are dominated by phosphoinositide 3-kinase (PI3K)/mechanistic target of rapamycin (mTOR) inhibitors, cyclin-dependent kinase 4/6 (CDK4/6) inhibitors, and fibroblast growth factor receptor (FGFR)/multikinase inhibitors. These predictions are computational and hypothesis-generating only; no compound was experimentally validated, and because the tumour-up programme is proliferation-dominated the nominated agents are broadly anti-proliferative rather than specific to gastric cancer.

File: `results/drug_repurposing_integrated/top_candidate_drugs.png`.

Supplementary Figure S4. Over-representation analysis of tumour-up-regulated genes.

Dot plots of Gene Ontology biological-process and Kyoto Encyclopedia of Genes and Genomes (KEGG) over-representation analysis for the tumour-up gene set. Dot size represents the number of genes in each term and colour represents the adjusted P value. Cell-cycle and nuclear-division terms, together with extracellular-matrix-receptor and cytokine terms, dominate. Over-representation analysis depends on an arbitrary significance threshold and is sensitive to gene-family size; the prominent olfactory and sensory-perception terms reflect the large olfactory-receptor gene family rather than tumour biology. The threshold-free gene-set enrichment analysis in Figure 1 is therefore the primary enrichment result, with this panel provided for completeness.

Files: `results/enrichment/dotplot_GO_BP_UP.png`, `results/enrichment/dotplot_KEGG_UP.png`.

Supplementary Figure S5. Mendelian-randomisation scatter plot for a representative exposure.

Scatter plot of single-nucleotide-polymorphism (SNP) effects on the exposure, anti-Helicobacter pylori immunoglobulin G seropositivity, against their effects on gastric-cancer risk. Each point represents one genetic instrument and error bars represent standard errors on both axes. The five overlaid lines are the slopes estimated by inverse-variance weighting (IVW), MR-Egger, weighted median, weighted mode and simple mode. All slopes are flat and statistically indistinguishable from zero (IVW odds ratio 0.96, 95% confidence interval 0.71-1.30). Equivalent plots for all six exposures in both ancestries are provided in the repository.

File: `results/mr_real/scatter_H__pylori_IgG_seropositivity.png`.

Supplementary Figure S6. Mendelian-randomisation leave-one-out analysis for a representative exposure.

Leave-one-out inverse-variance-weighted estimates for anti-Helicobacter pylori immunoglobulin G seropositivity. Each row shows the pooled causal estimate recomputed with one instrument removed; points represent the estimate and horizontal bars the 95% confidence interval. Every estimate straddles the null, confirming that no individual instrument drives the result.

File: `results/mr_real/loo_H__pylori_IgG_seropositivity.png`.

Supplementary Figure S7. Compositional differential abundance of tissue microbiota.

Bar plots of centred-log-ratio effect sizes for the two discovery-cohort contrasts, tested by Wilcoxon rank-sum with Benjamini-Hochberg correction. The centred-log-ratio transform is used because sequencing yields relative rather than absolute abundances, so raw proportions are not independent. The left panel shows control versus cancer-adjacent mucosa (44 of 61 genera at q < 0.05) and the right panel the paired cancer-adjacent versus tumour contrast (18 of 61 genera at q < 0.05). Bars are coloured by direction where q < 0.05, with red indicating enrichment and blue depletion, and the twelve most enriched and depleted genera are shown per panel. These genus-level shifts are reported for completeness only; as detailed in the Results, the tumour contrast is confounded with sequencing batch and does not replicate in an independent batch-clean cohort, so these taxa are not proposed as biomarkers.

Files: `results/microbiome_biomarker/04a_DA_control_vs_GCN.csv`, `04b_DA_GCN_vs_GCT_paired.csv`.

Supplementary Figure S8. Immune infiltration with deconvolution validated against pathology.

(a) Scatter plot of the deconvolution-derived T-cell score against the measured histological leukocyte fraction scored by a pathologist (Spearman rho = 0.67, P = 3.6 x 10-36), establishing that the expression-based estimates track true tissue composition.

(b) Comparison of immune compartments between tumour and normal tissue. Enrichment is dominated by the macrophage and monocyte lineage, with no net gain in CD8-positive T cells.

(c) Immune infiltration across the four molecular subtypes defined by The Cancer Genome Atlas. Epstein-Barr virus-positive tumours are the most infiltrated and chromosomal-instability tumours the least (Kruskal-Wallis P < 1 x 10-6).

(d) Association between CD8-positive T-cell score and overall survival. The score is not prognostic in this cohort (Cox hazard ratio 1.04, P = 0.41); this null result is reported as observed.

Files: `results/plots/Immune_*.png`.

Supplementary Figure S9. Clinical nomogram, calibration, and external decision-curve analysis.

(a) Nomogram combining clinical covariates with the signature score for prediction of 1-, 3- and 5-year overall survival. A nomogram converts each predictor into points on a common scale whose total maps to a predicted survival probability.

(b) Calibration plots at 1, 3 and 5 years, computed in the development cohort. Calibration compares predicted against observed survival; the diagonal represents perfect agreement.

(c) Decision-curve analysis performed out-of-sample in the Asian Cancer Research Group cohort. Decision-curve analysis plots net benefit against the threshold probability at which a clinician would act. Both the clinical model and the combined model carry net benefit relative to treat-all and treat-none strategies, but their curves are essentially superimposed across the 5-40% threshold range, indicating that adding the signature confers no incremental net benefit over standard staging. The nomogram is presented as an illustrative research tool, not as a validated decision instrument.

Files: `results/nomogram_combined/`, `results/external_utility_ACRG/DCA_external.png`.

Supplementary Figure S10. Mendelian-randomisation leave-one-out analysis for all six exposures.

Leave-one-out inverse-variance-weighted estimates for each of the six microbial exposures. Removing any single instrument leaves every pooled estimate straddling the null, confirming that no individual instrument drives any result and that the overall null is not an outlier artefact.

Files: `results/mr_real/loo_*.png`.

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


---

## Appendix: Data-Acquisition & Reproducibility Checklist

*Gastric-cancer multi-omics study. Every dataset used, its accession, access route,
version/build, the raw→processed entry point, and the pipeline script that consumes it.
Datasets are listed in the order the pipeline uses them. This lets a reader reacquire
every input from its primary repository and re-enter the pipeline at the correct stage.*

## A. Host transcriptome (differential expression, WGCNA, signature, immune, drug-repurposing)

- [ ] **TCGA-STAD** — primary tumour + solid-tissue normal RNA-seq.
  - Accession/portal: GDC / UCSC Xena (project TCGA-STAD).
  - Build: STAR gene counts, GENCODE v36; 448 samples (412 tumour, 36 normal).
  - Clinical: Lauren class, AJCC stage, grade, TCGA molecular subtype, driver mutations,
    TMB, measured leukocyte/lymphocyte %, overall survival.
  - Raw→processed entry point: count matrix → `results/rdata/tcga_processed.RData`.
  - Consuming scripts: DEG (`analysis/…DEG…`), WGCNA (`analysis/…wgcna…`),
    immune deconvolution (`08_immune_deconvolution.R`), enrichment,
    microbe-response GSEA (`31_microbe_response_enrichment.R`).
- [ ] **GTEx v10** — normal stomach reference (DE denominator).
  - Portal: gtexportal.org. Combined with TCGA normals for the integrated tumour-vs-normal contrast.

## B. External survival / DE validation cohorts (transcriptome)

- [ ] **GSE62254 (ACRG)** — 300 patients; Affymetrix. Primary external validation of the signature.
  - Access: GEO. Local: `data/geo/GSE62254.rda`. Consumer: `12_multicohort_validation.R`, `32_nested_cv_signature.R`.
- [ ] **GSE15459** — 191 patients. Access: GEO. Local: `data/geo/GSE15459_es.rds` + `GSE15459_outcome.xls`.
- [ ] **GSE84437** — 431 patients. Access: GEO. Local: `data/geo/GSE84437_es.rds`. (Reported negative — retained transparently.)
- [ ] **GSE27342, GSE63089** — independent DE-concordance cohorts. Access: GEO.
  Consumer: DEG cross-cohort concordance (`results/tables/DEG_GEO_*.csv`).

## C. Single-cell transcriptome (module localisation)

- [ ] **GSE134520** — premalignant-to-early gastric-cancer scRNA-seq.
  - Access: GEO. Processing: Seurat v4 (QC 200–6,000 genes, MT<20% → 43,992 cells;
    LogNormalize; 2,000 HVGs; PCA 30; SNN res 0.5; UMAP). 8 cell types by canonical markers.
  - Consumer: scRNA localisation → `results/scrna/`.

## D. Tissue microbiome (16S)

- [ ] **PRJDB20660 (Japan, discovery)** — 944 libraries, 16S V3–V4.
  - Access: DDBJ, **raw FASTQ**. Processing: DADA2 from raw → ASV table
    (897 samples × 314 genera after host/off-target removal).
  - Read-tracking validation: per-sample counts vs published Supplementary Table 3, Pearson r=0.98
    (`results/microbiome_biomarker/00_readtracking_concordance.csv`).
  - Consumer: `23_dada2_16S.R` → `24_microbiome_real.R`; reference DB `silva_nr99_v138.1`.
- [ ] **PRJNA641258 (Italy, validation)** — paired tumour/non-tumour, V3–V4. Access: NCBI/ENA. (Ravegnini 2020, PMC7766162.)
- [ ] **PRJNA413125 (Portugal, validation)** — gastritis vs carcinoma, V5–V6. Access: NCBI/ENA. (Ferreira 2018.)
  - **Batch caveat:** the discovery tumour/normal split is confounded with sequencing flowcell (§3.8);
    validation cohorts are the batch-clean test.

## E. GWAS summary statistics (Mendelian randomisation)

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

## F. Environment & reproducibility

- [ ] R 4.3.3; pinned package versions in `package_versions.csv`; full `sessionInfo()` with the code.
- [ ] Complete-case analysis throughout; no imputation.
- [ ] Bundled `r_env/` has a build-time hardcoded path — after relocation use the `R_HOME`-override
      invocation or a system R (documented in `PIPELINE.md`).
- [ ] Every `results/` output maps to its generating script in `PIPELINE.md`;
      figure→source mapping in `FIGURE_SOURCES.md`.

## G. Analyses NOT performed (stated to prevent over-reading)

- [ ] **No genuinely paired per-patient tissue-microbiome + host-transcriptome cohort exists** in this
      study (16S and transcriptome are from different, unpaired cohorts). A per-genus host-gene MaAsLin2
      association was run in a superseded pipeline on a bridged, non-verifiable sample set; its outputs are
      not reproducible from retained files, so it is excluded and no host-gene↔genus claim is made. The
      microbe-response result (§3.9) is a pathway-level GSEA on the TCGA tumour-vs-normal ranking only —
      hypothesis-generating context, not a host–microbiome axis.

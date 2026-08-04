Figure 1. Transcriptional programmes distinguishing gastric tumours from normal mucosa and diffuse-type from intestinal-type tumours.

(a) Bar plot showing Hallmark gene-set enrichment analysis (GSEA) results for the tumour-versus-normal comparison. The x-axis represents the normalised enrichment score (NES), which is the enrichment score scaled for gene-set size so that sets of different sizes are directly comparable; a positive NES indicates enrichment in tumours and a negative NES indicates enrichment in normal tissue. Bars are coloured by the −log10 adjusted P value, and only gene sets reaching a false-discovery-rate-adjusted P value below 0.05 are shown. Cell-cycle programmes (E2F targets, NES 3.75; G2M checkpoint, NES 3.64) and MYC targets (NES 2.67) are enriched in tumours, whereas oxidative phosphorylation (NES −2.54) and fatty-acid metabolism (NES −2.58) are depleted.

(b) Bar plot showing Hallmark GSEA results for the comparison of diffuse-type with intestinal-type tumours as defined by the Lauren classification, a histological scheme that separates gastric cancers into cohesive gland-forming (intestinal) and poorly cohesive infiltrative (diffuse) types. A positive NES indicates enrichment in diffuse tumours. Epithelial–mesenchymal transition (EMT) is the most strongly enriched programme (NES 3.17, adjusted P = 1.6 × 10⁻⁴⁰), accompanied by inflammatory, interleukin-6–JAK–STAT3 and tumour-necrosis-factor-α–NF-κB signalling.

Figure 2. Development and honest evaluation of the 25-gene LASSO–Cox prognostic signature.

(a) Kaplan–Meier survival curves for patients in The Cancer Genome Atlas stomach adenocarcinoma (TCGA-STAD) cohort stratified into high-risk and low-risk groups at the median signature score (n = 383, 156 deaths). The x-axis represents follow-up time in months and the y-axis represents overall survival probability. Shaded bands denote 95% confidence intervals and the P value is from the log-rank test.

(b) Forest plot showing the hazard ratio per one standard deviation increase in signature score in each external validation cohort, adjusted for age and stage, together with the Hartung–Knapp random-effects pooled estimate. Squares represent point estimates, horizontal lines represent 95% confidence intervals, and the diamond represents the pooled estimate (hazard ratio 1.19, 95% confidence interval 0.96–1.47; prediction interval 0.90–1.57). The pooled effect does not reach statistical significance.

(c) Bar plot comparing discrimination measured with and without information leakage. The leakage-free estimate is obtained from 20 repeats of 5-fold nested cross-validation, in which gene selection and model fitting occur only inside the training folds (Harrell C = 0.611, 95% confidence interval 0.562–0.659; Uno C = 0.573), whereas the apparent estimate re-uses the same samples for selection and evaluation (C = 0.72). The difference of approximately 0.11 quantifies optimism.

(d) Line plot showing the time-varying hazard ratio for the signature in the Asian Cancer Research Group (ACRG) cohort. The x-axis represents time since surgery in months and the y-axis represents the hazard ratio. The signature is prognostic early (hazard ratio 1.49 at 12 months, 95% confidence interval 1.23–1.80) and attenuates towards the null by 36–60 months, indicating that the proportional-hazards assumption does not hold across the full follow-up period.

Figure 3. Weighted gene co-expression network analysis identifies a stromal module associated with survival.

(a) Cluster dendrogram of genes with module assignments shown as coloured bars beneath. Genes are clustered by topological overlap, a similarity measure that counts shared network neighbours in addition to direct correlation, so that genes in the same module are densely interconnected. Each colour denotes one module; grey denotes genes not assigned to any module.

(b) Heatmap of module–trait relationships. Rows represent module eigengenes, defined as the first principal component of each module's expression matrix and therefore a single summary profile for the module, and columns represent clinical traits. Each cell reports the Pearson correlation coefficient with its associated P value, with the colour scale ranging from negative correlation (blue) to positive correlation (red).

(c) Forest plot of Cox proportional-hazards results for the stromal module eigengene in three independent validation cohorts. The module remains prognostic in GSE84437 (hazard ratio 1.24 per standard deviation, 95% confidence interval 1.08–1.42, P = 0.0021), the cohort in which the 25-gene signature failed, indicating that the module captures prognostic information the sparse signature does not.

Figure 4. Immune-cell deconvolution and tissue-microbiome comparison.

(a) Lollipop plot showing the Spearman correlation between each computationally deconvolved immune-cell score and the corresponding pathologist-measured value in TCGA-STAD tumours. Points represent correlation coefficients and horizontal segments connect each point to the null value of zero. Filled points denote correlations significant at P < 0.05.

(b) Box plots showing immune scores stratified by Lauren histological subtype. Boxes span the interquartile range, the central line denotes the median, whiskers extend to 1.5 times the interquartile range, and individual points denote samples beyond that range.

(c) Bar plot summarising tumour-versus-control separation of the tissue microbiome in three cohorts. In the Japanese discovery cohort the apparent separation collapses after adjustment for sequencing flowcell (Bray–Curtis R² 0.065 before adjustment, 0.011 after), in the Italian cohort no separation is detected (R² = 0.018, P = 0.80), and in the Portuguese cohort reduced diversity replicates (Shannon P = 0.004; Bray–Curtis R² = 0.145, P = 0.001). These analyses are presented as a cautionary assessment because batch structure and disease status are partially confounded in the discovery cohort.

Figure 5. Construction and quality assessment of the gene co-expression network.

(a) Cluster dendrogram showing hierarchical clustering of genes by topological overlap, with the resulting module assignments displayed as coloured bars beneath the tree.

(b) Heatmap of correlations between module eigengenes and clinical traits. Rows represent modules and columns represent traits; the colour scale ranges from negative (blue) to positive (red) correlation.

(c) Line plots showing selection of the soft-thresholding power. The left panel plots the scale-free topology model fit (signed R²) against candidate powers and the right panel plots mean connectivity. The soft-thresholding power raises correlations to a power so that strong correlations are emphasised and weak ones suppressed, approximating the scale-free topology observed in biological networks; the chosen power is the lowest value at which the model fit plateaus.

Figure 6. Differential expression between gastric tumours and normal mucosa and its reproducibility.

(a) Volcano plot showing differentially expressed genes in TCGA-STAD. The x-axis represents log2 fold change (log2FC) and the y-axis represents −log10(adjusted P value). Red points indicate significantly upregulated genes, blue points indicate significantly downregulated genes, and grey points represent genes not meeting the significance thresholds of |log2FC| > 1 and adjusted P < 0.05 (2,134 upregulated and 2,362 downregulated of 21,446 tested features).

(b) Heatmap of the 30 most significantly differentially expressed genes, comprising the 15 most upregulated and 15 most downregulated. Rows represent genes and columns represent samples. Values are z-scores of expression, so the colour scale ranges from low relative expression (blue) to high relative expression (red).

(c) Scatter plots comparing the integrated tumour-versus-normal ranking with three reference rankings. Each point represents one of 12,899 shared genes and the black line denotes a linear fit. The left panel compares against TCGA alone and therefore re-uses the discovery tumours, so it assesses internal consistency rather than replication (r = 0.73); the middle and right panels compare against the independent cohorts GSE27342 (r = 0.62) and GSE63089 (r = 0.58).

Figure 7. Prognostic signature coefficients and external validation.

(a) Bar plot showing the LASSO–Cox regression coefficients of the 25 retained genes. A positive coefficient indicates that higher expression is associated with shorter survival and a negative coefficient indicates the opposite.

(b) Forest plot showing hazard ratios per standard deviation of signature score across validation cohorts, with squares denoting point estimates and horizontal lines denoting 95% confidence intervals.

(c) Kaplan–Meier curves for high-risk and low-risk groups defined by the median signature score in each external cohort, with P values from the log-rank test.

Figure 8. Two-sample Mendelian randomisation of microbial exposures on gastric-cancer risk.

Scatter plots showing the association between single-nucleotide-polymorphism effects on each microbial exposure (x-axis) and the corresponding effects on gastric-cancer risk (y-axis). Each point represents one genetic instrument, with horizontal and vertical bars denoting standard errors. Fitted lines represent the inverse-variance-weighted, MR-Egger, weighted-median, simple-mode and weighted-mode estimators. No exposure shows evidence of a causal effect, with the smallest inverse-variance-weighted P value being 0.35. These results should be interpreted as exploratory because the number of genetic instruments per exposure is limited (8 to 23), which reduces statistical power and destabilises the MR-Egger estimate.
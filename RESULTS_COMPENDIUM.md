# Gastric-Cancer Multi-Omics — Complete Results Compendium

*Every analysis in the project, its result, and what it means. All numbers are copied
from the result files under `results/`; each block names its source file. This is the
plain-language companion to `PAPER.md` — the paper states these results formally; this
document explains each one.*

---

## PART 1 — TRANSCRIPTOMIC CORE (the primary story)

### 1.1 Differential expression: tumour vs normal
**Result.** Integrating TCGA tumours (n=412) against a combined TCGA+GTEx normal reference
(n=443), 3,722 genes were up-regulated and 4,025 down-regulated (of 12,899 tested; BH
p<0.05). Top up: mitotic/cell-cycle genes (KIF14, ECT2, CENPF, TPX2, ASPM); top down:
differentiated gastric genes (ATP4A/ATP4B parietal-cell pumps, GPX3, AQP4, ADH7).
*Source: `results/tables/DEG_integrated_TCGA_GTEx.csv`.*

**What it means.** Gastric tumours switch on a proliferation machine and switch off the
specialised machinery of normal stomach lining (acid-secreting parietal cells, metabolism).
This is the expected core biology of carcinoma — it confirms the data behave correctly
before anything subtler is asked of them.

**The honesty check that matters here.** The test-statistic inflation is large (λ=17.3),
which in a GWAS would signal batch artefact. Two controls prove it is real biology, not
batch: (i) a **permutation null** (shuffling tumour/normal labels, B=100) collapses λ to
≈1.06 — an artefact would have stayed inflated; (ii) the **TCGA-only** contrast (no GTEx)
independently reproduces 100% of the top 100–200 integrated genes and ≥98.6% of the top
500 (99.5% overall concordance). *Sources: `DEG_permutation_null_lambda.csv`,
`DEG_TCGAonly_replication.csv`.* → The inflation reflects a genuinely huge tumour-vs-normal
difference, not a technical bias.

### 1.2 Functional enrichment (what pathways are on/off)
**Result.** GSEA on the tumour-vs-normal ranking: **E2F targets NES +3.75**, G2M checkpoint
+3.64, MYC targets +2.67 (all up); **oxidative phosphorylation NES −2.54**, fatty-acid
metabolism −2.58 (down). *Source: `results/enrichment_integrated/fgsea_Hallmark_integrated.csv`.*

**What it means.** The up-programmes are all "cell division" pathways (E2F and G2M are the
transcription factors and checkpoints that drive the cell cycle). The down-programmes are
energy metabolism — tumours rewire away from normal oxidative metabolism. Textbook cancer
metabolism, quantified.

### 1.3 The diffuse subtype is EMT-driven
**Result.** Diffuse-vs-intestinal GSEA: **EMT is the single top Hallmark set in diffuse
tumours (NES 3.17, adjusted p=1.6×10⁻⁴⁰)**, with TGF-β and stromal/inflammatory programmes;
intestinal tumours are dominated by proliferative E2F/MYC. *Source:
`results/enrichment/GSEA_Hallmark_DiffuseVsIntestinal.csv`.*

**What it means.** Lauren's two histological types of gastric cancer have distinct biology:
diffuse tumours are mesenchymal/stromal (EMT = cells losing epithelial identity and gaining
migratory, invasive character), intestinal ones are proliferation-led. This is why the
stromal/CAF story (below) concentrates in diffuse tumours.

### 1.4 Immune microenvironment, validated against pathology
**Result.** Deconvolution (MCP-counter/xCell) validated against **measured** histological
leukocyte fraction: T-cell score vs leukocyte % Spearman **ρ=0.67 (p=3.6×10⁻³⁶)**. Tumours
are enriched for the **macrophage/monocyte** compartment (p=2.5×10⁻⁴), not CD8 T-cells.
EBV/MSI molecular subtypes are the most immune-infiltrated (KW p<10⁻⁶). *Source:
`results/immune/`.*

**What it means.** You can estimate immune-cell content from bulk expression, but it must be
validated — and here it is, against the pathologist's actual leukocyte count (ρ=0.67 is a
strong agreement). The tumours are "macrophage-hot, T-cell-cold," and the immune-hot
subtypes (EBV/MSI) match established immunotherapy biology. An honest null is reported: CD8
score is not prognostic here (HR 1.04, p=0.41).

---

## PART 2 — THE PROGNOSTIC SIGNATURE (and why it is reported honestly)

### 2.1 The 25-gene LASSO-Cox signature
**Result.** A 25-gene signature trained on TCGA stratifies survival (apparent C-index 0.72), but the **leakage-free nested cross-validation C-index is 0.611 (95% CI 0.562–0.659)** — Uno's C 0.573. *Sources: `results/validation/signature_coefficients.csv`, `results/nested_cv/performance.csv`.*

**What it means — the single most important methodological point in the paper.** If you pick genes and tune a model on the whole dataset and then score it on that same data, you get an *optimistic* number (0.72) because the model has effectively seen the answers. The honest way is nested CV: rebuild the *entire* pipeline (gene screen + LASSO) inside each fold so it never sees its own test data. That gives 0.61. The 0.11 gap **is** the leakage most papers silently keep. (C-index 0.5 = coin flip, 1.0 = perfect; 0.61 = modest but real.)

### 2.2 External validation in 3 independent cohorts
| Cohort | Platform | N (events) | C-index | HR high-vs-low (95% CI) | p |
|---|---|---|---|---|---|
| ACRG / GSE62254 | Affymetrix | 300 (152) | 0.61 | 1.90 (1.37–2.62) | 8.6e-5 |
| GSE15459 | Affymetrix | 191 (95) | 0.58 | 1.68 (1.11–2.54) | 0.014 |
| GSE84437 | Illumina | 431 (207) | 0.53 | 1.11 (0.84–1.46) | 0.46 (null) |
*Source: `results/validation_multi/cindex_HR_summary.csv`.*

**What it means.** The signature reproduces in 2 of 3 external cohorts — genuinely independent patients and platforms. Validating externally at all puts it ahead of most signature papers; it is also honest about the third cohort failing.

### 2.3 The meta-analysis is NOT significant (key honesty point)
**Result.** Pooling the three external cohorts by a conservative Hartung–Knapp random-effects meta-analysis: **pooled HR 1.19 (95% CI 0.96–1.47), p=0.073 — not significant** (I2=19%). *Source: `results/meta_HK/meta_result.csv`.*

**What it means.** Combined under a conservative model, the effect does not clear significance. Most papers would report only the two positive cohorts; reporting the non-significant pooled estimate chooses rigour over a cleaner headline.

### 2.4 Why GSE84437 failed — tested, not assumed
**Result.** GSE84437 is overwhelmingly deeply-invasive disease (67% pT4, 89% pT3–T4). Stratifying by pT stage does **not** rescue the signal: C-index <0.5 in every stratum (early pT1–T3 C=0.44; pT4 C=0.48; pT2–T3 C=0.44). *Source: `results/validation_multi/GSE84437_Tstage_stratified.csv`.*

**What it means.** A stromal signature discriminates by measuring *variation* in desmoplastic content; a nearly all-late-stage cohort has little such variation (range restriction). That is the hypothesis — but instead of asserting it, the paper tested it by stratifying and reports honestly that stratification does not recover discrimination.

### 2.5 The signature effect is time-limited (non-proportional hazards)
**Result.** In ACRG the proportional-hazards assumption is violated (cox.zph p=0.003): prognostic **early** (HR/SD 1.49 at 12 mo) and **attenuating to null** (1.03 at 36 mo, 0.87 at 60 mo). *Source: `results/timevarying_ACRG/`.*

**What it means.** The signature flags patients who die *sooner*, but the gap closes over time — an early-hazard marker, not a durable lifelong risk gradient.

### 2.6 Signature stability — no single gene is robust
**Result.** Bootstrap stability selection (B=200): 13/25 genes selected in >50% of resamples, **none >80%** (median 0.505). Most stable: NETO2 0.77, EGF 0.75, SRMS 0.73, SERPINE1 0.70. *Source: `results/signature_stability/stability_summary.csv`.*

**What it means.** Resample the data and you get a *different* gene list each time — the exact genes are not reproducible, though the *module* they represent is. Hence the durable finding is the co-expression module, not the gene list.

### 2.7 Coefficients are penalised weights, not prognostic directions
**Result.** e.g. GPC3 carries a positive (higher-risk) LASSO weight though generally reported as a metastasis suppressor; HBB (haemoglobin) likely reflects blood contamination. *Source: `signature_coefficients.csv`.*

**What it means.** LASSO weights are shaped by correlations among all predictors, so a gene's sign can flip relative to its univariate biology — flagged so no reader mis-reads a coefficient as a standalone claim.

---

## PART 3 — THE STROMAL/CAF PROGRAMME (the primary biological finding)

### 3.1 WGCNA co-expression module
**Result.** Ten co-expression modules; the **red** module is the most prognostic (univariable HR 1.31/SD, 95% CI 1.12–1.53, p=9.3e-4). Its hub genes are canonical CAF/stromal-EMT markers (CDH11, COL1A2, FNDC1, SPARC, LUM, BGN, FAP, POSTN). *Source: `results/wgcna_real/ME_survival_cox.csv`.*

**What it means.** WGCNA groups genes that move together across patients into "modules." One module — full of fibroblast/collagen genes — tracks survival. This is the unsupervised discovery of the stromal-prognosis biology, independent of the supervised signature.

### 3.2 The honest stage caveat
**Result.** Adjusting for stage, age, leukocyte fraction and Lauren subtype, the module attenuates to HR 1.35/SD (95% CI 0.95–1.90, p=0.090) — no longer significant; stage dominates (HR 1.82, p=5e-4). *Source: `results/wgcna_real/ME_survival_cox_adjusted.csv`.*

**What it means.** The stromal programme is entangled with tumour stage — it does not clearly add prognostic information *independent* of stage. Disclosed rather than buried.

### 3.3 It is a stromal signal, not a purity artefact
**Result.** The risk score correlates with stromal content (xCell StromaScore rho=+0.39, p=1.7e-16) and negatively with tumour purity (ABSOLUTE rho=-0.20, p=8.9e-4), yet its prognostic effect **survives purity adjustment** (HR 2.97/SD, p=7.6e-16; purity itself non-significant, p=0.35). *Source: `results/purity/`.*

**What it means.** One could object "this just measures low tumour purity." The purity-adjusted Cox rules that out: the signal remains strong with purity in the model. It co-localises with stroma but is not reducible to it.

### 3.4 The module is externally preserved — the paper's most robust result
**Result.** WGCNA module-preservation in 3 external cohorts: **Zsummary 15.9 / 16.8 / 17.1** (all > the Z=10 "strong" threshold). The eigengene is independently prognostic in all three (HR/SD 1.27 ACRG, 1.55 GSE15459, 1.24 GSE84437, all p<0.005). *Source: `results/module_preservation/`.*

**What it means.** The module is not a TCGA quirk: its wiring reproduces in independent cohorts (Z>10 is a strong-preservation standard), and it predicts survival in each. This is what turns "an observation in one dataset" into "an externally-replicated stromal prognostic programme" — the strongest claim in the paper.

### 3.5 Single-cell localisation (non-circular)
**Result.** In 43,992 single cells (GSE134520), 28/29 hub genes are fibroblast-dominant. After removing the 6 genes used to annotate the fibroblast cluster, **all 23 remaining hub genes stay fibroblast-dominant** (median fraction 0.96). *Source: `results/scrna/gene_dominant_celltype_noncircular.csv`.*

**What it means.** Bulk data can't tell you *which cell* expresses a gene; single-cell can. The stromal module lights up specifically in fibroblasts. The circularity check (removing annotation genes) proves this isn't a definitional loop — the localisation holds on independent genes. Honest limits: FAP (an activated-CAF marker) didn't localise to fibroblasts, and this is early-stage tissue (fibroblasts only ~4% of cells), so "fibroblast/stromal" is claimed, not "activated CAF."

---

## PART 4 — CLINICAL TRANSLATION (honest about limited value)

### 4.1 Nomogram (complete-case clinical model)
**Result.** On complete cases (N=199, no imputation), backward-AIC selected Age, Stage, Grade, TMB; optimism-corrected C-index **0.636 (95% CI 0.563–0.704)**, EPV≈8. Only Stage was a stable predictor (selected 93% of bootstraps). *Source: `results/nomogram/`.*

**What it means.** A standard clinical prognostic tool, built honestly (selection repeated inside every bootstrap so instability is penalised). Discrimination is modest and only stage is a reliable predictor — presented as illustrative, not a validated instrument.

### 4.2 Does the signature add value over staging? (the decisive test)
**Result.** *In-sample* (TCGA) the combined model looks great: ΔC +0.114, LRT p=3e-19. But **out-of-sample in ACRG** the added value collapses: clinical C 0.711 vs combined 0.716, **ΔC +0.005**, and **decision-curve analysis shows no net benefit**. *Source: `results/external_utility_ACRG/`.*

**What it means.** This is the most clinically important honesty in the paper. In-sample, the signature seems to add a lot — but that is circularity (it was trained there). Tested on independent data, it adds essentially nothing a clinician could act on (ΔC +0.005 is negligible; DCA = no decision benefit). The paper states plainly: **staging remains standard of care; this is a research instrument, not a clinical test.**

---

## PART 5 — SECONDARY ARM: TUMOUR MICROBIOME (cautionary/negative)

### 5.1 Raw reprocessing and the batch confound
**Result.** 944 libraries reprocessed from raw FASTQ (DADA2/SILVA) → 897 samples × 314 genera; denoising validated (read-tracking Pearson r=0.98). But **tumour and normal libraries were sequenced on separate flowcells** — of 216 tumour/adjacent-normal pairs, only 1 shares a flowcell. *Sources: `results/microbiome_biomarker/00_readtracking_concordance.csv`, `01_confound_crosstab_final.csv`.*

**What it means.** The pipeline is real and faithful (r=0.98). But tumour status is almost perfectly confounded with sequencing batch — so *any* tumour-vs-normal microbial difference could be batch, not biology. Everything downstream is interpreted through this caveat.

### 5.2 What survives batch vs what is artefact
**Result.** Reduced **richness** along the gastritis→cancer cascade survives (Jonckheere–Terpstra Z=−4.94, p=7.7e-7; Cliff's δ=−0.27), but **Simpson (evenness) diversity does not** (δ=−0.07, p=0.17). β-diversity separation **collapses** under flowcell adjustment (Bray R² 0.065→0.011), and the AUC-0.916 cancer classifier is **batch-driven** (same features predict flowcell at 78% vs 55% baseline). *Sources: `02_alpha_effectsizes_cascade.csv`, `02_beta_permanova.csv`, `05_rf_metrics_and_batch_sanity.csv`.*

**What it means.** Loss of rare taxa (richness) is a real, batch-robust signal; the impressive-looking cancer classifier is a batch artefact. Distinguishing the two is the whole point of the arm.

### 5.3 Differential abundance (CLR)
**Result.** Compositional CLR differential abundance: 44/61 genera differ control-vs-cancer-adjacent, 18/61 in the paired cancer-adjacent-vs-tumour contrast. *Sources: `04a_DA_control_vs_GCN.csv`, `04b_DA_GCN_vs_GCT_paired.csv` (Supplementary Figure S14).*

**What it means.** Using the correct compositional method (CLR, not raw proportions), many genera differ — but per the batch caveat and the non-replication below, these are shown for completeness, not proposed as biomarkers.

### 5.4 Independent-cohort test — decisive
**Result.** Batch-clean cohort (PRJNA641258, 20 vs 20): **no diversity difference (p=0.25), no compositional difference (Bray R²=0.018, p=0.80), 0/344 genera significant.** A second cohort (PRJNA413125, Portugal) reproduces reduced diversity in carcinoma (Shannon p=0.004) but with oral taxa *depleted* — the "H. pylori paradox." *Sources: `validation_IT/`, `validation_PT/`.*

**What it means.** The tumour "oralization" biomarker **does not replicate** in clean data. The Portugal cohort proves the pipeline detects real signal when it exists (so the null is a true negative), and shows the oral-taxa direction is comparator-dependent. Conclusion: reduced richness is the only durable signal; the rest is batch + comparator choice. A genuinely cautionary contribution.

---

## PART 6 — SECONDARY ARM: MENDELIAN RANDOMISATION (null)

### 6.1 Causal test of 6 microbial exposures
**Result.** Two-sample MR of anti-H. pylori seropositivity + 5 gut genera on GC risk: **every IVW estimate null** (smallest p=0.35). Instruments strong (mean F 20.9–22.8, min F 19.3–20.3, all >10). *Source: `results/mr_real/MR_per_exposure_instruments_REAL.csv`.*

| Exposure | nSNP | IVW OR (95% CI) | p |
|---|---|---|---|
| Anti-H. pylori IgG | 17 | 0.96 (0.71–1.30) | 0.79 |
| Streptococcus | 15 | 1.10 (0.90–1.34) | 0.35 |
| Fusobacterium | 23 | 1.04 (0.79–1.36) | 0.79 |
| Prevotella | 15 | 0.98 (0.84–1.14) | 0.76 |
| Veillonella | 8 | 1.04 (0.84–1.29) | 0.69 |
| Lactobacillus | 10 | 0.96 (0.85–1.09) | 0.51 |

**What it means.** MR uses genetic variants as natural randomisation to test *causation* (not just correlation). No microbe shows a causal effect on GC risk. Instruments are strong (F>10), so this isn't weak-instrument failure — but see the power caveat.

### 6.2 Sensitivity diagnostics and MR-PRESSO
**Result.** MR-Egger intercepts null for 5/6 (Fusobacterium nominal, intercept 0.038 p=0.040); Cochran's Q null except H. pylori (p=0.035); **MR-PRESSO global** non-significant for 5/6, H. pylori borderline (RSSobs 32.6, p=0.052). East-Asian outcome (7,921 cases) sensitivity also null throughout. *Sources: `MR_PRESSO_global_REAL.csv`, `results/mr_real_eas/`.*

**What it means.** A full robustness suite — pleiotropy (Egger), heterogeneity (Q), outliers (PRESSO), directionality (Steiger), and a higher-powered ancestry check. All converge on the null. H. pylori is flagged as the one to read most cautiously (three diagnostics nudge together), but its estimate is null under every method.

### 6.3 The honest interpretation
**Result/caveat.** CIs (e.g. H. pylori 0.71–1.30) **span** the small effects other MR studies report (OR≈1.12). European GC outcome has only ~1,029 cases; gut (faecal) instruments are an imperfect proxy for the gastric-mucosal niche. *Source: §3.9 discussion.*

**What it means.** The null does not *disprove* a causal effect — the study is underpowered to confirm or exclude effects of the reported size, and the instruments are an imperfect proxy. So the microbiome–GC link is best regarded as **observational; causality neither established nor refuted.** This is the association-vs-causation discipline that is the paper's methodological signature.

---

## PART 7 — DRUG REPURPOSING (hypothesis-generating) + orthogonal support

### 7.1 Signature-reversal candidates
**Result.** Compounds reversing the tumour signature (both arms), ranked: resveratrol #1, PD-173074 #2, calcitriol #3, palbociclib #4, dovitinib #5, NVP-BEZ235 #6 — dominated by **PI3K/mTOR, CDK4/6 and FGFR inhibitors.** *Source: `results/drug_repurposing_integrated/candidate_drugs_ranked.csv`.*

**What it means.** Connectivity-map logic: find drugs whose expression footprint is the *opposite* of the tumour's. The hits are anti-proliferative classes — consistent with the proliferation-dominated tumour-up programme, but broad rather than GC-specific. Hypotheses, not leads.

### 7.2 Orthogonal genetic-dependency support (DepMap)
**Result.** In 35 gastric cell lines (DepMap CRISPR), the nominated target classes are genuine dependencies: **MTOR mean Chronos −1.18 (100% of lines dependent), CDK4 −0.82 (66%), PIK3CA −0.74 (57%)**. *Source: `results/depmap/gastric_dependency.csv`.*

**What it means.** An independent line of evidence: CRISPR knockout data show gastric cells actually *need* these genes to survive (Chronos < −0.5 = dependent). So the drug nominations aren't just expression-pattern coincidence — the targets are real vulnerabilities. Strengthens the hypothesis without overclaiming therapy.

---

## Summary — what the project actually shows

1. **One robust, externally-validated finding:** a stromal/cancer-associated-fibroblast co-expression programme underlies GC prognosis — preserved (Z 15.9–17.1) and prognostic in 3 external cohorts, localised to fibroblasts at single-cell resolution.
2. **A modest, honestly-bounded signature:** real (nested C 0.61) but adds negligible value over staging out-of-sample (ΔC +0.005) — a research instrument, not a clinical test.
3. **Two clean negatives, correctly framed:** the tumour microbiome "oralization" biomarker is batch-driven and doesn't replicate; MR finds no causal microbial effect (power-limited).
4. **The through-line:** every claim is scoped to its evidence; association is never dressed up as causation. That discipline is the paper's real contribution.

*Every number above is copied from the named result file. This document is the plain-language companion to `PAPER.md`.*

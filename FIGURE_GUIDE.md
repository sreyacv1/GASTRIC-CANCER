# FIGURE_GUIDE.md — every figure explained

For each figure: **why we made it, how we made it, the exact parameters, what
every term means, what we found, and what it does and does not prove.**

All numbers below are read from the committed result files, not from memory.
Each figure names the script that regenerates it.

---

## How to read the vocabulary (terms used across figures)

| Term | Plain meaning | Why it appears |
|---|---|---|
| **log2 fold change (log2FC)** | Doubling = +1, halving = −1 | Symmetric around zero, so a 2× rise and a 2× fall have equal visual weight |
| **adjusted P value** | P value corrected for testing ~21,000 genes | Without it, ~1,000 genes look "significant" by chance alone |
| **moderated t-statistic** | log2FC ÷ standard error, with the error shrunk toward the average across all genes | With only 36 normal samples, per-gene variance is unreliable; limma borrows information across genes |
| **z-score** | (value − row mean) ÷ row SD | Puts every gene on one colour scale; otherwise one high-expression gene dominates the heatmap |
| **NES** | Enrichment score scaled for gene-set size | Lets a 20-gene set and a 200-gene set be compared directly |
| **hazard ratio (HR)** | Relative rate of death per unit increase | HR 1.5 = 50% higher rate at any instant; not a probability |
| **C-index** | Probability the model ranks two patients' survival correctly | 0.5 = coin flip; 1.0 = perfect |
| **eigengene** | First principal component of a module | One number summarising a whole co-expression module per patient |
| **Cliff's δ** | Non-parametric effect size, −1 to +1 | Works on skewed count data where a difference in means would mislead |
| **R² (PERMANOVA)** | Fraction of community variation explained | 0.065 = the grouping explains 6.5% of microbiome differences |

---

## Figure 1 — Transcriptional programmes
**Script:** `make_figures.R` · **Source:** `results/enrichment/GSEA_Hallmark_*.csv`

**Why we made it.** Before any prediction modelling, the reader needs to know
that the tumours behave like gastric cancer. This is the orientation figure.

**How.** Genes were ranked by moderated t-statistic and tested against the 50
MSigDB Hallmark sets with fgsea (`minSize = 10`, `maxSize = 500`, `eps = 0`).
Gene sets, not single genes, are used because single-gene lists are unstable
across cohorts while pathway-level signals replicate.

**What we found.** Proliferation up (E2F targets NES 3.75, G2M 3.64, MYC 2.67);
energy metabolism down (oxidative phosphorylation −2.54, fatty-acid metabolism
−2.58). That pairing is the **Warburg effect** — tumours shift from efficient
mitochondrial respiration to fast aerobic glycolysis, trading ATP yield for
speed and biosynthetic precursors.

Panel (b) is the more interesting result: diffuse-type tumours are dominated by
**EMT** (NES 3.17, adjusted P = 1.6 × 10⁻⁴⁰). Epithelial cells lose adhesion and
acquire migratory mesenchymal character — the mechanistic explanation for why
diffuse tumours infiltrate rather than form discrete masses, and why they carry
worse prognosis.

**Limitation.** GSEA shows coordinated expression change, not causation or
protein-level activity.

---

## Figure 2 — The 25-gene signature, evaluated honestly
**Script:** `make_figures.R` · **Source:** `results/nested_cv/`, `results/meta/`

**Why we made it.** This is the figure most likely to be attacked in review, so
it is built to pre-empt the attack.

**How.** LASSO-penalised Cox regression (`alpha = 1`, `family = "cox"`,
`set.seed(1105)`). LASSO shrinks most coefficients to exactly zero, performing
selection and fitting in one step — necessary when genes outnumber patients.

**The critical design choice.** Discrimination is estimated by **20 repeats of
5-fold nested cross-validation** (`N_REPEAT = 20`, `K_OUTER = 5`, `K_INNER = 5`).
Gene selection happens *inside* the training folds only. If genes are chosen on
all patients and then evaluated on the same patients, the model has already seen
the answers.

**What we found — and why it matters.**

| Estimate | C-index |
|---|---|
| Apparent (selection sees evaluation data) | 0.72 |
| **Nested CV (honest)** | **0.611** (95% CI 0.562–0.659) |
| Uno C | 0.573 |

The ~0.11 gap *is* the finding. Many published signatures report the 0.72-style
number. For calibration, we re-ran a **previously published 5-gene gastric
signature** on identical data with identical correction: **C = 0.481 (value from R 4.3.3/`rms` 6.8-1; re-running under `rms` 8.1-1 gives 0.52 — both at chance, apparent C 0.54 unchanged) — below
chance.**

Panel (b): pooled HR 1.19 (95% CI 0.96–1.47) — **not significant**. We report
this rather than hide it. Panel (d): the signature works early (HR 1.49 at 12
months) and decays to null by 36–60 months, so proportional hazards fails and a
single HR would misrepresent it.

---

## Figure 3 — Stromal module: preservation, prognosis, localisation *(the paper's centre)*
**Script:** `make_figures.R` (lines 195-232) · **Source:** `results/module_preservation/preservation_summary_RED.csv`, `results/module_preservation/module_eigengene_cox_external.csv`, `results/scrna/gene_dominant_celltype.csv`

**Panels:** (a) preservation Zsummary bars, (b) eigengene Cox forest, (c) single-cell localisation of hub genes. The WGCNA network *construction* (dendrogram, module-trait heatmap, soft-power selection) is Figure 5; Figure 3 shows what the module *does*.

**Why.** A 25-gene list is fragile. A **module** — a whole co-expressed
programme — is more robust, because losing any one gene barely moves the
eigengene.

**How (network built once, shown in Figure 5).** `blockwiseModules(networkType = "signed hybrid", TOMType = "signed",
corType = "bicor", maxPOutliers = 0.1, deepSplit = 2, minModuleSize = 30,
mergeCutHeight = 0.25, maxBlockSize = 6000, randomSeed = 1105)`.

- **signed** — genes moving *oppositely* are not treated as similar
- **bicor** — biweight midcorrelation, resistant to outlier samples
- **minModuleSize = 30** — smaller modules are noise
- **mergeCutHeight = 0.25** — modules with eigengene correlation > 0.75 are merged

**What we found.** A stromal/fibroblast module, preserved across all three
external cohorts (Zsummary 15.9, 16.8, 17.1 — Z > 10 is strong preservation) and
prognostic in each.

**The key result:** the module is prognostic in **GSE84437 (HR 1.24/SD,
P = 0.0021) — the cohort where the 25-gene signature completely failed**
(HR 1.11, P = 0.46). This is why the paper's title leads with the stromal
programme, not the signature.

**Biology.** Cancer-associated fibroblasts build the collagen-rich matrix that
stiffens tissue, shields tumour cells from drugs, and secretes pro-invasive
signals. Consistent with Figure 1's EMT result.

---

## Figure 4 — Immune deconvolution and the microbiome caution
**Script:** `make_figures.R` · **Source:** `results/immune/`, `results/microbiome_biomarker/`

**How.** MCP-counter and xCell estimate cell-type abundance from bulk expression
using cell-type-specific marker genes.

**Validation is the point.** We checked estimates against **pathologist-measured**
values:

| Comparison | ρ | n | P |
|---|---|---|---|
| T cells vs measured leukocyte % | **0.666** | 272 | 3.6 × 10⁻³⁶ |
| CD8 T cells vs measured leukocyte % | 0.468 | 272 | 3.1 × 10⁻¹⁶ |
| CD8 T cells vs measured **lymphocyte** % | **0.020** | 274 | 0.74 |

Read honestly: deconvolution tracks **overall immune burden** but *not*
lymphocyte subset composition. We state this rather than quoting only ρ = 0.67.

**Panel (c) — the negative result we chose to publish.** In the Japanese cohort,
tumour-vs-control microbiome separation looked real (Bray–Curtis R² = 0.065,
P = 0.001) — but tumour samples sat on **separate sequencing flowcells**. After
flowcell adjustment, **R² collapses to 0.011**. A random forest predicts
*flowcell* at 78% (baseline 55%), and the top discriminating genera —
*Cutibacterium*, *Sphingomonas*, *Methylobacterium* — are **known reagent
contaminants**, not gastric flora.

Independent cohorts: Italy null (R² = 0.018, P = 0.80); Portugal replicates
reduced diversity (Shannon P = 0.004).

---

## Figure 5 — Network construction quality
**Script:** `analysis/24_main_figures_5to8.R`

Panel (c) justifies the **soft-thresholding power**, which raises correlations to
a power so strong correlations dominate and weak ones vanish — approximating the
scale-free topology of real biological networks, where a few hub genes have many
connections. Model fit R² ranged 0.865–0.886 across candidate powers; we chose
the lowest power at which fit plateaus, avoiding an over-sparse network.

---

## Figure 6 — Differential expression and its reproducibility
**Script:** `analysis/24_main_figures_5to8.R`

**Thresholds:** |log2FC| > 1 and adjusted P < 0.05 → **2,134 up, 2,362 down of
21,446 tested features**.

**On "21,446".** These are *features*, not protein-coding genes. TCGA quantifies
against GENCODE (60,660 features). After filtering, 21,446 remain: **16,164
protein-coding (75.4%)**, 1,567 ncRNA, 572 pseudogenes, 2,944 unmapped. So the
tested set contains *fewer* protein-coding genes than the ~19,000–20,000 human
total, not more.

**Batch robustness (asked by every reviewer).** Refitting with tissue-source-site
as a batch covariate: **4,129 of 4,456 DEGs retained (92.7%)**, logFC correlation
**r = 0.974**. Caveat stated: site and tumour status are partly confounded, so
adjustment cannot be perfect.

**Panel (c) — an important honesty point.** The left panel compares the
integrated ranking against TCGA alone. Since TCGA supplies all 412 tumours of the
integrated contrast, that r = 0.73 is **internal consistency, not replication**.
Only GSE27342 (r = 0.62) and GSE63089 (r = 0.58) are independent cohorts.

---

## Figure 7 — Signature coefficients and validation
**Script:** `analysis/24_main_figures_5to8.R`

25 genes, 16 positive and 9 negative coefficients. Validation: ACRG C = 0.608
(HR 1.90, 95% CI 1.37–2.62, log-rank P = 8.6 × 10⁻⁵); GSE15459 C = 0.575;
GSE84437 C = 0.530 (**not significant**, P = 0.46).

**Stability, reported honestly.** Across 200 bootstrap resamples, **13 of 25 genes
were selected in > 50% of resamples but none exceeded 80%** (median 0.505; most
stable NETO2 0.77, EGF 0.75, SRMS 0.73). Individual gene identities are not
stable — which is precisely why the module (Figure 3) is the stronger claim.

---

## Figure 8 — Mendelian randomisation
**Script:** `analysis/24_main_figures_5to8.R`

**Why.** Observational microbiome associations cannot separate cause from
consequence — a tumour changes its own microenvironment. MR uses genetic variants
as instruments: alleles are randomised at conception, so they cannot be caused by
the disease.

**Parameters:** instruments at P < 5 × 10⁻⁸ (relaxed to 1 × 10⁻⁵ where too few),
clumped at r² = 0.001 within 10,000 kb, harmonised with `action = 2`.
Instrument strength F = 19.3–20.3 minimum (F > 10 is the conventional
weak-instrument threshold).

**Result: no causal effect for any of six exposures**, smallest IVW P = 0.35.
MR-PRESSO global tests non-significant for all six (P = 0.052–0.676).

**Stated limitation.** With 8–23 instruments per exposure, power is low. This is
a **null result, not evidence of absence**.

---

## What makes this project publishable

1. **The honest-evaluation result.** Nested CV C = 0.611 versus apparent 0.72, with a published competitor scoring 0.481 (value from R 4.3.3/`rms` 6.8-1; re-running under `rms` 8.1-1 gives 0.52 — both at chance, apparent C 0.54 unchanged) — below chance — on identical data.
2. **The stromal module.** Externally validated across three cohorts, prognostic where the sparse signature failed.
3. **A properly diagnosed batch artefact.** Most papers would have published R² = 0.065 as a microbiome finding. We traced it to flowcell structure and reagent contaminants.
4. **Every negative reported.** Null MR, non-significant pooled HR, failed GSE84437 validation, unstable gene selection, deconvolution that fails on lymphocyte subsets.

The unifying claim is not that we built a great predictor. It is that **the
stromal compartment carries reproducible prognostic signal**, and that several
attractive-looking results do not survive honest evaluation.

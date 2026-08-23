# The Project, A to Z

**A multi-omics analysis of gastric cancer identifies an externally-validated
stromal/fibroblast prognostic program, with cautionary tumor-microbiome and
Mendelian-randomization assessments**

This is the master reference: what the project asks, what it did, what it found,
what it does *not* claim, and where every number lives.

**How to trust the numbers here.** Every quantitative claim below is checked by
`analysis/verify_a_to_z.py`, which re-parses the value from its source file under
`results/` (or, for cohort composition, the raw GEO series matrix) and confirms the
value as written appears in this document — in the Markdown, Word and PDF renderings
simultaneously. The script prints the number of checks and exits non-zero on any
mismatch, so you can re-run it yourself rather than take this paragraph on trust. Where
a number is quoted, the file it came from is named.

Four companion documents go deeper on specific angles and are not repeated here:

| Document | What it is for |
|---|---|
| `LEARN_MY_PROJECT.md` | Teaching from first principles — DNA, genes, expression, cancer biology |
| `PROJECT_END_TO_END.md` | Study design and measurement: why paired samples, what pairing does and does not fix |
| `PROJECT_WALKTHROUGH.md` | Statistical reasoning per analysis, with the traps named |
| `RESULTS_COMPENDIUM.md` | Every result with a plain-language explanation |

---

# PART 1 — THE QUESTION

## 1.1 What the project asks

Two questions, one primary and one secondary.

**Primary (transcriptomic).** In gastric cancer, which genes are expressed
differently in tumor versus normal stomach; do those genes organize into coherent
programs; and does any program carry information about how long a patient survives
that is not already contained in standard clinical staging?

**Secondary (microbial).** Is the community of bacteria living in gastric tissue
different in tumors, and if so, is there evidence that any specific microbe *causes*
gastric cancer rather than merely accompanying it?

## 1.2 Why the second question needs a different kind of evidence

An association between a microbe and a tumor has at least four explanations:

1. the microbe causes the cancer;
2. the cancer changes the local environment so the microbe thrives (reverse causation);
3. something else — diet, smoking, age — drives both (confounding);
4. the two groups of samples were processed differently and the "difference" is technical (batch effect).

Sequencing tissue can never distinguish these on its own. That is precisely why the
project adds Mendelian randomization, which uses inherited genetic variants — fixed at
conception, therefore not caused by a tumor that appears decades later — as instruments
for microbial exposure. The design is what licenses causal language, and the result
(Part 5) is a clean null.

---

# PART 2 — THE DATA

| Arm | Dataset | N | Role |
|---|---|---|---|
| Transcriptome (discovery) | TCGA-STAD | 412 tumor / 36 adjacent normal | Differential expression, WGCNA, signature training |
| Transcriptome (normals) | GTEx v10 stomach | 407 | Adds normal contrast for the integrated analysis |
| Survival validation | GSE62254 (ACRG) | 300 (152 events) | Independent survival cohort |
| Survival validation | GSE15459 | 191 (95 events) | Independent survival cohort |
| Survival validation | GSE84437 | 431 (207 events) | Independent survival cohort |
| DE replication | GSE27342, GSE63089 | — | Independent tumor-vs-normal replication |
| Single cell | GSE134520 | 43,992 cells | Which cell type expresses the module |
| Tissue microbiome | DDBJ PRJDB20660 | 944 samples | 16S discovery, raw FASTQ reprocessed |
| Microbiome validation | PRJNA641258, PRJNA413125 | — | Independent 16S cohorts (Italy, Portugal) |
| GWAS (MR) | IEU OpenGWAS | 6 exposures, 1,029 cases | Causal inference on microbial exposures |

Survival-cohort sample sizes and event counts above are those used in the analyses and
match `results/module_preservation/module_eigengene_cox_external.csv` and
`results/meta_HK/meta_inputs.csv`.

**One design fact worth stating plainly.** The transcriptome cohort and the microbiome
cohort share **zero patients** — 415 TCGA patients and 944 microbiome samples with no
overlapping identifiers. This is why the project is *multi-omics* in the sense of two
parallel arms answering different questions, and why genuine joint factor analysis
(MOFA) was assessed and correctly declined: with no shared samples there is nothing to
factorize jointly.

---

# PART 3 — THE TRANSCRIPTOMIC CORE

## 3.1 Differential expression

Of **21,446** genes tested in TCGA-STAD, **2,134 are up-regulated** and **2,362
down-regulated** in tumor versus adjacent normal at |log₂FC| > 1 and BH-adjusted
p < 0.05 (`results/tables/TCGA_DEG_results.csv`).

The up-regulated set is dominated by proliferation machinery (MKI67, TPX2, ECT2,
KIF14, CENPF, TOP2A) and stromal/matrix genes (COL1A1). The down-regulated set is
differentiated gastric tissue: parietal-cell and metabolic genes, AQP4, ADH7, KRT24,
CLEC3B, MYOC. Biologically this is the expected signature of a tissue losing its
specialized identity and acquiring proliferative, matrix-rich character.

**Why 21,446 and not ~20,000.** The raw matrix carries 60,660 GENCODE features. After
filtering, 21,446 remain, of which 16,164 (75.4%) are protein-coding; the rest are
long non-coding RNAs, pseudogenes and unmapped identifiers. The count is not an error —
it simply is not a count of protein-coding genes.

**Batch sensitivity was tested, not assumed.** Adding tissue source site as a covariate
changes the DEG list from 4,456 to 4,569 genes with 4,129 shared (92.7%), and log-fold
changes correlate at r = 0.974. The result does not depend on the batch model.

## 3.2 The stromal co-expression module (the central finding)

WGCNA on the 5,000 most variable genes at soft power 3 (scale-free R² = 0.88) detected
ten modules. The **red** module — **263 genes**
(`results/wgcna_real/hub_genes_prognostic_module.csv`) — is a
cancer-associated-fibroblast / extracellular-matrix program, and it is the project's
most robust result for one reason: **it replicates externally**.

Module preservation (`modulePreservation`, 200 permutations) gives
Zsummary = **15.9** (ACRG), **16.8** (GSE15459), **17.1** (GSE84437) — all far above the
Z > 10 "strong preservation" threshold.

**Three different gene counts appear around this module; keep them distinct.** The
module itself is 263 genes. Preservation is computed only on the module genes that the
external platform actually measures, which is **77** of 263 on the ACRG array and
**233** on the other two (the `moduleSize` column of
`results/module_preservation/preservation_stats_*.csv` is this per-cohort overlap, not
the module size). The single-cell localization test in §3.3 uses a further subset — the
23 non-circular hub genes. That the structure still scores Zsummary 15.9 in ACRG on
fewer than a third of its genes is evidence of robustness, not of a smaller module.

The red module's eigengene is prognostic in all three independent cohorts, always in
the same direction (more stromal activation → worse survival):

| Cohort | n (events) | HR per SD | 95% CI | p |
|---|---|---|---|---|
| ACRG/GSE62254 | 300 (152) | 1.274 | 1.091–1.487 | 0.00216 |
| GSE15459 | 191 (95) | 1.548 | 1.248–1.919 | 6.8×10⁻⁵ |
| GSE84437 | 431 (207) | 1.237 | 1.080–1.416 | 0.00211 |

Source: `results/module_preservation/module_eigengene_cox_external.csv`.

**The honest caveat, stated in the paper.** Within TCGA the module eigengene is
prognostic unadjusted (HR/SD 1.307, 95% CI 1.116–1.532, p = 0.00093) but **loses
significance after adjusting for stage and age** (HR 1.347, 95% CI 0.954–1.901,
p = 0.090; `results/wgcna_real/ME_survival_cox_adjusted.csv`). The program tracks tumor
stage and is not cleanly separable from it. The finding is that the program is real,
reproducible and prognostic — not that it is independent of staging.

**It does not depend on an arbitrary parameter choice.** Rebuilding the network at
powers 3, 6, 9 and 12 keeps the module survival-associated throughout
(HR/SD 1.283–1.308, all p < 0.0024) with 11–12 of 12 hub genes staying inside it.

## 3.3 Which cell expresses it

Single-cell data answers the "which cell" question that bulk expression cannot. Across
43,992 cells in eight annotated types (`results/scrna/celltype_composition.csv`), **all 23 of 23** hub genes in the non-circular
test are fibroblast-dominant, with the dominant-type fraction ranging 0.527–0.998
(median 0.96; `results/scrna/gene_dominant_celltype_noncircular.csv`).

**Why "non-circular" matters.** If you select genes *because* they are stromal and then
show they are stromal, you have proven nothing. The non-circular test removes the genes
whose selection depended on stromal identity and asks whether the remainder still
localize to fibroblasts. They do. This is the difference between a real result and a
tautology, and it is worth understanding as a general lesson.

---

# PART 4 — THE PROGNOSTIC SIGNATURE, REPORTED HONESTLY

## 4.1 The number that matters

A 25-gene LASSO–Cox signature separates survival in TCGA with an **apparent C-index of
0.72**. That figure is optimistic and the paper says so, because gene screening and
LASSO ran once on the whole cohort — the model saw the data it was then scored on.

Under **20 repeats of 5-fold nested cross-validation**, with standardization, screening
and tuning rebuilt inside every training fold and performance taken only from
untouched out-of-fold predictions, the honest discrimination is:

| Metric | Value |
|---|---|
| Harrell C (ensemble) | **0.6112** |
| Uno C (ensemble) | 0.5727 |
| Harrell C (per-repeat) | 0.5981 |
| Integrated Brier score | 0.1789 |

Source: `results/nested_cv/performance.csv`.

**The optimism is ~0.11 of C-index.** Reporting 0.72 would not be fabrication — it is a
real quantity — but presenting it as generalizable performance would be misleading, and
that gap is the single most transferable statistical lesson in the project.

## 4.2 External validation: two of three

| Cohort | n (events) | C-index | Log-rank p | HR high-vs-low (95% CI) |
|---|---|---|---|---|
| TCGA (training) | 383 (156) | 0.719 | 2.05×10⁻¹² | — |
| ACRG/GSE62254 | 300 (152) | 0.608 | 8.6×10⁻⁵ | 1.90 (1.37–2.62) |
| GSE15459 | 191 (95) | 0.575 | 0.014 | 1.68 (1.11–2.54) |
| GSE84437 | 431 (207) | 0.530 | 0.46 | 1.11 (0.84–1.46) |

Sources: `results/validation/cindex_comparison.csv`,
`results/validation_multi/cindex_HR_summary.csv`.

**GSE84437 is a genuine negative, and the project treats it as one.** The proposed
explanation is range restriction: GSE84437 is 88.7% pT3–T4 disease (384 of 433 staged
samples: T1 11, T2 38, T3 92, T4 292; `data/geo/GSE84437_series_matrix.txt.gz`), and a stromal
signature discriminates by capturing *variation* in desmoplastic content, which
increases with invasion depth. That hypothesis was then **tested directly** by
stratifying on pT stage — and it failed: C-index stayed below 0.5 in every stratum.
The paper reports this as a tested negative rather than an explained-away artifact.
Testing your own escape hatch and reporting that it failed is exactly the behavior
reviewers look for.

**The pooled effect is not significant.** A per-1-SD, age/stage-adjusted random-effects
meta-analysis (REML with Hartung–Knapp) across the three external cohorts gives a
pooled **HR 1.188 (95% CI 0.962–1.466, p = 0.073)**, I² = 19.2%
(`results/meta_HK/meta_result.csv`). Directionally consistent, modest, and short of
significance under a conservative model.

## 4.3 Time-varying and predictive structure

In ACRG the signature is prognostic after adjustment for stage and age (HR 1.76,
p = 7.4×10⁻⁴), but the proportional-hazards assumption is **violated**
(cox.zph p = 0.0018; `results/timevarying_ACRG/coxzph.csv`). A time-varying model shows
the effect is early (HR/SD 1.49 at 12 months) and attenuates to null later (1.03 at 36
months, 0.87 at 60). It marks earlier mortality, not a durable gradient — a distinction
most papers with this data pattern silently skip.

The hazard is also modified by stage (risk × stage LRT p = 0.044), concentrated in
advanced and mesenchymal disease; the risk × molecular-subtype interaction is *not*
significant (p = 0.094; `results/predictive/interaction_tests.csv`).

## 4.4 The benchmark against a published signature

A previously published 5-gene gastric-cancer signature, re-tested on identical data
under the same leakage-controlled pipeline, scores apparent C = 0.5446 and
optimism-corrected **C = 0.5193** — chance level, with the risk direction inverted
(HR high-vs-low 0.627).

**One subtlety the project documents rather than hides.** The bootstrap optimism
correction is implementation-sensitive: the same script and seed give 0.4807 under
R 4.3.3 with `rms` 6.8-1, and 0.5193 under R 4.5.3 with `rms` 8.1-1. Both are at chance
(0.50 is the no-information value) and the apparent C is identical in both, so the
conclusion is unchanged — but both environments are archived
(`results/base_paper_replication/sessionInfo.txt` for R 4.5.3 / `rms` 8.1-1 and
`sessionInfo_rms6.8-1_R4.3.3.txt` for R 4.3.3 / `rms` 6.8-1) so the discrepancy is
inspectable rather than mysterious.

---

# PART 5 — THE MICROBIAL ARM (CAUTIONARY, AND THAT IS THE POINT)

## 5.1 What the 16S data shows

Differential abundance on CLR-transformed counts finds **44 of 61** genera differing
between control and cancer-adjacent mucosa, and **18 of 61** in the paired
cancer-adjacent versus tumor comparison
(`results/microbiome_biomarker/04a_*.csv`, `04b_*.csv`).

Alpha diversity declines from control to cancer-adjacent mucosa, with the effect size
reported rather than only a p-value:

| Metric | Median control | Median cancer-adjacent | Cliff's δ | Wilcoxon p |
|---|---|---|---|---|
| Observed richness | 47 | 30 | −0.27 | 3.3×10⁻⁷ |
| Shannon | 2.176 | 1.923 | −0.122 | 0.021 |
| Simpson | 0.781 | 0.756 | −0.072 | 0.172 |

Note the internal discordance, which the paper discloses: **Simpson diversity does not
significantly decline**. Richness falls; evenness-weighted diversity does not clearly
follow. Reporting the metric that disagrees is not a weakness — it is what stops a
reviewer from finding it for you.

## 5.2 The batch problem, diagnosed rather than ignored

A random-forest classifier separates cancer from control at **AUC 0.916** (95% CI
0.896–0.936). The same approach predicts **sequencing flowcell** at 77.6% accuracy
against a 54.5% majority baseline
(`results/microbiome_biomarker/05_rf_metrics_and_batch_sanity.csv`).

That second number is the honest finding. Tumor samples were sequenced on separate
flowcells, so a classifier that appears to detect cancer may be detecting the run.
Compounding it, tumor tissue is the **richest** group in this dataset — an inversion of
the expected biology that is itself evidence of a technical artifact rather than a
discovery. This is why the arm is framed as *cautionary*: the analysis is sound, and
what it establishes is that the apparent signal cannot be cleanly separated from batch.

## 5.3 Mendelian randomization: a clean null

Two-sample MR of six microbial exposures against ancestry-matched European gastric
cancer (1,029 cases) finds **no significant causal effect anywhere**:

| Exposure | nSNP | IVW OR (95% CI) | p |
|---|---|---|---|
| Anti-*H. pylori* IgG seropositivity | 17 | 0.96 (0.71–1.30) | 0.79 |
| *Streptococcus* (genus) | 15 | 1.10 (0.90–1.34) | 0.35 |
| *Fusobacterium* A | 23 | 1.04 (0.79–1.36) | 0.79 |
| *Prevotella 9* | 15 | 0.98 (0.84–1.14) | 0.76 |
| *Veillonella* | 8 | 1.04 (0.84–1.29) | 0.69 |
| *Lactobacillus* | 10 | 0.96 (0.84–1.09) | 0.51 |

Source: `results/mr_real/MR_per_exposure_instruments_REAL.csv`.

**Read this as "not established at this power," not "disproved."** No exposure reached
three genome-wide-significant SNPs, so all six use a locus-wide-suggestive threshold of
p < 1×10⁻⁵ — a choice that biases toward the null through weak-instrument effects,
although the per-SNP F-statistics (minimum 19.3) show the retained instruments are not
themselves weak. With 1,029 cases the study is underpowered for modest causal effects.
A null under acknowledged low power is a different claim from a null under adequate
power, and conflating them is a common error.

Sensitivity diagnostics ran properly: MR-Egger intercepts, Cochran's Q, weighted-median
and mode estimators, Steiger filtering, leave-one-out, and MR-PRESSO global tests
(all six global p > 0.05, lowest 0.052 for *H. pylori*).

---

# PART 6 — IMMUNE AND DRUG-REPURPOSING ARMS

## 6.1 Deconvolution validated against pathology

Rather than trusting a deconvolution algorithm, its output was checked against
pathologist-scored histology (`results/immune/validation_vs_measured.csv`):

| Estimate | Measured against | n | Spearman ρ | p |
|---|---|---|---|---|
| T cells (MCP) | Leukocyte % | 272 | 0.666 | 3.6×10⁻³⁶ |
| ImmuneScore (xCell) | Leukocyte % | 272 | 0.652 | 2.2×10⁻³⁴ |
| CD8 T cells (MCP) | Leukocyte % | 272 | 0.468 | 3.1×10⁻¹⁶ |
| CD8 T cells (MCP) | Lymphocyte infiltration % | 274 | 0.020 | 0.744 |

**Validation is measure-specific, and that is the finding.** The estimates track overall
immune burden well and lymphocyte *subset composition* not at all. So the deconvolution
should be read as capturing how much immune infiltrate is present, not which subsets
compose it. Tumors are enriched for the macrophage/monocyte compartment with no net
CD8⁺ gain, and the CD8 score is **not prognostic** here (Cox HR 1.042, 95% CI
0.944–1.150, p = 0.411; `results/immune/CD8_survival_summary.csv`) — reported as
observed rather than dropped for being null.

## 6.2 Drug repurposing (hypothesis-generating only)

Signature-reversal enrichment against LINCS/CMap and DSigDB nominates CDK4/6, PI3K/mTOR
and FGFR inhibitors — palbociclib, NVP-BEZ235, dovitinib, PD-173074 — with orthogonal
support from DepMap dependency scores (PIK3CA −0.742, MTOR −1.184, CDK4 −0.825,
CDK6 −0.548; `results/depmap/gastric_dependency.csv`, 35 gastric lines). This arm is explicitly hypothesis-generating: no compound was tested in a
cell line or animal in this project, and the paper claims nothing beyond a ranked
candidate list.

---

# PART 7 — WHAT THE PROJECT DOES NOT CLAIM

Stating this precisely is what makes the rest credible.

1. **Not a novel biological discovery.** The stroma/CAF–prognosis link is well
   established. The contribution is external replication and honest quantification.
2. **Not a clinically usable test.** Honest C ≈ 0.61 is modest, the pooled external
   effect is non-significant, and the module is not separable from stage.
3. **No causal microbial claim.** MR is null and underpowered.
4. **No demonstrated therapeutic activity.** Drug candidates are computational.
5. **Not a joint multi-omics factorization.** The two arms share no patients.

---

# PART 8 — REPRODUCIBILITY

- Every reported number maps to a named script in `analysis/` and an output file under
  `results/` (see `PIPELINE.md` for the full mapping). Script and commit counts are
  deliberately not quoted here — they change with every commit, so any number written
  into this document would be stale by the time you read it. Run `git rev-list --count
  HEAD` and `ls analysis/ | wc -l` for the current values.
- Every result table and figure is committed alongside the code that produced it.
- **Claim verification is executable.** `analysis/verify_a_to_z.py` re-parses each
  quantitative claim in this document from its source file and confirms the value
  appears in the Markdown, Word and PDF renderings simultaneously. Run it from the
  repository root; it exits non-zero on any mismatch.
- **17 figures** (8 main + 9 supplementary), each verified for edge clipping, aspect
  ratio and text legibility at print size.
- Pinned package versions and `sessionInfo()` dumps archived, including both R/`rms`
  environments for the version-sensitive benchmark.
- Known gap, stated: two scripts (`17_external_utility_ACRG.R`,
  `09_functional_enrichment.R`) depend on packages unavailable in the current
  environment (`dcurves`, `clusterProfiler`), so their committed figures cannot be
  re-derived here. Their outputs are committed; independent re-execution is unverified.
- Known gap, stated: the harmonized MR SNP data was never persisted, so Figure 8 and
  Supplementary S8 panels cannot be re-plotted without a live OpenGWAS token.

## Outstanding author actions

1. Author names, ORCIDs, corresponding author and CRediT contributions — marked
   placeholders in the manuscript.
2. Repository URL and a minted DOI (Zenodo or equivalent) at submission.
3. Optional: five of the 17 embedded figures sit above the reference paper's aspect
   range (height/width > 1.00) and would need re-layout to match — Figure 4 (`Fig7.png`)
   at 1.38, Supplementary S2 at 1.14, S8 at 1.13, S6 at 1.07 and S7 at 1.04; measured
   directly from the embedded PNGs.

---

# PART 9 — THE ONE-PARAGRAPH VERSION

Gastric tumors lose differentiated gastric identity and gain proliferative and
matrix-rich programs. One of those programs — a fibroblast/extracellular-matrix module —
is preserved in three independent cohorts (Zsummary 15.9–17.1) and prognostic in all
three, and single-cell data localizes all 23 of its testable hub genes to fibroblasts.
A 25-gene signature built from this biology achieves an honest cross-validated C-index
of 0.611 (against an apparent 0.72), validates in two of three external cohorts, and
does not reach significance when pooled — while a published competitor scores at chance
on identical data. In parallel, gastric tissue microbiota differ between tumor and
normal, but a classifier predicts sequencing flowcell at 77.6% and tumor tissue is
implausibly the richest group, so the signal cannot be separated from batch; and
Mendelian randomization finds no causal effect for any of six exposures at the available
power. The project's contribution is a replicated stromal prognostic program, reported
with its optimism quantified and its negative results intact.

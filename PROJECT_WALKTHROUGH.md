# Your Project, A to Z — Biology, Statistics, and Technical Pipeline

A complete teaching walkthrough of the gastric-cancer multi-omics study. Read it top to
bottom once; afterwards each section stands alone as a reference. Three lenses run through
every section: **BIOLOGY** (what it means in a tumour), **STATS** (why the method is valid),
**TECH** (how the code does it).

---

## 0. The one-paragraph version

You asked whether the biology of gastric-cancer prognosis can be read out of public
transcriptomes, and whether two fashionable "add-ons" — the tumour microbiome and a causal
(Mendelian-randomisation) claim about microbes — hold up when tested rigorously. Your answer
is honest and three-part: (1) **Yes for the biology** — prognosis is written in the tumour
**stroma / cancer-associated fibroblasts (CAFs)**, a signal you find four independent ways and
validate in external cohorts. (2) **A signature can be built but it is not a clinical tool** —
it adds essentially nothing over TNM stage out-of-sample. (3) **The microbiome and MR arms are
negative/uninterpretable once you control for confounding** — and you say so, instead of
dressing up noise as discovery. The scientific contribution is as much *methodological honesty*
(separating association from causation, batch from biology) as it is the CAF finding itself.

---

## 1. The biology you need first

**Gastric cancer (GC).** Cancer of the stomach lining. Globally a top-5 cause of cancer death,
especially in East Asia. Two big ideas frame your study:

- **Lauren classification** (ref 20): GC splits histologically into **intestinal** type
  (gland-forming, arises through a chronic-inflammation → atrophy → metaplasia → dysplasia
  cascade, often *H. pylori*-driven) and **diffuse** type (poorly cohesive, scattered signet-ring
  cells, worse prognosis, strong EMT). Your data reproduce the textbook: **EMT is enriched
  specifically in the diffuse subtype** (§3.2), and your CAF score is highest in diffuse tumours
  (`purity/lauren_association.csv`, p=1.9×10⁻⁸).

- **TNM / AJCC stage.** The clinical gold standard for prognosis: how deep the tumour invades
  (T), lymph-node spread (N), metastasis (M). Everything your signature does has to be judged
  *against stage*, because stage is what an oncologist already has. This is why your honest
  finding — "adds little beyond stage" — is the scientifically correct headline.

**The tumour microenvironment (TME).** A tumour is not just cancer cells. It is cancer cells
plus immune cells plus blood vessels plus **fibroblasts embedded in extracellular matrix (ECM)**.
In many carcinomas the fibroblasts become **cancer-associated fibroblasts (CAFs)** — activated,
matrix-depositing, signalling cells that promote invasion and correlate with poor outcome. The
markers you keep seeing — **COL1A1, COL1A2, COL3A1, DCN, LUM, POSTN, FAP, BGN, FN1, SPARC** — are
collagens and matrix proteins: the molecular fingerprint of desmoplastic (fibrous) stroma.

**Why "stromal prognosis" makes biological sense.** A tumour that has recruited abundant
activated stroma is one that is invading and remodelling tissue — a more aggressive tumour. So a
bulk-tissue gene-expression signal dominated by CAF/ECM genes is, in effect, measuring how much
desmoplastic stroma is in the biopsy, which tracks invasion and therefore prognosis. The subtlety
your paper is careful about: this also means the signal is partly **stromal admixture** (how much
stroma happened to be in the sampled piece), which is why it is not cleanly independent of stage.

---

## 2. The data (TECH + why each dataset is there)

| Dataset | Modality | Role | Why |
|---|---|---|---|
| **TCGA-STAD** | bulk RNA-seq, ~412 tumours | discovery | the reference GC transcriptome + clinical/survival |
| **GTEx** v10 | bulk RNA-seq, normal stomach | normal reference for DE | TCGA has few normals; GTEx supplies healthy tissue |
| **GSE62254 (ACRG)** | microarray, 300 | external validation | the definitive Asian GC cohort with rich clinical data |
| **GSE15459, GSE84437** | microarray | external validation | independent survival cohorts (different platforms) |
| **GSE27342, GSE63089** | microarray | DE replication | independent tumour-vs-normal contrasts |
| **GSE134520** | single-cell RNA-seq | cell-type localisation | resolves *which cell* expresses the module |
| **PRJDB20660** | 16S rRNA (tissue) | microbiome discovery | 944 libraries, reprocessed from raw FASTQ |
| **PRJNA641258 / PRJNA413125** | 16S | microbiome validation | independent cohorts to test replication |
| **IEU OpenGWAS** | GWAS summary stats | Mendelian randomisation | genetic instruments for microbes + GC outcome |

**Key technical decision — combining TCGA + GTEx.** Tumour vs normal DE needs normals. Mixing two
platforms/labs creates a **batch effect** so severe it dominates real biology. You harmonise with
**ComBat** and scale-standardise, but — crucially — you *diagnose that it is not fully fixed*
(next section). That honesty is the backbone of the whole paper.

---

## 3. Differential expression and the λ=17 problem (STATS — the most important lesson)

**What DE does.** For each of 12,899 genes, test whether mean expression differs between
tumour (n=412) and normal (n=443), using **limma** moderated-t statistics (ref 34). Output:
log fold-change + p-value + BH-adjusted p per gene. You find 3,722 up / 4,025 down.

**The λ (genomic inflation factor).** Here is the statistical heart of the paper. If your null
hypotheses were properly calibrated, the test statistics across 12,899 genes should mostly follow
their null distribution, and **λ** (median observed χ² ÷ median expected χ²) should be ≈1. Yours
is **λ=17.3**. That is enormous — test statistics are ~17× over-dispersed. It means the
TCGA+GTEx harmonisation left a massive residual structure (platform/batch/normal-source), so the
BH-adjusted "significant gene" counts are **not trustworthy as literal discoveries**.

**Why the paper is still valid despite λ=17.** This is the subtle, correct move: inflation shifts
and scales the whole distribution of statistics, but it **largely preserves the ranking** of genes
by effect. So instead of trusting the p-value *threshold*, the paper pivots to methods that only
use the **rank**:
- **GSEA / fgsea** (ref 33): asks whether a gene *set* (e.g. Hallmark EMT) sits
  disproportionately at the top or bottom of the ranked list — a rank-based test, robust to a
  global inflation.
- It reports the raw count *with the caveat*, rather than deleting it or pretending λ≈1.

**The lesson for you:** always compute λ after any cross-dataset DE. If it's >~1.2, do not report
adjusted-p gene counts as findings; move to rank/enrichment methods and say why. Reviewers who
know this (yours did) will trust a paper that diagnoses its own inflation far more than one that
hides it.

**The concordance safeguard.** Because GTEx normals could be *driving* the DE (a normal-tissue
artefact rather than tumour biology), you re-ran DE **TCGA-only** and showed the top up/down genes
replicate (top-100/200 = 100%). That proves the ranking reflects tumour biology, not the GTEx
graft. (`tables/DEG_TCGAonly_replication.csv`.)

---

## 4. Functional enrichment — what the tumour is *doing* (BIOLOGY + STATS)

Two flavours, and knowing the difference matters:
- **ORA (over-representation analysis)** — take your significant-gene *list*, ask which GO/KEGG
  terms are over-represented by a hypergeometric test (clusterProfiler, ref 32). Weakness: depends
  entirely on where you drew the significance line — and with λ=17 that line is shaky.
- **GSEA (gene-set enrichment)** — rank *all* genes by the statistic and test whether a set
  clusters at an extreme (fgsea, ref 33). No threshold needed; robust to inflation. **This is your
  primary method**, correctly.

**What you find biologically:**
- **Up in tumour:** cell-cycle / proliferation / DNA-repair programs — FOXM1, TOP2A, MKI67, CDK1,
  the mitotic machinery. A tumour is a proliferation engine.
- **EMT (epithelial-mesenchymal transition):** enriched, and **specifically in the diffuse
  Lauren subtype** — epithelial cells losing adhesion and becoming migratory/mesenchymal, the
  molecular basis of the diffuse phenotype's aggressiveness.
- **Immune microenvironment:** macrophage-weighted infiltration (next section).

This proliferation signature is exactly what makes the **drug-repurposing** arm (§11) nominate
CDK4/6 and PI3K/mTOR inhibitors — drugs that would *reverse* a proliferation-up signature.

---

## 5. Immune deconvolution — and the validation trick you should copy (STATS)

**The problem.** Bulk RNA-seq is a blender: each sample is a mix of cell types. **Deconvolution**
(MCP-counter ref 30; xCell ref 31) estimates the abundance of immune/stromal populations from the
bulk mixture — giving you "how many T cells / macrophages / fibroblasts" per tumour, in silico.

**The trick.** Deconvolution outputs are notoriously method-dependent. So you did the thing most
papers skip: you **validated the in-silico estimates against a real measurement**. TCGA has
pathologist-scored **leukocyte fraction** for many samples. Your MCP-counter T-cell estimate
correlates with measured leukocyte % at **ρ=0.67, p=3.6×10⁻³⁶** (`immune/validation_vs_measured.csv`).
That single correlation converts "here are some deconvolution numbers" into "our deconvolution is
measuring something real." **Copy this pattern** whenever you deconvolve.

(Note the specificity control: the same estimate does *not* correlate with lymphocyte-*infiltration*
score, ρ≈0 — i.e. it tracks total leukocyte content, the thing it should track, not an unrelated
score. Good discriminant validity.)

---

## 6. The prognostic signature — LASSO-Cox and the leakage trap (STATS, central)

**Goal.** Build a gene-expression score that predicts overall survival.

**Cox proportional-hazards model.** The workhorse of survival analysis. It models the **hazard**
(instantaneous risk of death at time t) as a baseline hazard × exp(β·genes). A positive β means
higher expression → higher risk. The **hazard ratio (HR)=exp(β)**: HR=1.9 means ~90% higher risk
per unit (or per SD) of the score. It handles **censoring** (patients still alive at last
follow-up) correctly, which ordinary regression cannot.

**LASSO.** You have thousands of candidate genes and only hundreds of patients — classic p≫n.
**LASSO** (L1-penalised Cox) shrinks most coefficients to exactly zero, selecting a sparse
subset — your **25-gene signature**. The penalty strength λ (here different from the inflation λ!)
is tuned by cross-validation.

**The leakage trap — and why your C=0.72 vs 0.61 gap is the honest number.** If you select genes
and tune λ on the *whole* dataset and then report the C-index on that same data, you have let the
test set peek at training — **information leakage** — and the performance is optimistically
inflated. Your **apparent C-index is 0.72**. But under **nested cross-validation** — where feature
selection and tuning happen *inside* each training fold and the held-out fold is scored blind —
it drops to **C=0.61** (`nested_cv/performance.csv`). That gap (0.72→0.61) *is* the leakage. The
0.61 is the honest, generalisable number, and you lead with it. This is the single most common
way prognostic-signature papers fool themselves; yours doesn't.

**C-index (Harrell's concordance).** The survival analogue of AUC: the probability that, for a
random pair of patients, the one who died first had the higher risk score. 0.5 = coin flip,
1.0 = perfect. **0.61 = modest but real** discrimination. You also report **Uno's C** (0.57), an
IPCW-weighted version robust to the censoring distribution — reporting both is best practice.

---

## 7. External validation and meta-analysis — the honesty engine (STATS)

Internal CV is not enough; a signature must work in *other people's cohorts on other platforms*.

- **Per-cohort** (`validation_multi/cindex_HR_summary.csv`): **ACRG** C=0.61, HR **1.90**
  (1.37–2.62), p=10⁻⁴ ✓; **GSE15459** C=0.58, HR 1.68, p=0.015 ✓; **GSE84437** C=0.53, **null**.
  Two of three validate; one fails.

- **Why GSE84437 fails — and the honest reading.** It's a different platform (Illumina), and
  plausibly a stage-restricted cohort. You *flag the platform confound* rather than explaining the
  failure away, and offer the testable hypothesis that a stromal signature discriminates best in
  earlier-stage disease. That's the right way to report a non-validation.

- **Hartung–Knapp random-effects meta-analysis** (metafor, ref 35): pools the three cohorts'
  per-SD effects. **Pooled HR 1.19 (0.96–1.47), p=0.073 — not significant.** The **Hartung–Knapp**
  adjustment widens the CI appropriately when pooling few (k=3) studies; a naive random-effects
  model would have over-stated significance. You chose the *conservative* estimator and reported
  the *non-significant* pooled result. I²=19% = low heterogeneity.

The tension you present honestly: individually 2/3 validate, but the conservative pooled estimate
crosses 1. Both are true; you show both.

---

## 8. Added clinical value — why "significant" ≠ "useful" (STATS, the clinical lesson)

A signature can be statistically prognostic yet clinically worthless if it adds nothing to what
the clinician already has (**stage**). You test this properly, **out-of-sample in ACRG**:

- **ΔC-index (combined − clinical) = +0.005.** Adding 25 genes to stage+age improves concordance
  by half a percentage point. Negligible.
- **Likelihood-ratio test p=0.002** — "statistically detectable."
- **Decision-curve analysis (DCA):** no net-benefit gain across clinically relevant thresholds.
- **IDI +0.021, NRI +0.175** (borderline, p≈0.04).

**The lesson:** LRT significance with ΔC≈0 and no DCA benefit = a *real but tiny* effect that is
**not clinically actionable**. Your conclusion — "research instrument, not a clinical tool; AJCC
staging remains standard of care" — is exactly what this pattern licenses. Most papers would have
led with "significantly improves prognostication (p=0.002)"; that would be technically true and
practically misleading. Yours doesn't.

---

## 9. WGCNA — finding the stromal module without supervision (STATS + BIOLOGY)

The signature (§6) is *supervised* — built to predict survival. **WGCNA** (ref 28) is
*unsupervised* — it finds co-expression **modules** (groups of genes that rise and fall together
across samples) with no knowledge of outcome, then asks which module happens to be prognostic.
When a supervised signature and an unsupervised module point at the **same biology**, that's
convergent evidence.

**How WGCNA works (TECH):**
1. Correlate every gene with every gene (biweight midcorrelation — robust to outliers).
2. Raise correlations to a **soft power** β so the network becomes approximately **scale-free**
   (a few hub genes, many peripheral) — biologically realistic. You pick β=3; you show the
   scale-free fit R² is flat (0.865–0.886) across β=3–12 with no clear elbow, and prove the
   finding is **invariant to β** (the module survives at β=3,6,9,12). That robustness check
   pre-empts the obvious "you cherry-picked the power" objection.
3. Cluster genes into modules (colours: red, turquoise, …).
4. Summarise each module by its **eigengene** (first PC) and test *that* against survival.

**Result:** the **red module** (263 genes) is the most prognostic (univariable HR 1.31 per SD,
p=9.3×10⁻⁴), and its hub genes are the CAF/ECM markers — CDH11, COL8A1, COL1A2, FNDC1, SPARC,
LUM, BGN, FAP, POSTN. Independently of the LASSO signature, the unsupervised network says
**stroma drives prognosis**. And **3 of the 25 signature genes (SERPINE1, POSTN, MATN3) live in
this module** — the supervised and unsupervised analyses physically overlap.

**The stage-independence caveat (STATS honesty again).** In a multivariable Cox adjusting for
stage/age/leukocyte-fraction/Lauren, the red module attenuates to HR 1.35 (0.95–1.90), **p=0.090
— no longer significant**, while stage stays dominant. You interpret this correctly: the stromal
programme **tracks stage and stromal admixture** rather than adding stage-independent information.
The **tumour-purity** analysis nails it down: the signature correlates *negatively* with tumour
purity (ρ=−0.20) and *positively* with stromal score (ρ=+0.39) — literally measuring stromal
content.

**Module preservation** (ref 29): a Zsummary statistic asking whether the red module's structure
reproduces in external cohorts. **Z=15.9/16.8/17.1** (all ≫10 = "strong") — the module is a real,
reproducible biological unit, not a TCGA artefact. This is the strongest single piece of external
evidence in the paper.

---

## 10. Single-cell — *which cell* is the module, and the circularity you must understand

Bulk WGCNA says "these genes co-vary and predict survival," but bulk cannot say *which cell type*
they come from. **Single-cell RNA-seq (GSE134520)** resolves 43,992 cells into 8 types (Seurat,
ref 37). You then ask: for each red-module hub gene, which cell type expresses it most?

**Result:** nearly all hub genes are **fibroblast-dominant** — the module is a *fibroblast*
programme, confirming the CAF interpretation at single-cell resolution.

**The circularity — understand this precisely, because a reviewer will.** To *call* a cluster
"fibroblast," you annotated it using canonical markers (DCN, LUM, COL1A1, COL1A2, PDGFRB, FAP,
COL3A1). Six of those are *also* among the hub genes you then test. So for those six, "hub gene is
fibroblast-dominant" is **true by construction** — you defined the cluster with them. That is
circular reasoning for that subset.

**The fix (added this cycle):** re-run the localisation on the **23 hub genes NOT used for
annotation**. Result: **all 23/23 remain fibroblast-dominant** (median 0.96). Because the
non-annotation genes — which had no say in defining the cluster — still land in fibroblasts, the
localisation is **robust to the circularity**. (Amusingly, the one gene that *didn't* localise to
fibroblasts in the full panel, FAP → endothelial, was itself one of the removed annotation
markers.) `scrna/gene_dominant_celltype_noncircular.csv`.

**Two honest limits** you now state: (a) the dataset is premalignant/early-stage, so activated
CAFs (FAP⁺) are scarce — you can show fibroblast *identity* but not CAF *activation state*;
(b) enzymatic dissociation under-recovers ECM-embedded fibroblasts, which could itself depress
fibroblast/FAP capture. Hence the careful language: "**fibroblast/stromal programme**," reserving
"activated CAF" for where activation is actually demonstrated.

---

## 11. Drug repurposing (hypothesis-generating, correctly labelled)

**Logic (connectivity mapping):** if a tumour has an "up" gene program, find drugs whose own
transcriptional signature is the **reverse** (they push those genes *down*). Query the tumour
up/down gene lists against LINCS L1000 / DSigDB (ref 19) via Enrichr; rank drugs whose signature
reverses yours.

**Result:** top reversers are **CDK4/6 inhibitors (palbociclib)** and **PI3K/mTOR/FGFR inhibitors
(NVP-BEZ235, dovitinib, pd173074)** — anti-proliferative agents. This makes sense: your tumour-up
program is proliferation (FOXM1/TOP2A/MKI67/CDK1), so cell-cycle inhibitors reverse it. You label
this **hypothesis-generating**, not a therapeutic claim — correct, because the hits target generic
proliferation, not GC-specific or stroma-specific biology.

---

## 12. The tumour microbiome arm — a masterclass in *not* fooling yourself (STATS)

The fashionable claim: tumours harbour a distinct microbiome ("oralization" — oral taxa like
*Fusobacterium/Streptococcus* enriched in tumour). You **reprocessed all 944 16S libraries from
raw FASTQ** through DADA2 (ref 25) → amplicon sequence variants (ASVs) → SILVA taxonomy (ref 26),
rather than trusting a published table. Then you tested the claim honestly and it **fell apart**:

- **Batch confounding (the killer).** Tumour and normal libraries were sequenced on **separate
  flowcells** (of 216 tumour/normal pairs, only 1 shares a flowcell). So "tumour vs normal" is
  perfectly **confounded with sequencing batch** — you cannot tell biology from machine. Proof: a
  classifier that looks great (**AUC 0.916**) is shown to be **batch-driven** — its top features
  are skin/environmental contaminants, and they predict the *flowcell* at 78% accuracy. The
  "tumour microbiome signature" is largely a batch artefact.
- **β-diversity** (community composition, PERMANOVA via vegan, ref 27): phenotype separation
  R²=0.065 **collapses to 0.011 once flowcell is modelled**. Most of the apparent signal was batch.
- **Comparator-dependence.** The oral-taxa signal is *up* vs adjacent-normal, *null* vs matched
  controls, and *down* vs gastritis (the "*H. pylori* paradox") — i.e. it flips with your choice of
  comparison group. Not a fixed biological signature.
- **What DOES survive:** **reduced α-diversity (richness) in cancer** replicates. But even here you
  report **effect sizes**, not just p-values (with n≈900, p-values are near-guaranteed): richness
  Cliff's δ=−0.27, Shannon −0.12, but **Simpson (evenness) δ=−0.07, not significant**. So the
  durable signal is *loss of rare taxa*, not a uniform diversity collapse — a distinction you flag.

**The lesson:** in low-biomass tissue-microbiome work, **batch and comparator choice dominate**.
Your arm is presented not as a biomarker but as a cautionary demonstration — arguably more valuable
to the field than yet another uncontrolled "signature."

*Vocabulary:* **α-diversity** = diversity *within* a sample (richness = # taxa; Shannon/Simpson
add evenness). **β-diversity** = compositional difference *between* samples. **CLR/Aitchison** =
compositional-aware transforms (microbiome data are relative proportions, so standard distances
mislead; centred-log-ratio corrects for this "closure").

---

## 13. Mendelian randomisation — testing *causation*, and why yours is null (STATS)

**The question ORA/correlation can't answer:** is any microbe *causal* for GC, or just correlated?
**Mendelian randomisation (MR)** uses genetic variants (SNPs) associated with an exposure (e.g.
microbe abundance) as **instruments**. Because genotype is randomised at conception and fixed for
life, a SNP→outcome association that runs *through* the exposure implies causation — a "natural
randomised trial," immune to reverse causation and most confounding.

**Three assumptions (memorise these):**
1. **Relevance** — the SNP genuinely associates with the exposure (measured by the **F-statistic**;
   F>10 = not a weak instrument). Yours: F 19.3–20.3, fine.
2. **Independence** — the SNP is not confounded with the outcome.
3. **Exclusion restriction** — the SNP affects the outcome *only* through the exposure, not by
   other paths (**horizontal pleiotropy**). This is the assumption that can break, so you test it
   many ways.

**The methods, and what each guards against (TECH):**
- **IVW (inverse-variance weighted)** — the main estimate; assumes no pleiotropy.
- **MR-Egger** (ref 24) — allows directional pleiotropy; its *intercept* tests for it.
- **Weighted median / mode** — valid if ≤50% of instruments are bad.
- **Cochran's Q** — heterogeneity across instruments (a pleiotropy hint).
- **MR-PRESSO** (ref 23) — detects outlier SNPs and a **global** pleiotropy test.
- **Steiger filtering** — confirms the causal *direction* (exposure→outcome, not reverse).

**Your result: null everywhere.** No exposure (anti-*H. pylori* seropositivity + 5 gut genera)
shows a significant causal effect (smallest IVW p=0.35). The one thing to watch — and you now
report it honestly — is ***H. pylori***: its Cochran's Q is significant (p=0.035) **and** its
MR-PRESSO global test is borderline (p=0.052). Two convergent pleiotropy hints → read its null
with slightly more caution (not "no effect proven," but "this specific null is the least clean").

**The critical honesty (why null ≠ "no effect exists").** Your outcome GWAS has only ~1,029
European cases → **underpowered**. Your CIs (e.g. *H. pylori* OR 0.71–1.30) *span* the small
effects other MR studies report. So the correct statement is "**underpowered to confirm or
exclude**," not "microbes don't cause GC." You also flag a **validity** (not just power) problem:
the MiBioGen instruments are *faecal* gut-genus abundances — a poor proxy for the *gastric mucosal*
niche your 16S data measure. So the MR and tissue arms are **genetically decoupled**; neither
confirms the other. Stating that prevents a reader from over-reading the convergence.

---

## 14. The through-line — what makes this a good paper

Every arm follows the same discipline:
- **DE:** diagnose λ=17, pivot to rank methods, prove ranking ≠ GTEx artefact.
- **Signature:** report leakage-free C=0.61 (not apparent 0.72), and ΔC≈0 → "not a clinical tool."
- **WGCNA:** unsupervised convergence on stroma, but admit stage-dependence via purity.
- **Single-cell:** localise to fibroblasts, then *remove the circular markers and re-prove it*.
- **Microbiome:** show the AUC=0.916 "signature" is a batch artefact; keep only what replicates,
  with effect sizes.
- **MR:** run every sensitivity test, report the borderline *H. pylori* signal, and label the null
  as *underpowered*, not *negative*.

The unifying contribution is **separating association from causation, and biology from artefact**,
in a literature that routinely conflates them. The CAF-prognosis finding is real and externally
validated; everything around it is scoped to exactly what the data support. That restraint is the
paper's strength, and it is why the internal review panels — despite listing many items — never
questioned the *integrity* of the work, only the completeness of its reporting (now closed).

---

## 15. Glossary (fast reference)

- **HR (hazard ratio)** — multiplicative effect on instantaneous death risk; HR>1 = worse.
- **C-index** — survival concordance (AUC analogue); 0.5 chance, 0.61 = modest.
- **λ / genomic inflation** — over-dispersion of test statistics; ≈1 good, 17 = broken thresholds.
- **LASSO** — L1 penalty that zeroes most coefficients → sparse feature selection.
- **Nested CV** — CV with selection/tuning *inside* each fold → leakage-free performance.
- **WGCNA module / eigengene** — co-expressed gene cluster / its first-PC summary.
- **Deconvolution** — estimating cell-type proportions from bulk expression.
- **CAF** — cancer-associated fibroblast; activated stromal cell, matrix-depositing, pro-invasion.
- **EMT** — epithelial→mesenchymal transition; adhesion loss, migration; diffuse-GC hallmark.
- **α/β-diversity** — within-sample / between-sample microbial diversity.
- **ASV** — amplicon sequence variant (exact-sequence 16S feature, replaces OTUs).
- **CLR/Aitchison** — compositional-aware transform/distance for relative-abundance data.
- **MR / IVW / Egger / PRESSO** — causal inference from genetic instruments and its sensitivity tests.
- **F-statistic (MR)** — instrument strength; F>10 = not weak.
- **Pleiotropy** — a SNP affecting the outcome through paths other than the exposure (breaks MR).
- **ComBat** — empirical-Bayes batch-effect correction across datasets.



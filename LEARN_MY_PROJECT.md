# Understanding My Project — From DNA to the Finished Paper

A ground-up guide. Every number quoted here is from your own result files
(verified in `SUBMISSION_VERIFICATION.md`). Read in order; each part assumes only
the parts before it.

---

# PART 0 — THE ABSOLUTE BASICS

## 0.1 What is DNA?

DNA is a long chemical chain that stores instructions for building a living thing.
It is written in an alphabet of just **four letters** — A, T, G, C (adenine,
thymine, guanine, cytosine). These are called **bases** or **nucleotides**.

DNA is **double-stranded**: two chains twisted into a double helix. The strands
pair by a strict rule — **A always pairs with T, G always pairs with C**. So if one
strand reads `ATGC`, the other reads `TACG`. This complementarity is why DNA can be
copied: unzip the two strands, and each one is a template for rebuilding its partner.

Your DNA is packaged into 46 **chromosomes** (23 pairs — one set from each parent).
The total sequence is your **genome**: about 3 billion base pairs.

## 0.2 What is a gene?

A **gene** is a stretch of DNA that contains the recipe for one product — usually a
protein. Humans have roughly **20,000 protein-coding genes**. They occupy only about
2% of the genome; the rest includes regulatory switches and regions whose function
is still debated.

A gene has a name (e.g. `MKI67`, `COL1A1`) and an ID (e.g. `ENSG00000148773`).
In this project you will meet genes like:
- **COL1A1** — makes collagen type I, the main structural fibre in scar tissue
- **MKI67** — makes Ki-67, a protein present only in actively dividing cells
- **SERPINE1, POSTN, FAP** — made by fibroblasts remodelling tissue

## 0.3 The central dogma: DNA → RNA → protein

DNA is the archive; it stays safe in the nucleus. To use a gene, the cell makes a
working copy:

1. **Transcription** — the gene's DNA is copied into **messenger RNA (mRNA)**.
   RNA is single-stranded and uses U (uracil) instead of T.
2. **Translation** — a ribosome reads the mRNA three letters at a time (a
   **codon**), and each codon specifies one **amino acid**. The chain of amino
   acids folds into a **protein**.

Proteins do the actual work: they form structures, catalyse reactions, send
signals. **DNA is the blueprint; protein is the machine; RNA is the work order.**

## 0.4 Gene expression — the single most important idea in your project

Every cell in your body contains the *same* DNA. A skin cell and a stomach cell have
identical genomes. So why are they different?

Because they **express different genes**. Expression = how much mRNA is being made
from a gene right now. A muscle cell transcribes muscle genes heavily and keeps
liver genes switched off.

**Gene expression is a dial, not a switch** — it is a quantity. When we say "COL1A1
is up-regulated in tumours," we mean tumour cells are making *more* COL1A1 mRNA than
normal cells.

Your entire transcriptomic analysis measures one thing: **how much of each gene's
mRNA is present in each sample**. That table — genes as rows, samples as columns — is
called an **expression matrix**. Yours is **18,419 genes × 448 samples**.

## 0.5 What is cancer?

Cancer is what happens when cells stop obeying the rules that govern growth. Normal
cells divide only when instructed and die when damaged. Cancer cells acquire **DNA
mutations** that break these controls, so they:

- divide continuously (uncontrolled proliferation)
- ignore signals to stop or die (evading apoptosis)
- invade neighbouring tissue and spread to distant organs (**metastasis**)

Two important gene categories:
- **Oncogenes** — normally drive growth; when mutated they become stuck "on"
  (accelerator jammed down). Example: *PIK3CA*.
- **Tumour suppressors** — normally restrain growth; when broken, the brakes fail.
  Examples: *TP53*, *CDH1*.

Cancer is usually **not one mutation** but an accumulation over years.

## 0.6 Gastric cancer specifically

Gastric (stomach) cancer is among the leading causes of cancer death worldwide, with
particularly high incidence in East Asia. Key facts you need:

- **Main risk factor:** chronic infection with the bacterium *Helicobacter pylori*,
  which inflames the stomach lining for decades.
- **The Correa cascade** — the accepted stepwise path to gastric cancer:
  normal mucosa → chronic gastritis → atrophic gastritis → intestinal metaplasia →
  dysplasia → carcinoma. Your microbiome analysis samples points along this cascade.
- **Lauren classification** — the classical histological split:
  - **Intestinal type** — cells form gland-like structures; better prognosis
  - **Diffuse type** — cells scatter individually through the wall; worse prognosis;
    associated with loss of the adhesion protein **CDH1** (E-cadherin)
- **TNM staging** — the clinical gold standard for prognosis: **T** = how deep the
  tumour invades, **N** = lymph-node spread, **M** = distant metastasis.
  Stage is the benchmark any new prognostic marker must beat.

## 0.7 The tumour microenvironment — the key to your finding

A tumour is not a pure lump of cancer cells. It is an **ecosystem** containing:

- **Cancer cells** (epithelial origin in gastric cancer)
- **Fibroblasts** — connective-tissue cells that make collagen. Inside a tumour they
  become **cancer-associated fibroblasts (CAFs)**: activated, collagen-pumping cells
  that build a stiff scaffold, help invasion, and shield the tumour from immune cells
  and drugs.
- **Immune cells** — T cells, macrophages; sometimes attacking the tumour, often
  co-opted to help it
- **Blood vessels**, and the **extracellular matrix** (collagen, fibronectin)

The proportion of cancer cells is **tumour purity**. A biopsy that is 40% cancer
cells and 60% stroma is a *mixture*, and bulk sequencing measures the average of
everything present. **This matters enormously for your paper: your prognostic signal
comes from the fibroblast/stromal compartment, not from the cancer cells themselves.**

---

# PART 1 — HOW WE MEASURE EXPRESSION

## 1.1 Sequencing (RNA-seq)

To measure all genes at once:
1. Extract RNA from the tissue.
2. Convert it to DNA (**cDNA**) — more chemically stable.
3. Break into fragments and read millions of short stretches (**reads**, ~100 bases)
   on a sequencing machine.
4. **Align** each read to the reference genome to find where it came from.
5. **Count** reads per gene. More reads = more mRNA = higher expression.

**Normalisation** is essential. If sample A yielded 30 million reads and sample B
15 million, every gene looks "doubled" in A purely from sequencing depth. Methods
like **TPM**, **FPKM**, or DESeq2/limma size factors correct for this, and counts are
usually **log-transformed** because expression spans several orders of magnitude.

## 1.2 Microarrays

An older technology: a chip with probes for known genes; RNA binds and fluoresces
proportional to abundance. Your **GEO validation cohorts** (GSE62254, GSE15459,
GSE84437) are microarray-based, while **TCGA** is RNA-seq. Different platforms give
differently-scaled numbers, which is exactly why you **z-score within each cohort**
before applying your signature (Part 6.6).

## 1.3 Your data sources

| Source | What it is | Your use |
|---|---|---|
| **TCGA-STAD** | The Cancer Genome Atlas, stomach adenocarcinoma — RNA-seq + clinical + survival | Discovery/training: 412 tumours, 36 normals |
| **GTEx** | Genotype-Tissue Expression — healthy tissue from donors | Extra normal references (tumour-vs-normal needs enough normals) |
| **GEO** (GSE…) | Public repository of submitted datasets | 3 independent survival cohorts (n=922) + single-cell |
| **DDBJ / NCBI SRA** | Raw sequence archives | 16S microbiome FASTQ files |
| **IEU OpenGWAS** | Genome-wide association summary statistics | Mendelian-randomisation inputs |

Using public data is legitimate and standard — the science is in the analysis and the
rigour of validation.

---

# PART 2 — STATISTICS FROM SCRATCH

You cannot understand this project without these ideas. They are simpler than they look.

## 2.1 Population, sample, and why we need statistics

You can never measure every gastric-cancer patient. You measure a **sample** and try
to infer something about the **population**. Statistics is the discipline of stating
how much your sample can be trusted.

## 2.2 The null hypothesis and the p-value

Suppose gene X averages 7.7 in tumours and 5.1 in normals. Is that real, or luck?

- **Null hypothesis (H₀):** there is no true difference; the gap is chance.
- **p-value:** the probability of seeing a difference *this large or larger* **if the
  null hypothesis were true**.

Small p = the data would be surprising under "no difference," so we reject H₀.
Convention: p < 0.05.

**What a p-value is NOT** (examiners love this):
- It is *not* the probability that the null hypothesis is true.
- It is *not* a measure of effect size. A tiny, useless difference can have p < 10⁻¹⁵
  in a big sample.
- p = 0.06 does not mean "no effect"; it means "insufficient evidence at this threshold."

Your own example: **MKI67**, tumour mean 7.668 vs normal 5.071, Wilcoxon
**p = 1.32×10⁻¹⁵**. Interpretation: if there were truly no difference, seeing a gap
this large would be astronomically unlikely. Ki-67 marks dividing cells, so this says
tumours are proliferating far more than normal stomach — exactly as biology predicts.

## 2.3 Parametric vs non-parametric tests

- **t-test** assumes data are roughly normally distributed (bell-shaped).
- **Wilcoxon / Mann–Whitney** makes no such assumption; it compares **ranks**.
  Slightly less powerful when normality holds, but safe when it doesn't.

Gene expression is often skewed, so rank-based tests (and the Spearman correlation
below) are common in this field. You used Wilcoxon and Spearman for exactly this reason.

## 2.4 The multiple-testing problem — critical for your project

Test one gene at p < 0.05 and you accept a 5% false-positive risk. Now test
**18,419 genes**. Even if *nothing* were truly different, you would expect
0.05 × 18,419 ≈ **920 genes to look "significant" by chance alone.**

This is why raw p-values are useless in genomics. Two corrections:

- **Bonferroni** — multiply each p by the number of tests (or divide the threshold).
  Very strict; controls the chance of *any* false positive.
- **Benjamini–Hochberg / FDR** — controls the *proportion* of false positives among
  your hits. An **adjusted p (q-value / padj)** of 0.05 means "about 5% of the genes
  I call significant are expected to be wrong."

FDR is the standard in genomics, and it is what your `adj.P.Val` / `padj` columns are.
**Whenever you report a gene list, report adjusted p-values.**

## 2.5 Effect size and confidence intervals

A p-value says "probably not zero." An **effect size** says "how big." Report both.

- **log₂ fold change (log2FC)** — for expression. log2FC = 1 means **2× higher**;
  log2FC = 2 means 4×; log2FC = −1 means half.
- **Hazard ratio (HR)** — for survival (Part 6).
- **Cliff's delta** — a rank-based effect size you used for diversity: how often a
  value from one group exceeds one from the other, scaled to −1…+1.

A **95% confidence interval (CI)** is the range of values compatible with your data.
Its width tells you your precision. Crucially: **if a HR's 95% CI includes 1.0, the
effect is not statistically significant** — because HR = 1 means "no effect."

That single rule explains your most important honest disclosure: pooled
**HR 1.19, 95% CI 0.9616–1.4665**. The interval crosses 1, so despite a hazard ratio
above 1, the pooled result is **not significant** (p = 0.073).

## 2.6 Correlation

**Correlation** measures whether two variables move together, from −1 to +1.
- **Pearson r** — captures straight-line relationships; sensitive to outliers.
- **Spearman ρ** — correlates the **ranks**; captures any consistently increasing
  relationship and resists outliers.

Your immune validation: predicted T-cell content vs pathologist-measured leukocyte
percentage, **Spearman ρ = 0.6656, p = 3.57×10⁻³⁶ (n = 272)**. A strong positive
association: when the algorithm says "more immune cells," the pathologist independently
saw more immune cells.

**Correlation is not causation.** Two things can move together because one causes the
other, or because a third factor drives both (**confounding**), or by coincidence.
This limitation is precisely why Part 9 (Mendelian randomisation) exists in your paper.

## 2.7 Confounding and batch effects

A **confounder** is a hidden variable that influences both things you are comparing,
creating a fake association.

The version that nearly sank your microbiome arm is a **batch effect**: a technical
artefact from *how and when* samples were processed. If tumour samples were sequenced
on one flowcell and control samples on another, then any difference between the
sequencing runs masquerades as a difference between tumour and control.

Your evidence, from `05_rf_metrics_and_batch_sanity.csv`:
- A classifier separates cancer from control with **AUC = 0.9161** — impressive.
- But the same data predicts **which flowcell** a sample came from with
  **77.6% accuracy**, versus a **54.5%** baseline from always guessing the commonest
  flowcell.

The technical variable is nearly as predictable as the biology. **Therefore the
"biological" signal is substantially a batch artefact.** Recognising this is the most
scientifically mature moment in your paper.

---

# PART 3 — DIFFERENTIAL EXPRESSION

## 3.1 The question

Which genes differ between tumour and normal tissue?

## 3.2 The method

For each gene, compare expression across groups, get a log2FC and a p-value, then
FDR-adjust across all genes. Tools: **limma** (linear models, borrows information
across genes — well suited to microarrays and log-transformed data), **DESeq2** /
**edgeR** (negative-binomial models for raw counts).

Because you have thousands of genes but few samples, these tools use **empirical Bayes
shrinkage**: they stabilise each gene's variance estimate using information from all
other genes. This is why they outperform running 18,419 separate t-tests.

## 3.3 Your result

Integrating TCGA + GTEx: **3,722 genes up, 4,025 down** of 12,899 tested.
TCGA-only contrast: **2,134 up / 2,362 down** of 21,446.

## 3.4 The volcano plot

The standard visualisation (your Figure 6A): x-axis = log2FC (biological size),
y-axis = −log₁₀(adjusted p) (statistical confidence). Points at top-left and
top-right — big change *and* confident — are the interesting genes. A gene can be
highly significant but biologically trivial (bottom-right), which is why the plot
shows both axes at once.

## 3.5 The inflation problem you correctly disclosed

Comparing TCGA tumours against GTEx normals mixes two studies collected differently.
That inflates test statistics. You measured the inflation as **λ = 17.3** (λ ≈ 1 means
no inflation) and ran a **permutation null** — randomly shuffling labels 100 times to
see what λ arises by chance: mean 1.06, 95th percentile 1.58.

Since 17.3 is far above the permutation ceiling, the signal exceeds what shuffling
produces — but the inflation is real and you **say so**, and you confirm the gene
ranking against TCGA-only and two other GEO cohorts (concordance r = 0.62–0.81).
That is how an honest paper handles a known weakness: quantify it, disclose it,
triangulate around it.

---

# PART 4 — PATHWAYS AND GENE SETS

## 4.1 Why not just read the gene list?

A list of 3,722 genes is not a finding. Genes work in **pathways** — coordinated
groups. Interpretation means asking which *processes* changed.

## 4.2 Two approaches

**Over-representation analysis (ORA)** — take your significant genes; ask whether any
pathway appears more often than chance (Fisher's exact test). Weakness: depends
entirely on your arbitrary significance cutoff.

**Gene Set Enrichment Analysis (GSEA)** — better, and what your headline uses. Rank
*all* genes by change, then walk down the ranked list; a running score rises when you
hit a pathway member and falls when you don't. If a pathway's genes cluster near the
top, the score peaks strongly.
- **NES (normalised enrichment score)** — positive = pathway shifted up in tumours,
  negative = shifted down. Magnitude ≈ strength.
- GSEA needs no cutoff and detects coordinated modest shifts that ORA misses.

## 4.3 Your results (integrated TCGA+GTEx, `fgsea_Hallmark_integrated.csv`)

| Pathway | NES | Meaning |
|---|---|---|
| E2F targets | **+3.75** | Cell-cycle entry genes strongly up |
| G2M checkpoint | **+3.64** | Mitotic machinery up |
| MYC targets | **+2.67** | Growth-driving programme up |
| Oxidative phosphorylation | **−2.54** | Mitochondrial energy production down |
| Fatty-acid metabolism | **−2.58** | Normal metabolic function down |

**The biology:** E2F and G2M are the engines of cell division; MYC drives growth.
Their coordinated up-regulation is the transcriptional signature of proliferation —
the textbook cancer phenotype. Simultaneously, oxidative phosphorylation falls, which
matches the **Warburg effect**: tumours shift from efficient mitochondrial respiration
toward glycolysis, trading energy efficiency for fast biosynthesis. Down-regulated
fatty-acid metabolism reflects loss of the specialised functions of normal stomach lining.

**Diffuse vs intestinal:** EMT (epithelial–mesenchymal transition) is the top set in
diffuse tumours, **NES +3.17, padj 1.6×10⁻⁴⁰**. EMT is the programme by which
epithelial cells lose adhesion and become migratory — precisely the histology that
defines diffuse gastric cancer (scattered, non-cohesive cells, CDH1 loss). **Your
data independently rediscovers the classical Lauren classification from expression
alone.** That is a strong internal validity check.

---

# PART 5 — SURVIVAL ANALYSIS

## 5.1 Why survival needs its own statistics: censoring

You want to know whether high-risk patients die sooner. The complication: when the
study ends, many patients are still alive. You know patient A survived *at least* 5
years, not how long they will live. That is **right-censoring** — partial information.

Throwing censored patients away wastes data and biases results. Survival analysis
uses them properly: each patient contributes information for as long as they were
observed.

## 5.2 Kaplan–Meier curves

The **KM curve** plots the probability of surviving over time. It steps down at each
death; censored patients are tick marks. Split patients into high- and low-risk and
plot both curves: if they separate, risk groups have different survival.

The **log-rank test** asks whether two KM curves differ more than chance.
Your ACRG cohort: log-rank **p = 8.63×10⁻⁵** — clear separation.

## 5.3 Cox proportional-hazards regression

KM handles groups; **Cox regression** handles continuous predictors and adjusts for
multiple variables at once.

The **hazard** is the instantaneous risk of the event at a given moment. Cox models
the hazard as a baseline multiplied by predictor effects, giving a **hazard ratio (HR)**:

- **HR = 1** — no effect
- **HR = 2** — double the hazard per unit increase
- **HR = 0.5** — halved hazard (protective)

For continuous scores we report **HR per standard deviation (HR/SD)** so the number
means "risk change for a typical 1-SD increase in score."

**Two different HRs appear in your paper — do not confuse them:**

| Measure | ACRG value | What it means |
|---|---|---|
| **Median-split, high-vs-low** (unadjusted) | **1.8965 (1.3704–2.6245)**, p=1.13×10⁻⁴ | Patients above the median risk score have ~1.9× the hazard of those below |
| **Per SD, age/stage-adjusted** | **1.3012 (1.1004–1.5386)**, p=0.0021 | Each 1-SD rise in score raises hazard ~30%, after adjusting for age and stage |

The median-split number is larger because dichotomising contrasts the extremes of the
distribution; the per-SD adjusted number is the more conservative and more meaningful
estimate, and it is the one carried into the meta-analysis. In both cases the CI excludes
1, so the ACRG result is significant. (A third figure, **1.76**, is the age/stage-adjusted
*median-split* HR — same split, with adjustment.)

**The proportional-hazards assumption:** Cox assumes the HR is constant over time.
You tested it (`cox.zph`) and found it **violated** — p = 0.0018 for the signature alone,
p = 0.003 for the age/stage-adjusted model (the manuscript quotes the adjusted figure) — so you fitted a
time-varying model, which showed the signature is strongly prognostic early
(**HR 1.49 at 12 months**) and fades to nothing later (**1.03 at 36 months, 0.87 at
60 months**). Rather than hide an assumption violation, you characterised it and
reported that your marker is an *early*-risk marker. This is exactly the kind of
rigour reviewers respect.

## 5.4 Multivariable adjustment

A marker is only useful if it adds information **beyond what clinicians already have**.
So you adjust for **age** and **stage**: fit a Cox model containing all three and ask
whether the signature retains an independent effect. Your paper reports **both** the
unadjusted median-split HRs and the age/stage-adjusted per-SD HRs (Table 1), and it is the
*adjusted* ones that are pooled — a much stronger claim than an unadjusted one. When you
quote an HR, always say which of the two you mean.

## 5.5 The C-index — how to judge a prognostic model

**Harrell's concordance index (C-index)** asks: given two patients, how often does the
model correctly rank who dies first?

- **C = 0.5** — coin flip, no predictive ability
- **C = 0.7** — clinically useful
- **C = 1.0** — perfect

**Uno's C** is a variant that corrects for censoring distribution — usually slightly
lower and more conservative.

Your honest numbers: **Harrell C = 0.6112 (95% CI 0.5620–0.6589)**, Uno
**C = 0.5727**. Modest — better than chance, well short of clinically decisive. You
report this rather than the flattering number, and the next part explains why that
distinction is the heart of your paper.

---

# PART 6 — MACHINE LEARNING, OVERFITTING, AND THE MOST IMPORTANT LESSON IN YOUR PAPER

## 6.1 The goal

Build a score from gene expression that predicts survival.

## 6.2 Why ordinary regression fails here: p ≫ n

You have ~20,000 candidate genes and a few hundred patients. With more predictors
than patients, a model can fit the training data *perfectly* while learning nothing
generalisable. It memorises noise. This is **overfitting**.

## 6.3 Regularisation and LASSO

**Regularisation** penalises complexity. **LASSO** (L1 penalty) adds a penalty
proportional to the sum of absolute coefficients. Its special property: it drives many
coefficients **exactly to zero**, performing automatic **feature selection**.

**LASSO-Cox** = LASSO applied to a Cox survival model. The penalty strength **λ** is
chosen by cross-validation. Your model kept **25 genes** (16 positive coefficients =
higher expression → worse survival; 9 negative = protective).

## 6.4 Cross-validation

Split data into k folds (say 10). Train on 9, test on the held-out one, rotate, and
average. This estimates performance on unseen data using only your training set.

## 6.5 Data leakage and why nested CV is essential — THE KEY IDEA

Here is the trap that invalidates a great many published signatures.

Suppose you use all your data to pick which genes to include, *then* cross-validate
the model built from those genes. The gene-selection step already saw the test folds.
Information has **leaked**. The resulting C-index is optimistically biased —
sometimes wildly.

**Nested cross-validation** fixes it. There are two loops:
- **Inner loop:** on the training portion only, do everything data-dependent —
  select genes, tune λ.
- **Outer loop:** evaluate the finished pipeline on data the inner loop never touched.

Every choice is made without seeing the evaluation data, so the estimate is honest.

**Your numbers tell the whole story:**
- Apparent C (the leaky way): **0.72**
- Nested-CV C (the honest way): **0.6112**

**The optimism is ~0.11 of C-index — a large inflation.** Many published papers
report the 0.72-style figure. You report **0.61** and say why. When you present this,
say it plainly: *the difference between those two numbers is the difference between an
honest paper and an over-claimed one.*

Independent confirmation: you re-ran a **previously published 5-gene gastric signature**
on your identical data and pipeline. Its optimism-corrected C was **0.4807 — below
chance** (apparent 0.5446). This is powerful evidence that leakage-uncontrolled
signatures routinely fail to generalise, and it makes your modest 0.61 look genuinely
better rather than merely smaller.

## 6.6 External validation — the real test

Cross-validation still uses one dataset. **External validation** applies the frozen
model to completely independent cohorts. Because platforms differ, you **z-score
within each cohort** (subtract that cohort's mean, divide by its SD) so the score is
comparable — and crucially, no information crosses between cohorts.

Your results (`cindex_HR_summary.csv`):

| Cohort | n | Events | C-index | HR high-vs-low, median split (95% CI) | HR per SD, age/stage-adj (95% CI) | Verdict |
|---|---|---|---|---|---|---|
| ACRG/GSE62254 | 300 | 152 | 0.6079 | 1.8965 (1.3704–2.6245) | 1.3012 (1.1004–1.5386) | **validated** |
| GSE15459 | 191 | 95 | 0.5752 | 1.6757 (1.1067–2.5373) | 1.1996 (0.9736–1.4782) | **validated** (median split) |
| GSE84437 | 431 | 207 | 0.5297 | 1.1087 (0.8440–1.4565) | 1.1096 (0.9726–1.2658) | **failed** |

**Two of three validated.** You report the failure, and you investigated it: GSE84437
is 89% pT3–T4 (advanced disease), so the cohort has little of the early-stage
variation your signature detects. You even tested stage-stratified subgroups and it
still failed — reported as a clean negative rather than buried.

## 6.7 Meta-analysis — combining cohorts honestly

**Meta-analysis** pools effects across studies. A **random-effects** model assumes
true effects vary between studies; the **Hartung–Knapp** adjustment widens the
interval to account for having only a few studies — the conservative, correct choice
for k = 3.

Your pooled result: **HR 1.1875, 95% CI 0.9616–1.4665, p = 0.073, I² = 19.2%**
(I² measures between-study heterogeneity; 19% is low).

**The CI includes 1, so the pooled effect is not statistically significant.** You
state this. A less careful author would have highlighted only the significant ACRG
result. This is the single most defensible sentence in your manuscript.

## 6.8 Does it help clinically? Decision-curve analysis

Statistical significance ≠ clinical usefulness. The right question: does adding the
signature to stage improve decisions?

- **ΔC (change in C-index)** when adding signature to a clinical model: **+0.0045**
- Likelihood-ratio test: **p = 0.0020** (statistically detectable)
- **IDI@3y = 0.021**, **continuous NRI@3y = 0.175**
- **Decision-curve analysis (DCA)**: plots *net benefit* across decision thresholds.

Interpretation: the signature adds *statistically detectable* but *clinically
negligible* information over stage. **ΔC of 0.005 is not a clinical advance**, and
your paper says so. Reporting a significant p-value alongside a trivial effect size —
and interpreting it correctly — is exactly the maturity reviewers look for.

---

# PART 7 — NETWORKS: WGCNA AND YOUR PRIMARY FINDING

## 7.1 The idea

Instead of one gene at a time, ask which genes **move together** across samples.
Genes that are co-expressed are often co-regulated and functionally related.
**WGCNA** (weighted gene co-expression network analysis) finds such groups, called
**modules**.

## 7.2 How it works

1. Correlate every gene with every other gene. You used **bicor** (biweight
   midcorrelation), an outlier-robust alternative to Pearson.
2. Raise correlations to a power **β** (soft-thresholding) to emphasise strong links.
   You used **β = 3**, with scale-free fit R² ≈ 0.88.
3. Convert to **topological overlap** — two genes are close if they also share
   neighbours (more robust than pairwise correlation alone).
4. **Cluster** into modules, conventionally labelled by colour (your "red" module).
5. Summarise each module by its **eigengene** — the first principal component, i.e. a
   single number per sample capturing the module's overall activity.
6. Correlate eigengenes with traits (stage, survival, subtype).

## 7.3 Module preservation — why this is your strongest result

A module found in one dataset could be an artefact of that dataset. **Module
preservation** tests whether the same co-expression structure exists in independent
cohorts, summarised by **Zsummary**: Z > 2 = weak, Z > 10 = **strong** preservation.

Your red module (`preservation_summary_RED.csv`):

| Cohort | Zsummary | Verdict |
|---|---|---|
| ACRG/GSE62254 | **15.853** | STRONG |
| GSE15459 | **16.797** | STRONG |
| GSE84437 | **17.083** | STRONG |

All three far exceed 10. And the module's eigengene independently predicts survival
in all three (`module_eigengene_cox_external.csv`): **HR/SD 1.274** (p = 0.0022),
**1.548** (p = 6.8×10⁻⁵), **1.237** (p = 0.0021).

**This is your primary finding, and note something important:** the *module* replicates
in all three cohorts — including GSE84437, where the 25-gene signature failed. The
underlying **biology is more robust than the specific gene-score built from it**. That
is a genuinely interesting scientific point, and it is why the paper is framed around
the stromal programme rather than around the signature.

## 7.4 What the module is

Its genes are collagens and matrix remodellers — **COL1A1, COL1A2, COL3A1, POSTN,
FAP, LUM, DCN, THBS2, FN1, VCAN, BGN, SPARC, CDH11**. This is the transcriptional
fingerprint of **activated fibroblasts building extracellular matrix**: CAF biology.

---

# PART 8 — DECONVOLUTION AND SINGLE-CELL: PROVING WHERE THE SIGNAL COMES FROM

## 8.1 The bulk-tissue problem

Bulk RNA-seq averages every cell in the biopsy. A "high collagen" tumour might have
more fibroblasts, or the same number working harder. Bulk data cannot distinguish these.

## 8.2 Deconvolution

**Deconvolution** algorithms (CIBERSORT, MCP-counter, xCell, ESTIMATE) estimate cell-type
proportions from bulk expression using reference profiles of known cell types.

These are **estimates**, so they need validating. You did the right thing: compared
predictions against **pathologist-measured leukocyte percentage** from the same
tumours — **Spearman ρ = 0.6656, p = 3.57×10⁻³⁶ (n = 272)**. An independent,
non-computational ground truth. You also honestly report that CD8 predictions
correlate with leukocyte % (ρ = 0.468) but **not** with lymphocyte-infiltration %
(ρ = 0.0198, p = 0.74) — disclosing where the method is weak.

## 8.3 Single-cell RNA-seq

**scRNA-seq** measures expression in **individual cells**, so you can see which cell
type expresses a gene rather than inferring it. Pipeline: isolate cells → barcode each
→ sequence → cluster cells by expression → label clusters using **marker genes**
(e.g. EPCAM = epithelial, PTPRC/CD45 = immune, DCN/LUM/COL1A1 = fibroblast).

Your dataset (GSE134520): **43,992 cells**, 8 cell types — **68.6% epithelial,
4.2% fibroblast**.

## 8.4 Avoiding circular reasoning — a subtle and impressive step

Here is a trap you avoided. If you use COL1A1 to *label* a cluster "fibroblast," then
announce "COL1A1 is fibroblast-specific," you have proven nothing — the conclusion was
built into the labelling. That is **circularity**.

So you removed every gene used for annotation and tested only the remainder. Result
(`gene_dominant_celltype_noncircular.csv`): of the **23** hub genes not used for
annotation, **23/23 remain fibroblast-dominant**, median fraction **0.960**.

**Independent, non-circular proof that your prognostic module is fibroblast-derived.**
Three methods — bulk co-expression, deconvolution, single-cell — converge on the same
answer. That convergence is what makes the finding credible.

---

# PART 9 — THE MICROBIOME ARM AND HOW YOU CAUGHT AN ARTEFACT

## 9.1 16S rRNA sequencing

Bacteria all carry the **16S ribosomal RNA** gene. It has regions conserved across all
bacteria (so universal primers bind) interspersed with **variable regions** (V3–V4,
V5–V6) that differ between species. Amplify a variable region, sequence it, and the
sequences tell you which bacteria are present.

Cheap and robust, but limited: usually only **genus**-level resolution, and it counts
gene copies, not functions.

## 9.2 Your pipeline (DADA2)

Trim primers → truncate at quality drop-off (260/220) → learn the run's error model →
infer exact sequence variants (**ASVs**) → remove chimaeras → assign taxonomy against
**SILVA v138.1** → remove mitochondrial/chloroplast/non-bacterial reads.

From 944 libraries: **13,487,331 starting reads**, 11,795 ASVs → **897 samples**,
**314 genera** retained. You validated your reprocessing against the original study's
published read-tracking table: **Pearson r = 0.983**. That is a reproducibility check
most papers skip.

## 9.3 Compositional data — a subtlety many papers get wrong

Sequencing returns **relative** abundances summing to 100%. If one taxon blooms,
everything else's percentage falls *even if unchanged in absolute terms*. Standard
statistics on such data produce spurious correlations.

The fix is a **log-ratio transform**, e.g. **CLR (centred log-ratio)**: divide each
value by the sample's geometric mean and take the log. You used CLR, plus **Aitchison
distance** (Euclidean distance on CLR values) alongside Bray–Curtis. Using
compositionally-aware methods is a mark of a careful microbiome analysis.

## 9.4 Diversity

- **Alpha diversity** — diversity *within* a sample. **Observed richness** = number of
  taxa. **Shannon** weights by evenness. **Simpson** weights dominant taxa more.
- **Beta diversity** — dissimilarity *between* samples (Bray–Curtis, Aitchison),
  tested with **PERMANOVA**, which reports **R²** (variance explained).

Your findings (`02_alpha_effectsizes_cascade.csv`):

| Metric | Median (non-ulcer → cancer) | Cliff's δ | p |
|---|---|---|---|
| Observed richness | 47 → 30 | −0.27 | 3.3×10⁻⁷ |
| Shannon | 2.176 → 1.923 | −0.122 | 0.0211 |
| Simpson | 0.781 → 0.756 | −0.072 | **0.172 (ns)** |

Richness clearly falls; Shannon modestly; **Simpson is not significant**. You report
this discordance rather than citing only the significant metrics — the honest reading
is that **rare taxa are lost while dominant-taxon structure is largely preserved**.

## 9.5 The batch-effect discovery

As covered in 2.7: classifier AUC **0.9161** for cancer-vs-control, but flowcell
predictable at **77.6%** vs **54.5%** baseline. Adjusting for flowcell collapsed the
beta-diversity effect (**Bray R² 0.0647 → 0.0109**). Of 216 tumour/adjacent-normal
pairs, only 1 shared a flowcell — so batch and biology are almost perfectly
**confounded** and cannot be statistically separated.

Your random-forest figure (S20) shows the smoking gun: the top discriminating genera
are **Dietzia, Serinicoccus, Methylobacterium, Microbacterium, Sphingomonas, Serratia,
Cutibacterium** — classic **environmental and skin contaminants**, not gastric flora.
A real biological signal would be led by gut organisms; a reagent/batch artefact looks
exactly like this.

**Replication attempt:** in an independent Italian cohort — null (Bray R² 0.018,
p = 0.80). In a Portuguese cohort — reduced diversity replicated (Shannon p = 0.0044;
Bray R² 0.145, p = 0.001). So the diversity finding has some external support, while
the composition finding does not.

**This is the correct scientific outcome**: you did not delete the analysis, and you
did not over-claim it. You characterised its limits.

---

# PART 10 — CAUSALITY: MENDELIAN RANDOMISATION

## 10.1 The problem

Bacteria are more abundant in cancer tissue. Did the bacteria **cause** the cancer, or
did the cancer create an environment where they thrive (**reverse causation**), or does
something else drive both (**confounding**)? Observational data cannot distinguish these.

## 10.2 The idea

A randomised trial would settle it, but you cannot randomise people to a gut microbe.
**Mendelian randomisation (MR)** uses a natural experiment: **genetic variants are
randomly assigned at conception** and fixed for life.

If variants that raise your abundance of *Streptococcus* also raise gastric-cancer
risk, that supports causation — because the genotype cannot have been caused by the
cancer, and genotype is largely independent of lifestyle confounders. Genetic variants
used this way are **instrumental variables**.

## 10.3 The three assumptions (memorise these)

1. **Relevance** — the variant genuinely affects the exposure. Checked with the
   **F-statistic**; F > 10 is the rule of thumb. Yours: F = 19.3–22.8 ✓
2. **Independence** — the variant is not associated with confounders.
3. **Exclusion restriction** — the variant affects the outcome *only* through the
   exposure. Violated by **horizontal pleiotropy** (a variant affecting the outcome by
   another route). This is the hardest to guarantee.

## 10.4 Methods and sensitivity analyses

- **IVW (inverse-variance weighted)** — the main estimate; a weighted average.
- **MR-Egger** — its intercept tests for directional pleiotropy (intercept ≠ 0 = trouble).
- **Weighted median** — valid if up to 50% of instruments are invalid.
- **Weighted mode** — relies on the commonest effect estimate.
- **MR-PRESSO** — detects outlier variants and tests global pleiotropy.
- **Leave-one-out** — drop each SNP in turn; a result that hinges on one SNP is fragile.
- **Cochran's Q** — heterogeneity among instrument estimates.

Agreement across methods with different assumptions is what makes an MR result credible.

## 10.5 Your results — a clean null

Six exposures (anti-*H. pylori* IgG, *Streptococcus*, *Fusobacterium*, *Prevotella*,
*Veillonella*, *Lactobacillus*) against gastric cancer (1,029 cases / 475,087 controls):

**All six IVW estimates non-significant; the smallest p = 0.348.** MR-PRESSO global
tests non-significant for five of six. Replicated in an East-Asian outcome cohort.

**Interpretation — be precise here.** This is **not** "microbes definitely don't cause
gastric cancer." With 1,029 cases and modest instrument strength, power is limited.
The correct statement is: **at current instrument and outcome power, these data provide
no evidence for a causal effect.** Absence of evidence ≠ evidence of absence. Your
manuscript words it this way, which is why it is defensible.

Note the honest tension you disclose: *H. pylori* is an established gastric-cancer
cause epidemiologically, yet your MR is null for anti-*H. pylori* IgG. The resolution
is that the instrument measures *antibody response*, not lifetime infection burden —
a real limitation, stated rather than hidden.

---

# PART 11 — DRUG REPURPOSING AND DEPENDENCY

## 11.1 Signature reversal

If a tumour's expression programme is disease, a drug that produces the *opposite*
expression change might treat it. Databases (**LINCS L1000**, CMap) record expression
changes after drug treatment in cell lines. Query your tumour signature and rank drugs
whose profile reverses it.

Your top classes: **PI3K/mTOR inhibitors** (NVP-BEZ235, PI-103, GDC-0941, torin-2),
**CDK4/6 inhibitors** (palbociclib), **FGFR/multikinase** (PD-173074, dovitinib), and
resveratrol.

**Honest caveat you state:** because your up-signature is proliferation-dominated,
these hits are broadly **anti-proliferative rather than gastric-specific** — and are
hypothesis-generating only, with no experimental validation.

## 11.2 DepMap — are the targets real vulnerabilities?

Predicting a drug is cheap; showing the target matters is better. **DepMap** performs
genome-wide **CRISPR knockout** screens across hundreds of cancer cell lines and
reports a **gene-effect (Chronos/CERES) score**: how much killing the gene harms the
cell. More negative = more essential. Below about **−0.5** = dependent; **−1.0** ≈ a
median common-essential gene.

Your gastric lines (n = 35), from `gastric_dependency.csv`:

| Gene | Gene effect | Fraction dependent | Verdict |
|---|---|---|---|
| **PIK3CA** | −0.742 | 57% | dependent **and selective** for gastric |
| **MTOR** | −1.184 | 100% | dependent (common-essential) |
| **CDK4** | −0.825 | 66% | dependent |
| **CDK6** | −0.548 | 51% | dependent |
| FGFR1–4 | −0.04 to −0.22 | ≤11% | **not** dependent |

**PI3K/mTOR and CDK4/6 are genuine gastric dependencies; FGFR is not.** This upgrades
the drug section from pure in-silico speculation to experimentally-anchored — and the
FGFR negative shows you reported the result rather than the hoped-for answer.

---

# PART 12 — WHAT YOUR PROJECT ACTUALLY SHOWS

## 12.1 The findings, honestly stated

1. **Primary (strong).** An externally-validated **stromal/CAF programme** underlies
   gastric-cancer prognosis. Preserved in three independent cohorts (Zsummary
   15.9–17.1), prognostic in all three (HR/SD 1.24–1.55), and localised to fibroblasts
   by non-circular single-cell analysis (23/23, median 0.96). **This replicates where
   the signature does not.**
2. **Secondary (modest, honestly bounded).** A **25-gene signature** with nested-CV
   **C = 0.61** (not the apparent 0.72), validating in 2 of 3 cohorts, pooled
   **HR 1.19 (0.96–1.47), non-significant**, prognostic early and fading with time, and
   adding **negligible clinical value over stage (ΔC +0.005)**.
3. **Cautionary negative 1.** Apparent tumour-microbiome dysbiosis is substantially a
   **sequencing-batch artefact** (contaminant-led classifier; 78% flowcell-predictable).
4. **Cautionary negative 2.** **MR finds no causal microbial effect** at current power.
5. **Supporting biology.** Proliferation programmes up / oxidative metabolism down;
   EMT specific to the diffuse subtype; macrophage-weighted immune microenvironment
   with deconvolution validated against pathology (ρ = 0.67).
6. **Benchmark.** A prior published signature scores **C = 0.48** (chance) on the same
   data — context that makes 0.61 meaningful.

## 12.2 The intellectual contribution

Your paper's contribution is **not** "we found a new biomarker." It is:

> **A rigorous, multi-modal triangulation of gastric-cancer prognostic biology that
> separates what replicates from what does not, and association from causation.**

You have one robust replicated finding, one honestly-bounded modest predictor, and two
well-characterised negatives. In a literature crowded with over-claimed signatures and
uncontrolled microbiome "dysbiosis," being the paper that measured its own optimism
(0.72 → 0.61), caught its own batch artefact, and tested causality instead of assuming
it is a genuine contribution.

## 12.3 Honest weaknesses (know these before anyone asks)

- Modest discrimination; not clinically actionable.
- One of three validation cohorts failed.
- Non-significant pooled meta-analysis.
- The CAF–prognosis link is **confirmatory** — established biology, rigorously
  re-demonstrated, not newly discovered.
- Signature genes are unstable across resamples (13/25 selected >50%, none >80%) —
  the *programme* is robust, the *specific gene list* is not.
- Microbiome arm is confounded and only partly replicated.
- MR is underpowered (1,029 cases).
- All retrospective, public, secondary data; no wet-lab validation.

**Weaknesses that are measured and disclosed are strengths of the paper.** Undisclosed
ones are what get papers rejected or retracted.

---

# PART 13 — VIVA / DEFENCE QUESTIONS

**Q: Why is your C-index only 0.61?**
Because that is the honest, leakage-free nested-CV estimate. The apparent value was
0.72; the ~0.11 gap is optimism from letting gene selection see the evaluation data.
Most published signatures report the 0.72-style number. For calibration, a prior
published signature scores 0.48 on identical data.

**Q: What is the difference between apparent and nested-CV performance?**
Apparent performance evaluates on data used to build the model, so selection and tuning
have already seen it — biased upward. Nested CV puts all data-dependent choices in an
inner loop and evaluates on outer-fold data never used for any decision.

**Q: Your meta-analysis isn't significant. Doesn't that sink the paper?**
It bounds the signature claim, which is why the paper is framed around the *module*,
not the signature. The module preserves strongly in all three cohorts (Z 15.9–17.1)
and is prognostic in all three — including the cohort where the signature failed.

**Q: Why did GSE84437 fail?**
It is 89% pT3–T4, so it lacks the early-stage variation the signature detects. We
tested stage-stratified subgroups; it still failed, and we report that as a negative.

**Q: How do you know the microbiome result is an artefact and not biology?**
Three lines. (1) The data predicts sequencing flowcell at 78% vs 55% baseline.
(2) Flowcell adjustment collapses the effect (Bray R² 0.065 → 0.011). (3) The top
discriminating genera are environmental/skin contaminants (Dietzia, Serinicoccus,
Sphingomonas), not gastric flora. Also, of 216 tumour/normal pairs only 1 shares a
flowcell, so batch and biology are nearly inseparable by design.

**Q: Does your MR prove microbes don't cause gastric cancer?**
No. It shows no evidence of causation at current power — 1,029 cases with moderate
instruments. Absence of evidence is not evidence of absence.

**Q: Isn't the CAF finding already known?**
Yes, and we say so. Our contribution is the rigour: external module preservation,
non-circular single-cell localisation, leakage-controlled prognostic modelling, and
explicit separation of what replicates from what does not.

**Q: Why is stroma prognostic rather than the cancer cells?**
CAFs build a collagen-rich matrix that promotes invasion, stiffens tissue, impedes
drug delivery, and excludes immune cells. A stroma-rich tumour is a more permissive
tumour. Bulk expression captures this because the biopsy contains both compartments —
which is also why we adjusted for tumour purity (signature survived, HR 2.97;
purity itself non-significant, p = 0.35).

**Q: What would you do next?**
Protein-level validation (immunohistochemistry for POSTN/FAP on tissue), CAF
subclustering in single-cell data to find which fibroblast subtype carries the signal,
better-powered MR as GWAS grow, and prospective testing of whether the module adds to
stage in an early-stage-enriched cohort.

---

# GLOSSARY

**Adjusted p (padj/q)** — p corrected for multiple testing (FDR).
**ASV** — amplicon sequence variant; exact 16S sequence.
**AUC** — area under ROC curve; 0.5 = chance, 1 = perfect classification.
**Batch effect** — technical artefact from processing conditions.
**CAF** — cancer-associated fibroblast.
**Censoring** — outcome unobserved by end of follow-up.
**C-index** — probability the model ranks two patients' survival correctly.
**CLR** — centred log-ratio; compositional transform.
**Confounder** — hidden variable driving both compared quantities.
**Cox regression** — survival model giving hazard ratios.
**DCA** — decision-curve analysis; net clinical benefit.
**Deconvolution** — estimating cell-type proportions from bulk expression.
**DEG** — differentially expressed gene.
**Eigengene** — first principal component summarising a module.
**EMT** — epithelial–mesenchymal transition; adhesion loss, motility gain.
**FDR** — false discovery rate; expected proportion of false positives.
**GSEA** — gene set enrichment analysis; NES = normalised enrichment score.
**Hazard ratio (HR)** — multiplicative effect on instantaneous risk; 1 = no effect.
**I²** — % of meta-analysis variability from between-study heterogeneity.
**Instrumental variable** — genetic proxy for an exposure in MR.
**LASSO** — L1-penalised regression; shrinks coefficients to zero (feature selection).
**Leakage** — evaluation data influencing model construction; inflates performance.
**Log2FC** — log₂ fold change; 1 = twofold up.
**MR** — Mendelian randomisation; genetic instruments for causal inference.
**Nested CV** — cross-validation with tuning/selection in an inner loop.
**Overfitting** — fitting noise; good training, poor generalisation.
**PERMANOVA** — permutation test for group differences in community composition.
**Pleiotropy** — a variant affecting the outcome outside the exposure pathway.
**Purity** — proportion of cancer cells in a sample.
**Spearman ρ** — rank correlation.
**TNM** — tumour/node/metastasis clinical staging.
**WGCNA** — weighted gene co-expression network analysis; finds modules.
**Zsummary** — module-preservation statistic; >10 = strong.

# Anticipated Reviewer Comments — Pre-Submission Response Map

*Internal working document. For each comment a reviewer is likely to raise, this maps
the answer already present in the manuscript (section / table / figure) so a formal
point-by-point response can be assembled quickly. Values are copied from the result
files, not from memory.*

---

## A. Near-certain comments (expect from almost every reviewer)

**A1. "No experimental / wet-lab validation."**
- Status: **Disclosed limitation, not a fixable gap.**
- Where addressed: Limitations §4.1 item (6) — "No experimental validation was performed; all inference is computational."
- Response stance: The work is a rigorous *in-silico* reanalysis whose central claim is
  externally validated across independent cohorts and modalities. All predictive claims
  are framed as hypothesis-generating. Wet-lab validation is beyond scope and named as
  the primary avenue for follow-up.

**A2. "Prognostic gain over staging is small (ΔC ≈ 0.005) — clinical value?"**
- Status: **Pre-empted in text.**
- Where addressed: Results §3.7 (ΔC +0.0045, LRT p=0.0020, IDI 0.021, NRI 0.175 out-of-sample in ACRG); Clinical-implications callout in Discussion; Limitations §4.1 item (2).
- Response stance: We explicitly state the signature is a research instrument, not a
  clinical tool, and that it does not change management today. The honest ΔC is reported,
  not hidden.

**A3. "Two arms (microbiome, MR) are null — do they belong?"**
- Status: **Framed as boundary/cautionary evidence.**
- Where addressed: Abstract; §3.8–3.9; Limitations §4.1 items (4)–(5).
- Response stance: The null MR (power-limited European instruments) and the
  cautionary microbiome arm (batch-confounded discovery; only reduced α-diversity
  richness generalises) are reported as boundary evidence that constrains, rather than
  inflates, the claims. Fallback if a reviewer insists: both arms can be moved fully to
  Supplementary without affecting the central stromal finding.

**A4. "Signature stability is weak — only SERPINE1 robustly selected."**
- Status: **Disclosed; reframed to the module.**
- Where addressed: §4.1 item (2) — bootstrap stability selection retained SERPINE1
  robustly (13/25 genes >50% of resamples, none >80%).
- Response stance: The robust, externally-preserved finding is the co-expression
  *module* (preservation Zsummary 15.9–17.1), not any individual gene; the gene list is
  explicitly not over-interpreted.

---

## B. Likely comments (tier-dependent, ~50–70%)

**B1. Batch scrutiny on the integrated TCGA+GTEx DEG (genomic inflation λ=17.3).**
- Status: **Strongly armed.**
- Anchors (verified): permutation null collapses to λ mean **1.061** (95th pct 1.575,
  max 2.318, B=100; `results/tables/DEG_permutation_null_lambda.csv`); TCGA-only
  replication 100% of top 100–200, ≥98.6% of top 500, 99.5% overall concordance
  (`results/tables/DEG_TCGAonly_replication.csv`).
- Where addressed: §3.1 (full paragraph).
- Response stance: λ here is not a QC statistic in the GWAS sense; it is expected when the
  global null is grossly false. The permutation null and GTEx-free TCGA-only replication
  jointly show the signal is neither batch- nor GTEx-driven.

**B2. "Microbiome discovery cohort confounds tumour/normal with batch."**
- Status: **Fully analysed and closed** — this is the most batch-controlled arm of the paper.
- Anchors (verified):
  - Confound is near-total: GCT libraries sit almost entirely on flowcells carrying no other
    phenotype (of 216 tumour/adjacent-normal pairs, only 1 shares a flowcell;
    `results/microbiome_biomarker/01_confound_crosstab_final.csv`).
  - Compositional (closure-robust) analysis already done: **Aitchison/CLR** β-diversity +
    CLR-based differential abundance (§2.8).
  - Batch-adjusted PERMANOVA already computed: phenotype R² collapses under flowcell
    adjustment — Bray 0.065→0.011, Aitchison 0.049→0.016
    (`results/microbiome_biomarker/02_beta_permanova.csv`).
  - Batch sanity: classifier AUC 0.916 but the same features predict flowcell at 78% vs 55%
    baseline (`05_rf_metrics_and_batch_sanity.csv`).
  - Clean-cohort non-replication: PRJNA641258 PERMANOVA p=0.80, 0/344 genera significant.
- Where addressed: §2.8 (methods), §3.8 (results), §4.1 item (4), Discussion §4.
- Response stance: No additional analysis required — compositional methods, batch
  adjustment, batch sanity check, and independent clean-cohort replication are all present.
  A reviewer asking for any of these can be pointed directly to the above.

**B3. "GSE84437 validation failed — explain."**
- Status: **Explained + stratified sensitivity done.**
- Anchors (verified, `results/validation_multi/GSE84437_Tstage_stratified.csv`):
  C-index <0.5 in every T-stage stratum — All 0.465, Early(T1–T3) 0.444, Late(T4) 0.483,
  T2–T3 0.442.
- Where addressed: §3.4 tested-negative paragraph; Supplementary Table S6.
- Response stance: The null is attributed to platform (microarray) and stage composition
  (~67% pT4); T-stage stratification does not rescue it, and this is reported as a
  tested-negative result rather than omitted.

**B4. "Nomogram requires TMB, not routinely available."**
- Status: **Disclosed deployability limitation.**
- Where addressed: §4.1 item (3).

---

## C. Presentation / formatting (journal-specific, do at submission)

- **C1. Reference style & section order** — reformat to the chosen journal's template.
- **C2. Methods detail** — cohort sizes, scRNA QC (nFeature 200–6000, %mt<20, 2000 HVGs,
  30 PCs, resolution 0.5), and MR instrument selection (r²<0.001 / 10,000 kb, p<5e-8→1e-5)
  are now stated in §2.5, §2.9, §2.10 (added this revision).
- **C3. Graphical abstract** — produce if the target journal requests one.
- **C4. English-language polish** — a careful copyedit pass pre-empts the standard comment.

---

## D. Reviewer-proofing summary

The paper's defining feature — it concedes and pre-answers its own weaknesses — means the
usual reviewer ammunition (leakage-inflated C, hidden nulls, unstated batch risk) is
already spent. Most anticipated comments resolve to "point to section X," and the few
optional analyses (B2 microbiome sensitivity) are non-blocking. Expected outcome at a
well-matched specialty journal: **minor-to-moderate revision**, not a re-analysis.

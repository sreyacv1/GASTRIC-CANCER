# DEFECT REPORT — Gastric-Cancer Multi-Omics Manuscript (current-state audit)

**Scope.** Independent, current-state audit of `PAPER.md` (the revised manuscript) against
its `results/` source files, the two prior internal reviews (`VERIFICATION_AUDIT.md`;
`internal/REVIEW_PANEL_v5.md`), and the analysis code. Each item is classified
**FIXED** (already correct in the current `PAPER.md`), **OPEN — closed this cycle**
(a genuine defect I corrected with backing data), **OPEN — needs author** (cannot be
derived from the project), or **OPEN — optional rerun** (a strengthening item that needs
re-execution against large primary data; flagged for the authors' decision, not silently
skipped).

**Overall assessment.** The manuscript is numerically accurate and methodologically
above the norm for its target tier (IF 3–6 translational/oncology, e.g. *Journal of
Translational Medicine*, *BMC Cancer*, *Cancers*). Independent re-computation of headline
numbers from the result files matched to the reported precision in every case checked (see
§1). The remaining defects were **reporting/completeness** issues, not data or integrity
problems; the fixable ones are closed in this cycle with backing analyses now on disk.

---

## 1. Headline-number verification (independent re-computation)

Re-computed directly from the `results/` CSVs this session; all matched:

| Claim (§) | Manuscript value | Source file | Recomputed | Status |
|---|---|---|---|---|
| Nested-CV discrimination (§3.4) | Harrell C 0.611 (0.562–0.659); Uno 0.573 | `nested_cv/performance.csv` | 0.611 / 0.573 | ✓ |
| Hartung–Knapp meta (§3.4) | pooled HR 1.19 (0.96–1.47), p=0.073, I²=19% | `meta_HK/meta_result.csv` | 1.188 / 0.962–1.466 / 0.0726 / 19.2 | ✓ |
| Module preservation (§3.6) | Zsummary 15.9 / 16.8 / 17.1 | `module_preservation/preservation_summary_RED.csv` | 15.85 / 16.80 / 17.08 | ✓ |
| MR IVW, *H. pylori* (§3.9) | OR 0.96 (0.71–1.30), p=0.79, 17 SNP | `mr_real/MR_results_all_methods_REAL.csv` | 0.959 / 0.707–1.300 / 0.786 / 17 | ✓ |
| Immune validation (§3.3) | T-cell vs leukocyte ρ=0.67, p=3.6×10⁻³⁶ | `immune/validation_vs_measured.csv` | 0.6656 / 3.57×10⁻³⁶ | ✓ |
| scRNA composition (§3.10) | 43,992 cells; fibroblast 4.16% | `scrna/celltype_composition.csv` | 43,992 / 4.16% | ✓ |
| External utility (§3.7) | clinical C 0.711, combined 0.716, ΔC +0.005 | `external_utility_ACRG/*` | 0.711 / 0.716 / +0.0045 | ✓ |

No fabrication signature, no identical-statistic duplication, no unsourced headline number
found. This corroborates the prior audit's clean integrity verdict.

---

## 2. Defects already FIXED in the current PAPER.md

Verified present/correct in the current text (fixed by the prior audit cycle):

1. **ACRG HR table cell** — now **1.90 (1.37–2.62)** (unadjusted median-split, matching its
   own log-rank p and `cindex_HR_summary.csv`); the stage/age-adjusted 1.76 is retained only
   where the text explicitly describes the adjusted model. [was VERIFICATION_AUDIT §3.1 must-fix]
2. **Flowcell-prediction accuracy** — now **78%** consistently in §3.8 and the Discussion
   (file 0.776). [was §3.2 should-fix]
3. **MR F-statistic range** — now **19.3–20.3** (matches the per-exposure log), correcting the
   v5 panel's flagged "19.3–19.7". ✓ re-verified against `scratch/mr_real.log`.
4. **MR-Egger / *Fusobacterium* contradiction** — the current text states five of six exposures
   show no pleiotropy and *Fusobacterium* carries a **nominal significant** Egger intercept
   (0.038, p=0.040) that is disclosed rather than folded into a "no pleiotropy for any" summary.
5. **GTEx-confound concordance trace** — `results/tables/DEG_TCGAonly_replication.csv` now exists
   and backs the "≥98.6% top-500 / 99.5% across 6,952 genes" claim (was v5 "untraceable").
6. **Signature selection-stability** — `results/signature_stability/` now exists (B=200 bootstrap,
   per-gene selection frequency); §4.1 reports it (SERPINE1 retained; 13/25 genes >50%).
7. **Trait-naming** — MR exposures verified against raw `gwasinfo()` output (`MR_gwasinfo_RAW.md`);
   "Prevotella 9" / "Fusobacterium A" identities confirmed.

---

## 3. Defects OPEN at start of this cycle — CLOSED with backing data

Each corrected in `PAPER.md` this cycle; backing analysis written to `results/`.

### 3.1 MR threshold wording misstated which exposures needed relaxation — **[CLOSED]**
- **Defect (v5, major).** §2.9/§3.9 implied only the microbial genera needed the relaxed
  p<10⁻⁵ threshold, with *H. pylori* at genome-wide significance. The pipeline log shows
  **all six exposures**, including *H. pylori* (21 instruments at p<10⁻⁵), used the locus-wide
  threshold; none reached ≥3 genome-wide-significant SNPs.
- **Fix.** Text now states plainly that every exposure relied on locus-wide-suggestive
  (p<10⁻⁵) instruments. New supplementary table `results/mr_real/MR_per_exposure_instruments_REAL.csv`
  reports per-exposure threshold, pre-clump nSNP, harmonised nSNP-used, mean/min F, IVW OR/CI/p.

### 3.2 MR-PRESSO promised in Methods but never reported — **[CLOSED]**
- **Defect (v5, major).** MR-PRESSO named in §2.9 but results absent from Results/Discussion.
- **Backing.** Extracted global-test p from all six `presso_*.rds`:
  *H. pylori* **RSSobs 32.6, p=0.052** (borderline, the only one near 0.05, converging with its
  significant Cochran's Q p=0.035); *Streptococcus* 0.159; *Lactobacillus* 0.542; *Prevotella*
  0.605; *Fusobacterium* 0.611; *Veillonella* 0.676. Written to
  `results/mr_real/MR_PRESSO_global_REAL.csv`.
- **Fix.** §3.9 now reports the six PRESSO global p-values and notes *H. pylori*'s null should be
  read with marginally more caution given the convergent Q/PRESSO signal (no single outlier SNP
  flagged).

### 3.3 F-statistics computed pre-harmonisation — **[CLOSED, disclosed]**
- **Defect (v5, moderate).** Reported F characterises the fetched instrument set, not the harmonised
  SNPs entering the model (e.g. *H. pylori* 21→17, *Veillonella* 11→8).
- **Fix.** §3.9 now states the F-statistics are pre-harmonisation; the supplementary table reports
  both pre-clump and harmonised-used nSNP so the drop is explicit.

### 3.4 Alpha-diversity reported by p-value only; Simpson direction discordance omitted — **[CLOSED]**
- **Defect (v5, moderate; two reviewers).** §3.8 reported the diversity decline by p-value alone;
  with n≈900 large p-values are near-guaranteed, and Simpson (evenness) does not track Shannon.
- **Backing.** Effect sizes across the non-confounded cascade endpoints (`02_alpha_effectsizes_cascade.csv`):
  Observed richness Cliff's δ **−0.27** (p=3.3×10⁻⁷), Shannon **−0.12** (p=0.02), Simpson **−0.07**
  (p=0.17, **not significant**). The ordered Jonckheere–Terpstra trend is significant for
  richness/Shannon but not Simpson (`02_alpha_trend_tests.csv`).
- **Fix.** §3.8 now reports the median differences and Cliff's δ, states that the consistently
  replicated signal is **reduced richness** (not evenness), and notes the Simpson discordance
  explicitly — the same standard of transparency the paper applies to the DEG λ.

### 3.5 Single-cell "CAF programme" — partial circularity — **[CLOSED, non-circular re-analysis added]**
- **Defect (v5, major).** The fibroblast cluster was annotated with markers
  {DCN, LUM, COL1A1, COL1A2, PDGFRB, FAP, COL3A1}; six of these are among the 29 hub genes then
  reported as fibroblast-dominant — guaranteed by construction for that subset.
- **Backing.** Re-ran the localisation on the **23 hub genes not used for annotation**
  (`gene_dominant_celltype_noncircular.csv`): **all 23/23 remain fibroblast-dominant**, median
  fraction-in-dominant 0.96. The single non-fibroblast gene in the original 28/29 (FAP →
  endothelial) was itself one of the removed annotation markers.
- **Fix.** §3.10 now discloses the six-gene annotation overlap, reports the non-circular 23/23
  result as the primary evidence, adds the dissociation-bias caveat (enzymatic protocols
  under-recover fibroblasts/CAFs) and softens "CAF programme" toward "fibroblast/stromal
  programme" where activation status is not shown.

### 3.6 Structural / framing / discoverability — **[CLOSED]**
- Missing `## 4. Discussion` header (multiple reviewers) — **added**.
- No clinician-facing "Clinical implications" statement — **added** (signature/nomogram are
  research instruments, not validated clinical tools; AJCC staging remains standard of care).
- No TRIPOD statement for a prediction-model paper — **added** Methods line + supplementary
  checklist (`TRIPOD_checklist.md`).
- Keywords lacked CAF/stroma despite it being the "principal finding" — **added**.
- MiBioGen-gut vs gastric-tissue *Streptococcus* disclaimer in the Discussion — **added**.
- R-version string inconsistency (4.3.3 vs 4.3.x) and WGCNA soft-power range (0.877 vs 0.865
  at power 11) — **reconciled**.

---

## 4. Defects OPEN — require the authors (cannot be derived)

1. **Author identity** — full names, ORCIDs, corresponding-author name/email/postal address,
   and Authors' contributions. Left as clearly-marked `[PLACEHOLDER]` fields (line 3–4, 180).
2. **Code-repository URL / archival DOI** — the Data/Code-availability statements point to
   "[REPO URL on publication]". A Zenodo (or equivalent) DOI'd snapshot should be minted at
   submission and the URL inserted.

## 5. Defects OPEN — optional strengthening (need a rerun against primary data)

These are legitimate reviewer asks that require re-executing pipeline stages against large
primary downloads; they strengthen but are **not** blocking, and each has an honest caveat
already in the text. Flagged for the authors to decide before submission:

1. **GSE84437 T-stage stratification** (v5, major — cancer-genomics). Platform (Illumina) is
   confounded with the stage-restriction explanation for the one non-validating cohort. A within-
   cohort T1–T3 vs T4 C-index/HR split would separate the two. *Current mitigation:* the text
   flags the platform confound and hedges the range-restriction account as "well-supported
   explanation rather than proven cause."
2. **ORA on a more discriminating gene set** (v5, moderate — biostatistics). The GO/KEGG ORA runs
   on the permissive ~3,700–4,000-gene lists. *Current mitigation:* rank-based GSEA is foregrounded
   as primary; the fix here is to demote ORA to a qualitative consistency check in prose (a text
   edit — applied this cycle) or rerun on a top-|t| set (a rerun — deferred).
3. **CAF subclustering** (v5, major — single-cell). True myCAF/iCAF subclustering would convert the
   lineage check into an activation-state test. *Current mitigation:* framing softened to
   fibroblast/stromal programme; deferred as it needs re-processing GSE134520.
4. **Ulcer-vs-non-ulcer sensitivity within the Normal reference** (v5, major — microbiome). The
   near-total compositional collapse in ulcer tissue is now surfaced in the text as a heterogeneous-
   reference caveat; a full paired-subset re-run of the tumour-vs-normal contrast is the deferred
   strengthening step.
5. **Compositional-closure robustness for the *H. pylori*–*Streptococcus* co-abundance** (v5,
   major — microbiome): recompute with a closure-robust metric (proportionality ρ / SparCC).
   *Current mitigation:* a compositional-closure alternative-explanation sentence added to §3.8.
6. **Repository/environment lock** (v5, major — reproducibility): commit untracked integrated-DEG
   scripts/outputs, add renv.lock/sessionInfo, account for the orphan MaAsLin2 folder. Partially
   addressed this cycle (PIPELINE.md MaAsLin2 note, sessionInfo dump, R-version reconciliation);
   the git commit/tag + environment lock is an author action on the repo.

---

## 6. Bottom line

The manuscript's science is sound and its numbers reproduce. Every backed reporting defect
identified across the prior audit and the v5 panel is closed in this cycle, with the supporting
analyses saved under `results/`. What remains is (a) two fields only the authors can supply
(identity, repo DOI) and (b) a short menu of optional reruns that would further strengthen —
but do not block — a submission at the target tier.

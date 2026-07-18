# Submission Checklist — Gastric-Cancer Multi-Omics Manuscript

Status after this revision cycle. Items are grouped by who must act.

## A. Done in this cycle (no further action)

- [x] All headline numbers independently re-verified against `results/` (see `DEFECT_REPORT.md` §1).
- [x] MR threshold wording corrected: all six exposures at p<1×10⁻⁵ (Methods §2.9, Results §3.9, Abstract).
- [x] MR-PRESSO global tests reported (H. pylori borderline p=0.052); Supplementary Table S2.
- [x] Per-exposure MR instrument table (pre/post-harmonisation nSNP, F): Supplementary Table S1.
- [x] F-statistics disclosed as pre-harmonisation.
- [x] Alpha-diversity effect sizes + Simpson discordance reported (§3.8); Supplementary Table S4.
- [x] Single-cell circularity disclosed; non-circular 23/23 re-analysis added (§3.10); Supplementary Table S5.
- [x] Dissociation-bias and MiBioGen-vs-tissue caveats added.
- [x] `## 4. Discussion` header added; drug-repurposing arm tied into narrative.
- [x] Clinical-implications callout added (not a clinical tool; AJCC standard of care).
- [x] TRIPOD statement (Methods §2.7) + completed checklist (Supplementary Table S3).
- [x] Keywords expanded (CAF, tumour stroma).
- [x] Reference list expanded to full method/biology citations (39 refs, Vancouver style).
- [x] R-version reconciled to 4.3.3 across manuscript, PIPELINE.md, README; sessionInfo dump referenced.
- [x] WGCNA scale-free R² range corrected (0.865–0.886).
- [x] MaAsLin2 orphan folder documented in PIPELINE.md.
- [x] Exposure trait names tightened (*Prevotella 9*, *Fusobacterium* A).

## B. Requires the authors (cannot be derived from the project)

- [ ] **Author names + ORCIDs** — fill `⟦PLACEHOLDER — AUTHOR NAMES + ORCIDs⟧` (line 3).
- [ ] **Corresponding author** — name, institutional email, postal address (line 4).
- [ ] **Authors' contributions** — assign CRediT roles (Declarations).
- [ ] **Repository URL + archival DOI** — mint a version-tagged Zenodo (or equivalent) snapshot
      of the repo at submission and insert the DOI (Code availability).
- [ ] **Funding** — confirm the "no specific grant" statement or edit.
- [ ] **Cover letter + journal choice** — the paper is framed for a BMC-family translational/
      oncology venue (e.g. *Journal of Translational Medicine*, *BMC Cancer*, *Cancers*).

## C. Repository actions before archival (author, on the git repo)

- [ ] Commit the currently-uncommitted `PAPER.md`, the integrated-DEG scripts and the new
      supplementary tables under `results/`.
- [ ] Confirm `renv.lock` (or the pinned `package_versions.csv` + `sessionInfo.txt`) fully
      captures the environment; the reviewer asked for an environment lock.
- [ ] Tag the release that the DOI will point to.

## D. Optional strengthening reruns (author decision — NOT blocking; see DEFECT_REPORT §5)

These are legitimate reviewer asks that each need a rerun against large primary data. The
manuscript already carries an honest caveat for every one; deciding to run them is a
cost/benefit call, not a correctness requirement.

- [ ] GSE84437 within-cohort T-stage stratification (separates platform from stage-restriction).
- [ ] ORA rerun on a top-|t| gene set (or demote ORA to qualitative — text already foregrounds GSEA).
- [ ] CAF myCAF/iCAF subclustering in GSE134520 (activation-state test).
- [ ] Ulcer-vs-non-ulcer paired-subset sensitivity for the tumour-vs-normal contrast.
- [ ] Closure-robust (proportionality ρ / SparCC) recomputation of any tissue co-occurrence claim.

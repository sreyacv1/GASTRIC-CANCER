# Review Resolution — 2026-07-23

Resolution of three review findings against `PAPER.md` / `FIGURE_SOURCES.md`.
All changes are working-tree edits (not yet committed).

## 1. Data-availability accessions PRJNA641258 / PRJNA413125 — *false positive, resolved*

**Finding (fail):** tissue-16S validation accessions flagged as fabricated —
no traceable origin in session history.

**Why the rebuttal didn't close it:** "already in git / used across scripts"
proves the accessions are *pre-existing*, not *real*. Pre-existing fabrication
is still fabrication, so the finding correctly stayed open. Provenance is
verified against the source databases, not session history.

**Verification (NCBI):**
- **PRJNA641258** — real BioProject, CNR Italy (reg. 2020-06-23), 40 SRA
  experiments, stomach microbiome 16S. Associated paper **Ravegnini et al.,
  *Int J Mol Sci* 2020, PMC7766162**.
- **PRJNA413125** — real BioProject, IPATIMUP (Porto), gastric microbial
  community profiling, 135 experiments — the Ferreira/Portugal group (2018).
- Repo corroboration: `analysis/28_validation_IT.R`, `29_validation_PT.R`,
  and `results/microbiome_biomarker/validation_cohort_search.md` record the
  real per-run accessions (`SRR12072987`, `SRR6151146`).

**Fix applied:** the one genuine error was a wrong first author — the
data-availability line credited *Palmieri et al.* for PMC7766162, which is
**Ravegnini et al.** Corrected in `PAPER.md`.

## 2. Broken figure embeds — *fixed (all 23, not just Figure 1)*

**Finding (fail):** Figure 1 embed `{{artifact:art_34ca7531-…}}` is malformed
and won't render.

**Root cause (broader than reported):** *all 23* figure embeds used claude.ai
artifact-store syntax `{{artifact:art_<uuid>}}`, which renders as literal text
in a git-hosted markdown file — none displayed on GitHub or any standard viewer.

**Fix applied:** rewrote every embed to a real repo-relative path, preferring
PNG over PDF so they render inline, and using the committed composites for
S15–S19. All 23 targets confirmed present on disk.

| Figures | Target |
|---|---|
| Fig 1–4 | `results/figures/Fig{1..4}.png` |
| S1–S14 | individual panel PNGs (per caption `File:` line) |
| S15–S19 | `results/composite_figures/s1{5..9}_*.png` |

## 3. FIGURE_SOURCES.md overreach — *reconciled*

**Finding (warn):** header claimed "every panel git-tracked … complete check,"
but only 13 of ~44 files were actually verified that session.

**Fix applied:**
- Ran the complete `git ls-files` check the document claimed: **every** source
  panel and all 5 composites are present and tracked — **0 missing, 0
  untracked** — so the statement is now literally true.
- Corrected two real path errors in the S18/S19 rows: shorthand filenames
  (`scatter_H_pylori`, `Streptococcus`) that don't exist on disk → actual names
  (`scatter_H__pylori_IgG_seropositivity`, `scatter_Streptococcus__genus_`).
- Added the composite display-file paths (`results/composite_figures/…`).

## Files changed
- `PAPER.md` — 1 author fix + 23 figure embeds rewritten
- `FIGURE_SOURCES.md` — verification statement + S15–S19 source rows

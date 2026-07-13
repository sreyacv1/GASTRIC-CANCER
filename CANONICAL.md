# CANONICAL PIPELINE — Gastric Cancer Multi-Omics v2.0

> **Last updated:** 2026-05-21 (post code-review patch)

## ▶ Run In This Order

```bash
# Step 1: Activate R environment
source r_env/bin/activate  # or: conda activate ./r_env

# Step 2: Run Part 1 (Steps 00–05: Data, QC, Batch Correction)
# Expected runtime: 45–90 min (mostly data I/O)
Rscript -e 'source("gastric_cancer_multiomics_v2_part1.R")'

# Step 3: Run Part 2 (Steps 06–12: Microbiome, Pathways, WGCNA, MR)
# Expected runtime: 60–120 min (WGCNA + DIABLO + MR)
Rscript -e 'source("gastric_cancer_multiomics_v2_part2.R")'
```

### Canonical Files

| File | Role |
|---|---|
| `gastric_cancer_multiomics_v2_part1.R` | **CANONICAL Part 1** — Steps 00–05 |
| `gastric_cancer_multiomics_v2_part2.R` | **CANONICAL Part 2** — Steps 06–12 |
| `README.md` | Full pipeline documentation |
| `GEMINI.md` | Engineering standards (read before editing) |

### All other `.R` files in this directory are ARCHIVED
They are preserved for reference but must **not** be used for publication analyses.
See `archive/` subdirectory.

---

## ⚠️ Before Running: Required Data

Ensure these files exist before running Part 2 (microbiome steps will hard-stop without them):

```
data/host/tcga_stad_rse.rds           ← TCGA-STAD (user-uploaded)
data/host/GTEx_v10_tpm.gct.gz         ← GTEx v10 (or gtex_stomach_tpm_log2.csv.gz)
data/microbiome/otu_table.csv         ← DDBJ PRJDB20660 / PRJNA830774
data/microbiome/taxonomy.tsv
data/microbiome/metadata_microbiome.tsv
```

> To use **demonstration data** for dry-run testing ONLY (NOT for publication):
> ```r
> options(allow_synthetic_microbiome = TRUE)
> source("gastric_cancer_multiomics_v2_part2.R")
> ```

---

## 🔄 Re-running a specific step

Each step saves a `.done` sentinel in `results/rdata/`. Delete it to force re-run:

```bash
# Example: re-run WGCNA
rm results/rdata/step08_wgcna.done
Rscript -e 'source("gastric_cancer_multiomics_v2_part2.R")'
```

---

## 📋 Key Parameters (v2.0 post code-review)

| Parameter | Value | Rationale |
|---|---|---|
| WGCNA minModuleSize | 15 | Relaxed from 30 to avoid 0-module collapse |
| WGCNA deepSplit | 1 | Reduced from 3 to avoid over-splitting |
| WGCNA mergeCutHeight | 0.30 | Relaxed from 0.20 to merge similar modules |
| WGCNA corType | bicor | Robust to outliers (unchanged) |
| ComBat mean.only | TRUE | Data is Z-scored; adjust means only (GEMINI.md) |
| MaAsLin2 shared samples | ≥20 by ID | Position-matching disabled (creates false associations) |
| GSVA API | GsvaParam() | Required for GSVA ≥1.50 (Bioconductor 3.17+) |
| PERMANOVA | ~ status + sex + age + hp_status | Marginal test with confounders |

---

*Pipeline v2.0 | Gastric Cancer Microbiome–Transcriptome Multi-Omics Study*

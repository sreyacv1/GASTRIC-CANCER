#!/usr/bin/env Rscript
# 26_depmap_dependency.R
# DepMap cross-reference for the in-silico drug-repurposing targets.
#
# The connectivity-map repurposing (16/21) nominated three target classes:
#   PI3K/mTOR : PIK3CA, MTOR, PIK3CB   (buparlisib/alpelisib/everolimus/sirolimus)
#   CDK4/6    : CDK4, CDK6             (palbociclib, rank 4)
#   FGFR      : FGFR1, FGFR2, FGFR3, FGFR4 (dovitinib, rank 5)
# Question: are gastric-cancer cell lines actually DEPENDENT on these genes?
#
# Data: DepMap 24Q4 Public CRISPR (Chronos) gene effect + Model annotations
#   (figshare article 27993248: CRISPRGeneEffect.csv=51064667, Model.csv=51065297)
# Chronos gene effect: more negative = more essential; ~ -1 = median common
#   essential; > ~ -0.5 typically not dependent. "Dependent" threshold: < -0.5.
#
# HARD RULE: real data only. Every number from the downloaded DepMap matrix.

suppressPackageStartupMessages({ library(data.table) })
set.seed(42)  # no randomness here; set for repo consistency

ge_file <- "scratch/depmap_tmp/CRISPRGeneEffect.csv"
md_file <- "scratch/depmap_tmp/Model.csv"
stopifnot(file.exists(ge_file), file.exists(md_file))

outdir <- "results/depmap"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

targets <- c("PIK3CA","MTOR","PIK3CB","CDK4","CDK6",
             "FGFR1","FGFR2","FGFR3","FGFR4")
target_class <- c(PIK3CA="PI3K/mTOR", MTOR="PI3K/mTOR", PIK3CB="PI3K/mTOR",
                  CDK4="CDK4/6", CDK6="CDK4/6",
                  FGFR1="FGFR", FGFR2="FGFR", FGFR3="FGFR", FGFR4="FGFR")

## 1. Map target symbols -> column POSITIONS (robust to " (entrez)" suffix) -----
hdr <- strsplit(readLines(ge_file, n = 1), ",")[[1]]
sym <- sub(" \\(.*$", "", hdr)                 # strip " (entrez)"
idx <- match(targets, sym)                     # column position of each target
missing <- targets[is.na(idx)]
if (length(missing)) cat("WARNING: target not in DepMap:", missing, "\n")
use_targets <- targets[!is.na(idx)]
sel_idx <- c(1L, idx[!is.na(idx)])             # col 1 = ModelID + target cols

## 2. Read only the needed columns (Chronos gene effect) -----------------------
ge <- fread(ge_file, select = sel_idx)
setnames(ge, c("ModelID", use_targets))
cat(sprintf("DepMap cell lines with CRISPR data: %d\n", nrow(ge)))

## 3. Model annotations -> gastric/stomach classification ----------------------
md <- fread(md_file)
# Gastric/stomach = Oncotree subtype containing "Stomach" (STAD and stomach
# histologies) within the Esophagus/Stomach lineage. Reported alongside the
# broader Esophagus/Stomach lineage for context.
md[, is_gastric := grepl("Stomach", OncotreeSubtype, ignore.case = TRUE)]
md[, is_esoph_stomach := OncotreeLineage == "Esophagus/Stomach"]

ge <- merge(ge, md[, .(ModelID, OncotreeLineage, OncotreePrimaryDisease,
                       OncotreeSubtype, is_gastric, is_esoph_stomach)],
            by = "ModelID", all.x = TRUE)

n_gastric <- sum(ge$is_gastric, na.rm = TRUE)
n_es      <- sum(ge$is_esoph_stomach, na.rm = TRUE)
n_pan     <- nrow(ge)
cat(sprintf("Gastric/stomach lines: %d | Esophagus/Stomach lineage: %d | pan-cancer: %d\n",
            n_gastric, n_es, n_pan))

## 4. Per-gene dependency: gastric mean vs pan-cancer mean ---------------------
summ <- rbindlist(lapply(use_targets, function(g) {
  v_all <- ge[[g]]
  v_gas <- ge[is_gastric == TRUE][[g]]
  data.table(
    gene              = g,
    target_class      = target_class[[g]],
    n_gastric         = sum(!is.na(v_gas)),
    gastric_mean_CERES= mean(v_gas, na.rm = TRUE),
    gastric_median    = median(v_gas, na.rm = TRUE),
    gastric_min       = min(v_gas, na.rm = TRUE),
    gastric_frac_dep  = mean(v_gas < -0.5, na.rm = TRUE),   # < -0.5 ~ dependent
    pancancer_mean    = mean(v_all, na.rm = TRUE),
    pancancer_frac_dep= mean(v_all < -0.5, na.rm = TRUE),
    delta_gastric_minus_pan = mean(v_gas, na.rm = TRUE) - mean(v_all, na.rm = TRUE))
}))
# selective dependency flag: gastric mean < -0.5 AND more essential than pan
summ[, gastric_dependent := gastric_mean_CERES < -0.5]
summ[, gastric_selective := gastric_dependent & (delta_gastric_minus_pan < -0.1)]

cat("\n== Gastric CRISPR dependency (Chronos gene effect) ==\n")
print(summ[, .(gene, target_class, n_gastric,
               gastric_mean_CERES = round(gastric_mean_CERES, 3),
               pancancer_mean = round(pancancer_mean, 3),
               gastric_frac_dep = round(gastric_frac_dep, 3),
               gastric_dependent, gastric_selective)])

fwrite(summ, file.path(outdir, "gastric_dependency.csv"))

## 5. SUMMARY ------------------------------------------------------------------
dep_genes <- summ[gastric_dependent == TRUE, gene]
sel_genes <- summ[gastric_selective == TRUE, gene]

sink(file.path(outdir, "SUMMARY.txt"))
cat("DepMap 24Q4 Public CRISPR (Chronos) dependency for drug-repurposing targets\n")
cat("=========================================================================\n")
cat(sprintf("Gastric/stomach cell lines: n=%d (Oncotree subtype ~ 'Stomach')\n", n_gastric))
cat(sprintf("Esophagus/Stomach lineage: n=%d | pan-cancer lines: n=%d\n\n", n_es, n_pan))
cat("Chronos gene effect: more negative = more essential; < -0.5 ~ dependent;\n")
cat("~ -1.0 = median common-essential gene.\n\n")
cat("Per-target (gastric mean vs pan-cancer mean):\n")
print(summ[, .(gene, target_class,
               gastric_mean = round(gastric_mean_CERES, 3),
               pancancer_mean = round(pancancer_mean, 3),
               gastric_frac_dep_lt_neg0.5 = round(gastric_frac_dep, 3),
               dependent = gastric_dependent, selective = gastric_selective)])
cat("\n== VERDICT ==\n")
if (length(dep_genes))
  cat("Gastric lines meet the dependency threshold (mean < -0.5) for:",
      paste(dep_genes, collapse = ", "), "\n") else
  cat("No nominated target reaches the gastric dependency threshold (mean < -0.5).\n")
if (length(sel_genes))
  cat("SELECTIVE gastric dependency (dependent AND more essential than pan-cancer):",
      paste(sel_genes, collapse = ", "), "\n") else
  cat("No nominated target shows SELECTIVE gastric dependency vs pan-cancer\n",
      "(where dependent, essentiality is pan-lineage, not gastric-specific).\n", sep = "")
sink()

cat("\nDepMap outputs written to", outdir, "\n")

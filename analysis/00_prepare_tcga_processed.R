#!/usr/bin/env Rscript
# 00_prepare_tcga_processed.R
# -----------------------------------------------------------------------------
# BRIDGE / REPRODUCIBILITY FIX.
# Regenerates results/rdata/tcga_processed.RData (objects: col_data, tcga_vst,
# res) that every downstream consumer loads (07,08,13,14,18,19, nomogram_*),
# but which NO tracked script previously created. This closes that hole so the
# pipeline is reproducible from a clean clone.
#
# PROVENANCE
#   Inputs (written by analysis/gastric_cancer_multiomics_v2_part1.R):
#     data/processed/TCGA_STAD_processed.RData  ->  tcga_raw  (raw counts,
#       60660 versioned-Ensembl x 448 samples), tcga_meta (colData + Lauren +
#       sample_type_simple), gene_map (Ensembl->symbol fallback annotation).
#   part1 itself reads data/host/tcga_stad_rse.rds (the GDC RangedSummarized-
#   Experiment); we re-read it ONLY for rowData()$gene_name, the symbol source
#   that best matches the committed tcga_vst. If absent we fall back to gene_map.
#
#   PREREQUISITE: data/processed/TCGA_STAD_processed.RData must exist. If it does
#   not, run analysis/gastric_cancer_multiomics_v2_part1.R first (heavy: needs
#   data/host/tcga_stad_rse.rds and the GTEx/GEO raw inputs). This script does
#   NOT trigger that download/preprocess itself by design.
#
# METHOD (matches the schema of the committed file: DESeq2, symbol-level VST)
#   col_data  = tcga_meta with sample_type_simple renamed to `status` and
#               TCGA_subtype dropped (448 x 114, key survival/clinical columns).
#   res       = DESeq2 Tumor-vs-Normal results (baseMean, log2FoldChange, lfcSE,
#               stat, pvalue, padj, gene, sig).
#   tcga_vst  = DESeq2 varianceStabilizingTransformation, collapsed to gene
#               symbols (genes x samples, "HAT1" present).
#
# VALIDATION / SAFETY
#   The exact pre-filter and symbol annotation the original author used are not
#   recoverable from the committed inputs (no simple count filter reproduces the
#   46969 res rows; only ~89% of the committed symbols exist in the RSE
#   rowData). We therefore reconstruct a SCHEMA-COMPATIBLE file and validate it
#   STRUCTURALLY against the existing one. The real committed file is NEVER
#   overwritten: if it already exists we write results/rdata/tcga_processed_REGEN.RData
#   and print a diff; we only write the canonical path when it is missing
#   (the clean-clone case this fix is for).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})
set.seed(42)

root <- Sys.getenv("GC_ROOT", getwd())
setwd(root)

BASE      <- "data/processed/TCGA_STAD_processed.RData"
RSE       <- "data/host/tcga_stad_rse.rds"
CANON     <- "results/rdata/tcga_processed.RData"
REGEN     <- "results/rdata/tcga_processed_REGEN.RData"
dir.create("results/rdata", recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n")

if (!file.exists(BASE)) {
  stop(sprintf(paste0(
    "Prerequisite missing: %s\n",
    "  Run analysis/gastric_cancer_multiomics_v2_part1.R first to create it\n",
    "  (heavy; needs data/host/tcga_stad_rse.rds + GTEx/GEO inputs)."), BASE))
}

say("Loading base preprocessing objects from %s ...", BASE)
b <- new.env(); load(BASE, envir = b)
stopifnot(all(c("tcga_raw", "tcga_meta") %in% ls(b)))
tcga_raw  <- b$tcga_raw
tcga_meta <- b$tcga_meta
gene_map  <- if ("gene_map" %in% ls(b)) b$gene_map else NULL

# ---- col_data ---------------------------------------------------------------
# status = Tumor/Normal (from sample_type_simple); drop the two non-committed
# helper columns, keep everything else (clinical + paper_* + Lauren).
col_data <- tcga_meta
col_data$status <- col_data$sample_type_simple
drop <- c("sample_type_simple", "TCGA_subtype")
col_data <- col_data[, setdiff(colnames(col_data), drop), drop = FALSE]
# place `status` just before `Lauren` to mirror the committed column order
if (all(c("status", "Lauren") %in% colnames(col_data))) {
  ord <- setdiff(colnames(col_data), c("status", "Lauren"))
  col_data <- col_data[, c(ord, "status", "Lauren"), drop = FALSE]
}
say("col_data: %d samples x %d columns", nrow(col_data), ncol(col_data))

# ---- DESeq2 Tumor vs Normal -------------------------------------------------
cnt <- round(as.matrix(tcga_raw))
mode(cnt) <- "integer"
stopifnot(identical(colnames(cnt), rownames(col_data)))

cd <- col_data
cd$status <- factor(cd$status, levels = c("Normal", "Tumor"))
keep_samp <- !is.na(cd$status)
cnt <- cnt[, keep_samp]; cd <- cd[keep_samp, ]

say("Building DESeqDataSet (%d genes x %d samples) ...", nrow(cnt), ncol(cnt))
dds <- DESeqDataSetFromMatrix(countData = cnt, colData = cd, design = ~ status)
# Loose pre-filter (exact original threshold unknown; independent filtering in
# results() removes remaining low-information genes via padj = NA).
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
say("After pre-filter (rowSums >= 10): %d genes", nrow(dds))

say("Running DESeq() ... (this is the heavy step)")
dds <- DESeq(dds)

res_raw <- results(dds, contrast = c("status", "Tumor", "Normal"))
res <- as.data.frame(res_raw)
res <- res[!is.na(res$padj), ]
res$gene <- rownames(res)
res$sig  <- ifelse(res$padj < 0.05 & res$log2FoldChange >  1, "Up",
             ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "Down", "NS"))
res <- res[, c("baseMean", "log2FoldChange", "lfcSE", "stat",
               "pvalue", "padj", "gene", "sig")]
res <- res[order(res$padj), ]
say("res: %d genes (Up=%d Down=%d NS=%d)", nrow(res),
    sum(res$sig == "Up"), sum(res$sig == "Down"), sum(res$sig == "NS"))

# ---- VST + Ensembl -> symbol ------------------------------------------------
say("Variance-stabilising transform ...")
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vmat <- assay(vsd)   # versioned-Ensembl x samples

# Symbol source: prefer the GDC RSE rowData(gene_name) (matches committed file
# best), else the org.Hs.eg.db-derived gene_map from part1.
sym <- NULL
if (file.exists(RSE)) {
  rse <- tryCatch(readRDS(RSE), error = function(e) NULL)
  if (!is.null(rse)) {
    rd <- rowData(rse)
    gn_col <- grep("gene_name", colnames(rd), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(gn_col)) {
      gn <- as.character(rd[[gn_col]]); names(gn) <- rownames(rse)
      sym <- gn[rownames(vmat)]
      say("Symbol mapping via RSE rowData$%s", gn_col)
    }
  }
}
if (is.null(sym) && !is.null(gene_map)) {
  key <- sub("\\..*", "", rownames(vmat))
  sym <- gene_map$hgnc_symbol[match(key, gene_map$ensembl_gene_id)]
  say("Symbol mapping via part1 gene_map (RSE unavailable)")
}
if (is.null(sym)) stop("No symbol annotation source available (RSE + gene_map both missing).")

ok <- !is.na(sym) & sym != ""
vmat <- vmat[ok, ]; sym <- sym[ok]
# collapse duplicate symbols: keep the highest-variance row per symbol
rv  <- matrixStats::rowVars(vmat)
ord <- order(rv, decreasing = TRUE)
vmat <- vmat[ord, ]; sym <- sym[ord]
keep_sym <- !duplicated(sym)
tcga_vst <- vmat[keep_sym, ]
rownames(tcga_vst) <- sym[keep_sym]
say("tcga_vst: %d symbols x %d samples | HAT1 present: %s",
    nrow(tcga_vst), ncol(tcga_vst), "HAT1" %in% rownames(tcga_vst))

# ---- structural validation against the committed file -----------------------
validated <- TRUE
if (file.exists(CANON)) {
  say("\n--- Validation vs existing %s ---", CANON)
  x <- new.env(); load(CANON, envir = x)
  chk <- function(label, cond) { say("  [%s] %s", ifelse(cond, "OK", "DIFF"), label); cond }
  v1 <- chk("object names {col_data,tcga_vst,res}",
            all(c("col_data","tcga_vst","res") %in% ls(x)))
  v2 <- chk(sprintf("col_data rows (new=%d old=%d)", nrow(col_data), nrow(x$col_data)),
            nrow(col_data) == nrow(x$col_data))
  key_cols <- c("status","vital_status","days_to_death","days_to_last_follow_up",
                "ajcc_pathologic_stage","tumor_grade","age_at_diagnosis","Lauren")
  v3 <- chk("col_data key columns present",
            all(key_cols %in% colnames(col_data)))
  v4 <- chk("tcga_vst orientation genes-in-rows / 448 samples in cols",
            ncol(tcga_vst) == nrow(col_data) &&
            identical(colnames(tcga_vst), rownames(col_data)))
  v5 <- chk("tcga_vst has HAT1 (as in committed file)",
            "HAT1" %in% rownames(tcga_vst))
  v6 <- chk("res columns match committed",
            identical(colnames(res), colnames(x$res)))
  # informational (not gating): exact dims are not recoverable
  say("  [INFO] tcga_vst genes new=%d vs old=%d ; res rows new=%d vs old=%d",
      nrow(tcga_vst), nrow(x$tcga_vst), nrow(res), nrow(x$res))
  validated <- all(v1, v2, v3, v4, v5, v6)
  say("  Structural validation: %s", ifelse(validated, "PASS", "FAIL"))
}

# ---- save -------------------------------------------------------------------
if (file.exists(CANON)) {
  save(col_data, tcga_vst, res, file = REGEN)
  say("\nExisting committed file preserved (not overwritten).")
  say("Regenerated schema-compatible objects saved to: %s", REGEN)
  if (!validated) say("WARNING: structural validation FAILED — inspect diff above.")
} else {
  save(col_data, tcga_vst, res, file = CANON)
  say("\nNo committed file present (clean clone): wrote %s", CANON)
}
say("Done.")

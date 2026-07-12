#!/usr/bin/env Rscript
# =============================================================================
#  22_deg_diagnostics.R
#  Reproducible, auditable defence of the INTEGRATED tumour-vs-normal DEG
#  (analysis/20) against the "batch not biology" reviewer critique.
#
#  Makes two ad-hoc numbers fully traceable:
#    TASK 1  TCGA-only replication of the integrated top-N up/down genes.
#    TASK 2  Permutation-null test for the median test-statistic inflation
#            (genomic-inflation factor lambda ~ 17 in the real contrast).
#
#  Discovery cohort, model and contrast are IDENTICAL to analysis/20:
#    discovery = TCGA + GTEx ; model ~ sample_type + dataset ;
#    contrast  = sample_typeTumor (Normal is the reference level).
#
#  Run:  ./r_env/bin/Rscript analysis/22_deg_diagnostics.R
#  Real data only. No fabricated numbers.
# =============================================================================
suppressPackageStartupMessages(library(limma))

OUT <- "results/tables"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# Genomic-inflation factor (median form), robust for tiny p-values.
# lambda = median(chisq_1df) / qchisq(0.5,1), chisq from upper tail of p.
lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0]
  chi <- qchisq(p, df = 1, lower.tail = FALSE)
  median(chi) / qchisq(0.5, df = 1)
}

## ===========================================================================
## TASK 1 — TCGA-only replication of the integrated top-N up / down genes.
##
##  Thresholds / matching (stated explicitly for the audit trail):
##   * Integrated significance : adj.P.Val < 0.05 ; direction = sign(moderated t).
##   * Ranking                 : top-N UP by t (desc), top-N DOWN by t (asc).
##   * Symbol matching         : exact gene symbol, integrated$gene ==
##                               TCGA$gene_symbol (both are HGNC symbols).
##   * "maps to TCGA-only DEG" : symbol present in the TCGA-only results table.
##   * "replicated"            : TCGA padj < 0.05 AND sign(log2FoldChange)
##                               equals the integrated direction.
##   * TCGA duplicate symbols  : kept the most significant row (min padj).
## ===========================================================================
intg <- read.csv(file.path(OUT, "DEG_integrated_TCGA_GTEx.csv"),
                 stringsAsFactors = FALSE)
tcga <- read.csv(file.path(OUT, "TCGA_DEG_results_symbols.csv"),
                 stringsAsFactors = FALSE)

# Dedup TCGA by symbol, keeping the most significant (smallest padj) row.
tcga <- tcga[order(is.na(tcga$padj), tcga$padj), ]
tcga <- tcga[!duplicated(tcga$gene_symbol), ]
tcga_lfc  <- setNames(tcga$log2FoldChange, tcga$gene_symbol)
tcga_padj <- setNames(tcga$padj,           tcga$gene_symbol)

sig  <- intg[intg$adj.P.Val < 0.05, ]
up   <- sig[sig$t > 0, ]; up <- up[order(-up$t), ]          # top UP by t
down <- sig[sig$t < 0, ]; down <- down[order(down$t), ]     # top DOWN by t

replic_row <- function(genes, dir_label, top_n) {
  g <- head(genes, top_n)
  n_sel    <- length(g)
  in_tcga  <- g %in% names(tcga_lfc)
  n_map    <- sum(in_tcga)
  gm       <- g[in_tcga]
  want_pos <- dir_label == "Up"
  padj_ok  <- !is.na(tcga_padj[gm]) & tcga_padj[gm] < 0.05
  dir_ok   <- if (want_pos) tcga_lfc[gm] > 0 else tcga_lfc[gm] < 0
  n_rep    <- sum(padj_ok & dir_ok)
  data.frame(
    direction               = dir_label,
    top_n                   = top_n,
    n_selected              = n_sel,
    n_mapped_TCGA           = n_map,
    pct_mapped_of_selected  = round(100 * n_map / n_sel, 1),
    n_replicated_TCGA       = n_rep,   # padj<0.05 & same direction
    pct_replicated_of_selected = round(100 * n_rep / n_sel, 1),
    pct_replicated_of_mapped   = round(100 * n_rep / n_map, 1),
    stringsAsFactors = FALSE)
}

rep_tab <- do.call(rbind, c(
  lapply(c(100, 200, 500), function(n) replic_row(up$gene,   "Up",   n)),
  lapply(c(100, 200, 500), function(n) replic_row(down$gene, "Down", n))
))

# Overall directional concordance across ALL genes significant in BOTH
# (integrated adj.P.Val<0.05 AND present in TCGA with padj<0.05).
both <- sig$gene[sig$gene %in% names(tcga_padj)]
both <- both[!is.na(tcga_padj[both]) & tcga_padj[both] < 0.05]
intg_t <- setNames(sig$t, sig$gene)
concord <- sign(intg_t[both]) == sign(tcga_lfc[both])
n_both  <- length(both); n_conc <- sum(concord)
overall <- data.frame(
  direction = "Overall_concordance_sig_in_both",
  top_n = NA_integer_, n_selected = n_both, n_mapped_TCGA = n_both,
  pct_mapped_of_selected = 100,
  n_replicated_TCGA = n_conc,
  pct_replicated_of_selected = round(100 * n_conc / n_both, 1),
  pct_replicated_of_mapped   = round(100 * n_conc / n_both, 1),
  stringsAsFactors = FALSE)

rep_tab <- rbind(rep_tab, overall)
write.csv(rep_tab, file.path(OUT, "DEG_TCGAonly_replication.csv"),
          row.names = FALSE)

cat("=== TASK 1: TCGA-only replication of integrated top-N ===\n")
cat(sprintf("Integrated sig genes (adj.P<0.05): %d up / %d down\n",
            nrow(up), nrow(down)))
cat(sprintf("TCGA-only table: %d unique symbols\n", length(tcga_lfc)))
print(rep_tab, row.names = FALSE)
cat(sprintf("Overall directional concordance (sig in both, n=%d): %.1f%%\n\n",
            n_both, 100 * n_conc / n_both))

## ===========================================================================
## TASK 2 — Permutation-null for the median test-statistic inflation lambda.
##
##  Observed: real contrast on TCGA+GTEx, model ~ sample_type + dataset.
##  Null    : permute tumour/normal labels WITHIN the TCGA stratum only
##            (GTEx stays all-normal, it has no tumours to swap). This
##            preserves the confound structure (all tumours in TCGA, the
##            TCGA-vs-GTEx mean shift, the dataset covariate, class margins)
##            while destroying the true tumour-vs-normal signal. Refit the
##            identical limma model per permutation and recompute lambda.
##
##  Reading: null lambda ~ 1 while observed ~ 17  =>  the inflation reflects
##  genuine, genome-wide tumour biology, not a batch/QC artifact. A high null
##  lambda would instead indicate residual batch/design structure.
## ===========================================================================
load("data/processed/combined_transcriptome.RData")  # combined_expr_bc, combined_meta
stopifnot(all(colnames(combined_expr_bc) == rownames(combined_meta)))
disc <- combined_meta$dataset %in% c("TCGA", "GTEx")
e  <- combined_expr_bc[, disc]
m  <- combined_meta[disc, ]
ds <- factor(m$dataset)
is_tcga <- which(m$dataset == "TCGA")

fit_lambda <- function(type_vec) {
  st  <- factor(type_vec, levels = c("Normal", "Tumor"))
  fit <- eBayes(lmFit(e, model.matrix(~ st + ds)))
  res <- topTable(fit, coef = "stTumor", number = Inf, sort.by = "none")
  lambda_gc(res$P.Value)
}

observed_lambda <- fit_lambda(m$sample_type_simple)

B <- 100  # >= 50 required
null_lambda <- numeric(B)
for (b in seq_len(B)) {
  perm <- m$sample_type_simple
  perm[is_tcga] <- sample(perm[is_tcga])   # shuffle labels within TCGA only
  null_lambda[b] <- fit_lambda(perm)
  if (b %% 10 == 0) cat(sprintf("  permutation %d/%d done\n", b, B))
}

null_mean <- mean(null_lambda)
null_p95  <- as.numeric(quantile(null_lambda, 0.95))
null_max  <- max(null_lambda)

conclusion <- if (observed_lambda > 5 && null_p95 < 2) {
  "biology: observed inflation reflects genome-wide tumour-vs-normal signal, not batch (null lambda ~1)"
} else if (null_p95 >= 2) {
  "batch/design: null lambda is also inflated, indicating residual structure"
} else {
  "inconclusive: observed not clearly separated from null"
}

null_df <- data.frame(
  observed_lambda = round(observed_lambda, 3),
  null_mean       = round(null_mean, 3),
  null_p95        = round(null_p95, 3),
  null_max        = round(null_max, 3),
  B               = B,
  conclusion      = conclusion,
  stringsAsFactors = FALSE)
write.csv(null_df, file.path(OUT, "DEG_permutation_null_lambda.csv"),
          row.names = FALSE)

cat("\n=== TASK 2: Permutation-null lambda ===\n")
cat(sprintf("Observed lambda        : %.3f\n", observed_lambda))
cat(sprintf("Null lambda mean/p95/max: %.3f / %.3f / %.3f (B=%d)\n",
            null_mean, null_p95, null_max, B))
cat(sprintf("Conclusion             : %s\n", conclusion))

cat("\n================ SAVED ================\n")
cat("results/tables/DEG_TCGAonly_replication.csv\n")
cat("results/tables/DEG_permutation_null_lambda.csv\n")

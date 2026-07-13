#!/usr/bin/env Rscript
# ==========================================================================
# 26_signature_stability.R
# Bootstrap stability-selection of the 25-gene TCGA LASSO-Cox signature.
# Addresses the reviewer critique that the signature was selected once with
# no stability check. The ENTIRE screening + LASSO pipeline of
# analysis/07_external_validation.R (identical input, pre-filter, glmnet Cox
# LASSO, lambda.min rule) is refit inside B=200 bootstrap resamples of the
# TCGA training data; per-gene selection frequency is recorded.
#
# HARD RULE: real data only. No simulation/imputation. Every number is a fit.
# ==========================================================================
suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})
set.seed(1105)

OUT <- "results/signature_stability"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
logf <- file(file.path(OUT, "run_log.txt"), open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, logf); flush(logf) }

zscore_rows <- function(m) t(scale(t(m)))
B <- 200

# ---------------------------------------------- TCGA training (mirror 07)
load("results/rdata/tcga_processed.RData")             # col_data, tcga_vst
stopifnot(identical(colnames(tcga_vst), rownames(col_data)))
is_tumor <- col_data$status == "Tumor"
cd <- col_data[is_tumor, ]
expr_t <- tcga_vst[, is_tumor, drop = FALSE]
OS_time  <- ifelse(cd$vital_status == "Dead", cd$days_to_death,
                   cd$days_to_last_follow_up)
OS_event <- as.integer(cd$vital_status == "Dead")
keep <- !is.na(OS_time) & OS_time > 0 & !is.na(OS_event)
expr_t <- expr_t[, keep, drop = FALSE]
OS_time <- OS_time[keep]; OS_event <- OS_event[keep]
say(sprintf("TCGA tumors with usable OS: %d (events=%d)",
            length(OS_time), sum(OS_event)))

# same pre-filter as 07: transferable candidate genes (TCGA n ACRG)
load("data/geo/GSE62254.rda")
candidates <- intersect(rownames(expr_t), rownames(GSE62254.expr))
expr_t <- expr_t[candidates, , drop = FALSE]
say(sprintf("Candidate transferable genes (TCGA n GSE62254): %d",
            length(candidates)))

# one pipeline pass -> selected gene names + model size
fit_pipeline <- function(expr, time, event) {
  Z <- zscore_rows(expr)
  Z <- Z[stats::complete.cases(Z), , drop = FALSE]     # drop zero-var genes
  surv <- Surv(time, event)
  uni_p <- apply(Z, 1, function(g) {
    f <- tryCatch(coxph(surv ~ g), error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    summary(f)$coefficients[1, "Pr(>|z|)"]
  })
  sig <- names(uni_p)[!is.na(uni_p) & uni_p < 0.05]
  if (length(sig) < 2) return(character(0))
  X <- t(Z[sig, , drop = FALSE])
  cv <- tryCatch(cv.glmnet(X, surv, family = "cox", alpha = 1, nfolds = 10),
                 error = function(e) NULL)
  if (is.null(cv)) return(character(0))
  co <- as.matrix(coef(cv, s = cv$lambda.min))
  rownames(co)[co[, 1] != 0]
}

# ---------------------------------------------- point estimate (full data)
sel_full <- fit_pipeline(expr_t, OS_time, OS_event)
say(sprintf("\nFull-data LASSO selected %d genes.", length(sel_full)))
sig_ref <- read.csv("results/validation/signature_coefficients.csv")$gene
say(sprintf("Published 25-gene signature; %d/%d reproduced in this full-data pass.",
            length(intersect(sel_full, sig_ref)), length(sig_ref)))

# ---------------------------------------------- B bootstraps
n <- length(OS_time)
sel_count <- setNames(integer(length(candidates)), candidates)
model_size <- integer(B)
for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  sg <- fit_pipeline(expr_t[, idx, drop = FALSE], OS_time[idx], OS_event[idx])
  model_size[b] <- length(sg)
  sel_count[sg] <- sel_count[sg] + 1L
  if (b %% 25 == 0) say(sprintf("  bootstrap %d/%d done (model size=%d)", b, B, length(sg)))
}
freq <- sel_count / B

# ---------------------------------------------- selection-frequency table
freq_df <- data.frame(gene = names(freq),
                      selection_freq = round(as.numeric(freq), 3),
                      in_signature = names(freq) %in% sig_ref,
                      row.names = NULL)
freq_df <- freq_df[order(-freq_df$selection_freq), ]
write.csv(freq_df, file.path(OUT, "selection_frequency.csv"), row.names = FALSE)

# signature-gene focused table (all 25, incl. any not in candidate pool)
coef_ref <- read.csv("results/validation/signature_coefficients.csv")
sig_freq <- data.frame(
  gene = coef_ref$gene, coefficient = round(coef_ref$coefficient, 4),
  in_candidate_pool = coef_ref$gene %in% candidates,
  selection_freq = round(as.numeric(freq[coef_ref$gene]), 3),
  row.names = NULL)
sig_freq$selection_freq[is.na(sig_freq$selection_freq)] <- NA
sig_freq <- sig_freq[order(-replace(sig_freq$selection_freq,
                                    is.na(sig_freq$selection_freq), -1)), ]
write.csv(sig_freq, file.path(OUT, "signature_gene_selection_frequency.csv"),
          row.names = FALSE)

# ---------------------------------------------- summary
sf <- sig_freq$selection_freq
n_over50 <- sum(sf > 0.50, na.rm = TRUE)
n_over80 <- sum(sf > 0.80, na.rm = TRUE)
stromal <- c("SERPINE1","POSTN","MATN3")
summ <- data.frame(
  metric = c("bootstraps_B", "candidate_pool_size",
             "median_model_size", "min_model_size", "max_model_size",
             "signature_genes_total", "signature_in_candidate_pool",
             "signature_selected_gt50pct", "signature_selected_gt80pct",
             "median_selection_freq_of_25", "SERPINE1_freq", "POSTN_freq",
             "MATN3_freq"),
  value = c(B, length(candidates),
            median(model_size), min(model_size), max(model_size),
            nrow(sig_freq), sum(sig_freq$in_candidate_pool),
            n_over50, n_over80,
            round(median(sf, na.rm = TRUE), 3),
            round(freq["SERPINE1"], 3), round(freq["POSTN"], 3),
            round(freq["MATN3"], 3)),
  row.names = NULL)
write.csv(summ, file.path(OUT, "stability_summary.csv"), row.names = FALSE)

say(sprintf("\nModel size across %d bootstraps: median=%d [min=%d, max=%d]",
            B, median(model_size), min(model_size), max(model_size)))
say(sprintf("Of the 25 signature genes: %d in candidate pool; %d selected in >50%%; %d in >80%%.",
            sum(sig_freq$in_candidate_pool), n_over50, n_over80))
say("Stromal anchors: SERPINE1=", round(freq["SERPINE1"],3),
    "  POSTN=", round(freq["POSTN"],3), "  MATN3=", round(freq["MATN3"],3))
say("\nPer-signature-gene selection frequency:")
for (i in seq_len(nrow(sig_freq)))
  say(sprintf("  %-10s coef=%+.3f  freq=%s%s", sig_freq$gene[i],
              sig_freq$coefficient[i],
              ifelse(is.na(sig_freq$selection_freq[i]), "NA(not in pool)",
                     sprintf("%.3f", sig_freq$selection_freq[i])),
              ifelse(sig_freq$gene[i] %in% stromal, "  <-stromal", "")))

save(freq, model_size, sig_freq, summ,
     file = file.path(OUT, "signature_stability.RData"))
say("\nDONE. Outputs in ", OUT)
close(logf)

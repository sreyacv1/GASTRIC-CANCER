#!/usr/bin/env Rscript
# 07_external_validation.R
# Prognostic gene-expression signature: TCGA-STAD training,
# independent external validation on ACRG cohort (GSE62254).
# Cross-platform (RNA-seq -> microarray): genes are z-scored WITHIN each
# cohort before scoring to mitigate the platform difference.
# HARD RULE: real data only. No simulation/imputation of values.

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(survminer)
})
set.seed(42)  # cv.glmnet fold assignment reproducibility only

outdir <- "results/validation"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

zscore_rows <- function(m) {  # z-score each gene (row) across samples
  t(scale(t(m)))
}

## ---------------------------------------------------------------------------
## 1. TCGA training data: tumors, OS, candidate genes
## ---------------------------------------------------------------------------
load("results/rdata/tcga_processed.RData")  # col_data, tcga_vst, res
stopifnot(identical(colnames(tcga_vst), rownames(col_data)))

is_tumor <- col_data$status == "Tumor"
cd <- col_data[is_tumor, ]
expr_t <- tcga_vst[, is_tumor, drop = FALSE]

OS_time  <- ifelse(cd$vital_status == "Dead",
                   cd$days_to_death, cd$days_to_last_follow_up)
OS_event <- as.integer(cd$vital_status == "Dead")
keep <- !is.na(OS_time) & OS_time > 0 & !is.na(OS_event)
cd <- cd[keep, ]; expr_t <- expr_t[, keep, drop = FALSE]
OS_time <- OS_time[keep]; OS_event <- OS_event[keep]
cat(sprintf("TCGA tumors with usable OS: %d (events=%d)\n",
            length(OS_time), sum(OS_event)))

## Validation cohort loaded now to define transferable candidate set
load("data/geo/GSE62254.rda")  # GSE62254.expr, GSE62254.subtype
stopifnot(identical(colnames(GSE62254.expr),
                    as.character(GSE62254.subtype$GEO_ID)))

candidates <- intersect(rownames(expr_t), rownames(GSE62254.expr))
cat(sprintf("Candidate transferable genes (intersection): %d\n",
            length(candidates)))
expr_t <- expr_t[candidates, , drop = FALSE]

## ---------------------------------------------------------------------------
## 2a. Univariate Cox screen on TCGA (z-scored genes), keep p<0.05
## ---------------------------------------------------------------------------
Z_tcga <- zscore_rows(expr_t)                    # genes x samples, z-scored
Z_tcga <- Z_tcga[stats::complete.cases(Z_tcga), , drop = FALSE]  # drop zero-var genes
surv_tcga <- Surv(OS_time, OS_event)

uni_p <- apply(Z_tcga, 1, function(g) {
  fit <- tryCatch(coxph(surv_tcga ~ g), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  summary(fit)$coefficients[1, "Pr(>|z|)"]
})
sig_genes <- names(uni_p)[!is.na(uni_p) & uni_p < 0.05]
cat(sprintf("Univariate-significant genes (p<0.05): %d\n", length(sig_genes)))

## ---------------------------------------------------------------------------
## 2b. LASSO-Cox on the significant genes
## ---------------------------------------------------------------------------
X <- t(Z_tcga[sig_genes, , drop = FALSE])        # samples x genes
y <- surv_tcga
cvfit <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)
lam <- cvfit$lambda.min
coef_lasso <- as.matrix(coef(cvfit, s = lam))
sel <- coef_lasso[coef_lasso[, 1] != 0, , drop = FALSE]
sig_signature <- rownames(sel)
coefs <- setNames(sel[, 1], sig_signature)
cat(sprintf("LASSO (lambda.min=%.4f) selected %d genes:\n", lam, length(coefs)))
print(round(coefs, 4))
stopifnot(length(coefs) >= 2)

sig_df <- data.frame(gene = names(coefs), coefficient = as.numeric(coefs),
                     row.names = NULL)
write.csv(sig_df, file.path(outdir, "signature_coefficients.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. TCGA risk score, median split, KM, C-index
## ---------------------------------------------------------------------------
score_cohort <- function(Zmat, genes, coefs) {
  # Zmat: genes x samples (already z-scored within cohort)
  as.numeric(coefs[genes] %*% Zmat[genes, , drop = FALSE])
}
rs_tcga <- score_cohort(Z_tcga, names(coefs), coefs)
grp_tcga <- factor(ifelse(rs_tcga > median(rs_tcga), "High", "Low"),
                   levels = c("Low", "High"))

df_tcga <- data.frame(time = OS_time, event = OS_event,
                      risk = rs_tcga, group = grp_tcga)
cox_tcga <- coxph(Surv(time, event) ~ risk, data = df_tcga)
cidx_tcga <- summary(cox_tcga)$concordance["C"]
sd_tcga   <- coxph(Surv(time, event) ~ group, data = df_tcga)
lr_tcga   <- survdiff(Surv(time, event) ~ group, data = df_tcga)
p_tcga    <- 1 - pchisq(lr_tcga$chisq, df = length(lr_tcga$n) - 1)
cat(sprintf("TCGA: C-index=%.3f  log-rank p=%.3e\n", cidx_tcga, p_tcga))

km_tcga <- survfit(Surv(time, event) ~ group, data = df_tcga)
p1 <- ggsurvplot(km_tcga, data = df_tcga, pval = TRUE, risk.table = TRUE,
                 palette = c("#377EB8", "#E41A1C"), legend.labs = c("Low", "High"),
                 xlab = "OS time (days)", title = "TCGA-STAD (training)",
                 conf.int = FALSE)
ggsave(file.path(outdir, "KM_TCGA.png"), plot = p1$plot,
       width = 6, height = 5, dpi = 150)

## ---------------------------------------------------------------------------
## 4. EXTERNAL VALIDATION on ACRG / GSE62254
## ---------------------------------------------------------------------------
st <- GSE62254.subtype
acrg_time  <- st$OS.m
acrg_event <- st$Death
keepA <- !is.na(acrg_time) & acrg_time > 0 & !is.na(acrg_event)
Z_acrg <- zscore_rows(GSE62254.expr)[names(coefs), keepA, drop = FALSE]
st_a <- st[keepA, ]
acrg_time <- acrg_time[keepA]; acrg_event <- acrg_event[keepA]
cat(sprintf("ACRG patients with usable OS: %d (events=%d)\n",
            length(acrg_time), sum(acrg_event)))
stopifnot(all(!is.na(Z_acrg)))  # all signature genes measured, non-constant

rs_acrg <- score_cohort(Z_acrg, names(coefs), coefs)
grp_acrg <- factor(ifelse(rs_acrg > median(rs_acrg), "High", "Low"),
                   levels = c("Low", "High"))
df_acrg <- data.frame(time = acrg_time, event = acrg_event,
                      risk = rs_acrg, group = grp_acrg,
                      Stage = factor(st_a$Stage, levels = c("I","II","III","IV")),
                      age = st_a$age)

cox_acrg <- coxph(Surv(time, event) ~ risk, data = df_acrg)
cidx_acrg <- summary(cox_acrg)$concordance["C"]
lr_acrg <- survdiff(Surv(time, event) ~ group, data = df_acrg)
p_acrg  <- 1 - pchisq(lr_acrg$chisq, df = length(lr_acrg$n) - 1)
cat(sprintf("ACRG: C-index=%.3f  log-rank p=%.3e\n", cidx_acrg, p_acrg))

km_acrg <- survfit(Surv(time, event) ~ group, data = df_acrg)
p2 <- ggsurvplot(km_acrg, data = df_acrg, pval = TRUE, risk.table = TRUE,
                 palette = c("#377EB8", "#E41A1C"), legend.labs = c("Low", "High"),
                 xlab = "OS time (months)", title = "ACRG / GSE62254 (validation)",
                 conf.int = FALSE)
ggsave(file.path(outdir, "KM_ACRG.png"), plot = p2$plot,
       width = 6, height = 5, dpi = 150)

## ---------------------------------------------------------------------------
## 5. Multivariable Cox in ACRG: risk group + Stage + age
## ---------------------------------------------------------------------------
mv <- coxph(Surv(time, event) ~ group + Stage + age, data = df_acrg)
smv <- summary(mv)
mv_tab <- data.frame(
  term    = rownames(smv$coefficients),
  HR      = smv$conf.int[, "exp(coef)"],
  CI_low  = smv$conf.int[, "lower .95"],
  CI_high = smv$conf.int[, "upper .95"],
  p_value = smv$coefficients[, "Pr(>|z|)"],
  row.names = NULL)
write.csv(mv_tab, file.path(outdir, "multivariable_cox_ACRG.csv"),
          row.names = FALSE)
cat("Multivariable Cox (ACRG):\n"); print(mv_tab)
rg <- mv_tab[mv_tab$term == "groupHigh", ]
cat(sprintf("ACRG risk-group HR=%.2f (95%% CI %.2f-%.2f), p=%.3e\n",
            rg$HR, rg$CI_low, rg$CI_high, rg$p_value))

## ---------------------------------------------------------------------------
## C-index comparison table
## ---------------------------------------------------------------------------
cidx_tab <- data.frame(
  cohort  = c("TCGA-STAD (training)", "ACRG/GSE62254 (validation)"),
  n       = c(nrow(df_tcga), nrow(df_acrg)),
  events  = c(sum(df_tcga$event), sum(df_acrg$event)),
  C_index = c(cidx_tcga, cidx_acrg),
  logrank_p = c(p_tcga, p_acrg))
write.csv(cidx_tab, file.path(outdir, "cindex_comparison.csv"), row.names = FALSE)
cat("\nC-index comparison:\n"); print(cidx_tab)

cat("\n=== DONE. Outputs in results/validation/ ===\n")

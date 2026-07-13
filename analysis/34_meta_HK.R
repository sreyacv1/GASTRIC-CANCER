#!/usr/bin/env Rscript
# 34_meta_HK.R
# UPGRADE 3: Hartung-Knapp REML meta-analysis of the CONTINUOUS per-1-SD
# signature effect across the 3 external cohorts, replacing the fragile
# DerSimonian-Laird pooled HR (1.50) built on median-dichotomised effects.
#
# Common estimand: per-1-SD signature log(HR), age + stage adjusted, scored
# with the fixed 25-gene TCGA signature (z-scored within each cohort exactly
# as analysis/07 & 12). Cohorts: ACRG/GSE62254, GSE15459, GSE84437.
# Pooling: metafor::rma(method="REML", test="knha") -> pooled HR, 95% CI, p,
# I2, tau2, and the 95% PREDICTION interval.
#
# HARD RULE: real data only, complete-case, no simulation/imputation.

suppressPackageStartupMessages({
  library(survival); library(Biobase); library(readxl); library(metafor)
})
set.seed(1105)

outdir <- "results/meta_HK"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sig   <- read.csv("results/validation/signature_coefficients.csv",
                  stringsAsFactors = FALSE)
coefs <- setNames(sig$coefficient, sig$gene)

zscore_rows <- function(m) t(scale(t(m)))
collapse_by_symbol <- function(mat, symbols) {   # gene x sample, max-mean probe
  symbols <- sub("///.*$", "", trimws(symbols))
  keep <- !is.na(symbols) & symbols != "" & !is.na(rowMeans(mat, na.rm = TRUE))
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; symbols <- symbols[ord]
  mat <- mat[!duplicated(symbols), , drop = FALSE]
  rownames(mat) <- symbols[!duplicated(symbols)]
  mat
}

# per-1-SD adjusted log(HR): risk (unit SD) + covars, complete-case
score_effect <- function(name, expr, time, event, covars, adj_label) {
  ok <- !is.na(time) & time > 0 & !is.na(event)
  expr <- expr[, ok, drop = FALSE]; time <- time[ok]; event <- event[ok]
  covars <- covars[ok, , drop = FALSE]
  present <- intersect(names(coefs), rownames(expr))
  Z <- zscore_rows(expr)[present, , drop = FALSE]
  Z <- Z[stats::complete.cases(Z), , drop = FALSE]
  present <- rownames(Z)
  risk <- as.numeric(coefs[present] %*% Z)
  df <- data.frame(time = time, event = event,
                   risk = as.numeric(scale(risk)), covars)
  cc <- complete.cases(df)
  df <- df[cc, ]
  form <- as.formula(paste("Surv(time, event) ~ risk +",
                           paste(colnames(covars), collapse = " + ")))
  fit <- coxph(form, data = df)
  s <- summary(fit)$coefficients["risk", ]
  cat(sprintf("[%s] n=%d events=%d genes=%d/%d  logHR/SD=%.3f SE=%.3f\n",
              name, nrow(df), sum(df$event), length(present), length(coefs),
              s["coef"], s["se(coef)"]))
  data.frame(cohort = name, n = nrow(df), events = sum(df$event),
             genes_mapped = length(present), adjustment = adj_label,
             logHR = s["coef"], SE = s["se(coef)"],
             HR = exp(s["coef"]),
             HR_low = exp(s["coef"] - 1.96 * s["se(coef)"]),
             HR_high = exp(s["coef"] + 1.96 * s["se(coef)"]),
             p = s["Pr(>|z|)"], row.names = NULL)
}

rows <- list()

## ---- ACRG / GSE62254 (age + Stage I-IV) ----------------------------------
load("data/geo/GSE62254.rda")
st <- GSE62254.subtype
rows[["ACRG"]] <- score_effect(
  "ACRG/GSE62254", GSE62254.expr[names(coefs), , drop = FALSE],
  as.numeric(st$OS.m), as.integer(st$Death),
  data.frame(Age = as.numeric(st$age),
             Stage = factor(st$Stage, levels = c("I", "II", "III", "IV"))),
  "age + TNM stage")

## ---- GSE15459 (age + Stage 1-4 from outcome xls) -------------------------
es2 <- readRDS("data/geo/GSE15459_es.rds")
expr2 <- collapse_by_symbol(exprs(es2), fData(es2)$`Gene symbol`)
oc <- as.data.frame(read_excel("data/geo/GSE15459_outcome.xls"))
idx <- match(colnames(expr2), oc$`GSM ID`)
expr2 <- expr2[, !is.na(idx), drop = FALSE]; oc <- oc[idx[!is.na(idx)], ]
rows[["GSE15459"]] <- score_effect(
  "GSE15459", expr2,
  as.numeric(oc$`Overall.Survival (Months)**`), as.numeric(oc$`Outcome (1=dead)`),
  data.frame(Age = as.numeric(oc$Age_at_surgery),
             Stage = factor(oc$Stage, levels = c("1", "2", "3", "4"))),
  "age + TNM stage")

## ---- GSE84437 (age + T-stage; no TNM overall stage available) ------------
es1 <- readRDS("data/geo/GSE84437_es.rds")
pd1 <- pData(es1)
expr1 <- collapse_by_symbol(exprs(es1), fData(es1)$`Gene symbol`)
rows[["GSE84437"]] <- score_effect(
  "GSE84437", expr1,
  as.numeric(pd1$`duration overall survival:ch1`), as.numeric(pd1$`death:ch1`),
  data.frame(Age = as.numeric(pd1$`age:ch1`),
             Tstage = factor(pd1$`ptstage:ch1`, levels = c("T1", "T2", "T3", "T4"))),
  "age + T-stage (no TNM stage in GEO; T-stage proxy)")

meta_in <- do.call(rbind, rows)
row.names(meta_in) <- NULL
write.csv(meta_in, file.path(outdir, "meta_inputs.csv"), row.names = FALSE)
cat("\n== Per-cohort per-1-SD adjusted effects ==\n"); print(meta_in)

## ---------------------------------------------------------------------------
## Hartung-Knapp REML meta-analysis + prediction interval
## ---------------------------------------------------------------------------
res <- rma(yi = logHR, sei = SE, data = meta_in, method = "REML", test = "knha")
pr  <- predict(res)                               # t-based CI + PI (knha)
pooled_HR  <- exp(as.numeric(coef(res)))
ci_lo <- exp(res$ci.lb); ci_hi <- exp(res$ci.ub)
pi_lo <- exp(pr$pi.lb);  pi_hi <- exp(pr$pi.ub)
I2 <- res$I2; tau2 <- res$tau2; pval <- res$pval

meta_res <- data.frame(
  method = "REML + Hartung-Knapp (metafor rma, test=knha)",
  k = res$k, pooled_HR = pooled_HR, CI_low = ci_lo, CI_high = ci_hi,
  p_value = pval, I2 = I2, tau2 = tau2,
  PI_low = pi_lo, PI_high = pi_hi,
  significant = ifelse(ci_lo > 1 | ci_hi < 1, "YES", "NO (CI crosses 1)"))
write.csv(meta_res, file.path(outdir, "meta_result.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## SUMMARY
## ---------------------------------------------------------------------------
sink(file.path(outdir, "SUMMARY.txt"))
cat("HARTUNG-KNAPP REML meta-analysis of the CONTINUOUS per-1-SD signature\n")
cat("effect (age + stage adjusted), replacing the DerSimonian-Laird pooled\n")
cat("HR (1.50) built on median-dichotomised effects.\n\n")
cat("== Per-cohort inputs (per-1-SD adjusted HR) ==\n")
print(meta_in[, c("cohort", "n", "events", "genes_mapped", "adjustment",
                  "HR", "HR_low", "HR_high", "p")], row.names = FALSE)
cat("\n== Pooled (REML + Hartung-Knapp) ==\n")
cat(sprintf("Pooled HR = %.2f (95%% CI %.2f-%.2f), p = %.3f\n",
            pooled_HR, ci_lo, ci_hi, pval))
cat(sprintf("I2 = %.1f%%   tau2 = %.4f   (k = %d studies)\n", I2, tau2, res$k))
cat(sprintf("95%% prediction interval: %.2f - %.2f\n", pi_lo, pi_hi))
cat(sprintf("\nPooled effect significant at alpha=0.05: %s\n",
            ifelse(ci_lo > 1 | ci_hi < 1, "YES",
                   "NO -- the CI CROSSES the null (HR=1)")))
cat("\nINTERPRETATION: with only k=3 heterogeneous cohorts, the Hartung-Knapp\n")
cat("small-sample correction widens the CI relative to DerSimonian-Laird.\n")
cat("Report the pooled HR WITH its prediction interval and state the\n")
cat("significance honestly.\n")
sink()

cat("\n=== DONE. Outputs in", outdir, "===\n")
cat(sprintf("HK pooled HR=%.2f (%.2f-%.2f) p=%.3f; PI %.2f-%.2f; I2=%.1f%%\n",
            pooled_HR, ci_lo, ci_hi, pval, pi_lo, pi_hi, I2))

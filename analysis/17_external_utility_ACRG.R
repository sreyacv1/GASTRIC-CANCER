#!/usr/bin/env Rscript
# 17_external_utility_ACRG.R
# EXTERNAL clinical-utility validation of the TCGA-derived prognostic
# signature, computed entirely in the independent ACRG cohort (GSE62254).
#
# Reviewer must-fix: DCA / IDI / NRI / time-AUC for the combined
# clinical+signature model were previously computed only in the TCGA
# training cohort (circular, since the 25-gene signature was derived there).
# Here every utility statistic is computed OUT-OF-SAMPLE in ACRG (n=300).
#
# The 25-gene signature (genes + LASSO-Cox coefficients) is fixed from TCGA;
# it is applied to ACRG (z-scored within ACRG). The three Cox models are then
# fit in ACRG (Clinical = Stage+age; Signature = risk; Combined = risk+Stage
# +age). ACRG has NO tumour grade, so Grade is excluded (noted).
#
# HARD RULE: real data only, complete-case, no simulation/imputation.

suppressPackageStartupMessages({
  library(survival)
  library(survIDINRI)
  library(survivalROC)
  library(dcurves)
})
set.seed(42)

outdir <- "results/external_utility_ACRG"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

zscore_rows <- function(m) t(scale(t(m)))          # z-score each gene (row)
# Harrell C for a fixed linear predictor (higher lp = higher risk)
cidx <- function(time, event, lp)
  survival::concordance(Surv(time, event) ~ lp, reverse = TRUE)$concordance

## ---------------------------------------------------------------------------
## 1. ACRG cohort: OS (years), risk score from fixed TCGA signature
## ---------------------------------------------------------------------------
load("data/geo/GSE62254.rda")                      # GSE62254.expr, .subtype
stopifnot(identical(colnames(GSE62254.expr),
                    as.character(GSE62254.subtype$GEO_ID)))
sig   <- read.csv("results/validation/signature_coefficients.csv",
                  stringsAsFactors = FALSE)
coefs <- setNames(sig$coefficient, sig$gene)
stopifnot(all(names(coefs) %in% rownames(GSE62254.expr)))

st  <- GSE62254.subtype
ZA  <- zscore_rows(GSE62254.expr[names(coefs), , drop = FALSE])  # z within ACRG
stopifnot(all(!is.na(ZA)))                         # all sig genes non-constant
risk_acrg <- as.numeric(coefs[rownames(ZA)] %*% ZA)

acrg <- data.frame(
  time  = as.numeric(st$OS.m) / 12,                # months -> years
  event = as.integer(st$Death),
  Age   = as.numeric(st$age),
  Stage = factor(st$Stage, levels = c("I", "II", "III", "IV")),
  risk  = as.numeric(scale(risk_acrg)))            # z-scored risk within ACRG
cc <- complete.cases(acrg) & acrg$time > 0
acrg <- acrg[cc, ]
cat(sprintf("ACRG complete cases: n=%d, events=%d\n",
            nrow(acrg), sum(acrg$event)))

## ---------------------------------------------------------------------------
## 2. Three Cox models (fit in ACRG) + optimism-corrected Harrell C-index
##    Efron-Gong bootstrap optimism (B=300).
## ---------------------------------------------------------------------------
f_clin <- Surv(time, event) ~ Age + Stage
f_sig  <- Surv(time, event) ~ risk
f_comb <- Surv(time, event) ~ Age + Stage + risk

m_clin <- coxph(f_clin, data = acrg, x = TRUE, y = TRUE)
m_sig  <- coxph(f_sig,  data = acrg, x = TRUE, y = TRUE)
m_comb <- coxph(f_comb, data = acrg, x = TRUE, y = TRUE)

boot_c <- function(formula, data, B = 300) {
  full     <- coxph(formula, data = data)
  apparent <- cidx(data$time, data$event, predict(full, type = "lp"))
  se       <- as.numeric(summary(full)$concordance[2])  # SE of apparent C
  opt <- numeric(0)
  for (b in seq_len(B)) {
    idx <- sample(nrow(data), replace = TRUE)
    bd  <- data[idx, ]
    fit <- tryCatch(coxph(formula, data = bd), error = function(e) NULL)
    if (is.null(fit) || any(is.na(coef(fit)))) next
    cb <- tryCatch(cidx(bd$time,   bd$event,   predict(fit, newdata = bd,   type = "lp")),
                   error = function(e) NA)
    co <- tryCatch(cidx(data$time, data$event, predict(fit, newdata = data, type = "lp")),
                   error = function(e) NA)
    if (is.na(cb) || is.na(co)) next
    opt <- c(opt, cb - co)                          # optimism_b
  }
  optimism  <- mean(opt)
  corrected <- apparent - optimism
  c(apparent = apparent, optimism = optimism, corrected = corrected,
    lo = corrected - 1.96 * se, hi = corrected + 1.96 * se, nboot = length(opt))
}

set.seed(42); bc_clin <- boot_c(f_clin, acrg)
set.seed(42); bc_sig  <- boot_c(f_sig,  acrg)
set.seed(42); bc_comb <- boot_c(f_comb, acrg)

cindex_tab <- data.frame(
  model  = c("Clinical (Stage+age)", "Signature (risk)",
             "Combined (risk+Stage+age)"),
  cohort = "ACRG/GSE62254 (external)",
  n = nrow(acrg), events = sum(acrg$event),
  C_apparent  = c(bc_clin["apparent"],  bc_sig["apparent"],  bc_comb["apparent"]),
  C_corrected = c(bc_clin["corrected"], bc_sig["corrected"], bc_comb["corrected"]),
  C_lo = c(bc_clin["lo"], bc_sig["lo"], bc_comb["lo"]),
  C_hi = c(bc_clin["hi"], bc_sig["hi"], bc_comb["hi"]),
  B = 300, row.names = NULL)
cat("\n== External C-index (optimism-corrected, B=300) ==\n"); print(cindex_tab)
write.csv(cindex_tab, file.path(outdir, "cindex_external.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. External added value: Combined vs Clinical
##    dC-index, nested LRT, IDI + continuous NRI @3y (survIDINRI::IDI.INF)
## ---------------------------------------------------------------------------
dC  <- as.numeric(bc_comb["corrected"] - bc_clin["corrected"])
lrt <- anova(m_clin, m_comb, test = "LRT")          # nested Cox LRT
lrt_chisq <- lrt$Chisq[2]; lrt_df <- lrt$Df[2]; lrt_p <- lrt$`Pr(>|Chi|)`[2]
cat(sprintf("\ndeltaC (combined-clinical, corrected) = %.4f\n", dC))
cat(sprintf("Nested LRT chisq=%.2f df=%d p=%.3e\n", lrt_chisq, lrt_df, lrt_p))

t3 <- 3                                              # 3 years
X0 <- model.matrix(~ Age + Stage,        data = acrg)[, -1, drop = FALSE]
X1 <- model.matrix(~ Age + Stage + risk, data = acrg)[, -1, drop = FALSE]
outc <- as.matrix(acrg[, c("time", "event")])
set.seed(42)
ii <- IDI.INF(outc, X0, X1, t0 = t3, npert = 300)
o  <- as.matrix(IDI.INF.OUT(ii))   # rows: M1=IDI, M2=continuous NRI, M3=median
idi_nri <- data.frame(
  metric   = c("IDI@3y", "continuousNRI@3y", "medianDiff@3y"),
  estimate = o[, 1], lo = o[, 2], hi = o[, 3], p_value = o[, 4],
  method   = "survIDINRI::IDI.INF (npert=300)", row.names = NULL)
cat("\n== External added value (IDI / NRI @3y) ==\n"); print(idi_nri)

added <- data.frame(
  statistic = c("deltaC_corrected(combined-clinical)", "LRT_chisq", "LRT_df",
                "LRT_p", "IDI@3y", "IDI_p", "continuousNRI@3y", "NRI_p"),
  value = c(dC, lrt_chisq, lrt_df, lrt_p,
            idi_nri$estimate[1], idi_nri$p_value[1],
            idi_nri$estimate[2], idi_nri$p_value[2]))
write.csv(added, file.path(outdir, "added_value_external.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 4. External time-dependent AUC @ 1/3/5 years (survivalROC, KM method)
## ---------------------------------------------------------------------------
tvec <- c(1, 3, 5)
auc_of <- function(lp, lab) {
  aucs <- sapply(tvec, function(tt)
    survivalROC(Stime = acrg$time, status = acrg$event, marker = lp,
                predict.time = tt, method = "KM")$AUC)
  data.frame(model = lab, t_years = tvec, AUC = aucs)
}
auc_tab <- rbind(
  auc_of(predict(m_clin, type = "lp"), "Clinical"),
  auc_of(predict(m_sig,  type = "lp"), "Signature"),
  auc_of(predict(m_comb, type = "lp"), "Combined"))
cat("\n== External time-dependent AUC ==\n"); print(auc_tab)
write.csv(auc_tab, file.path(outdir, "timeAUC_external.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5. External Decision Curve Analysis @ 3 years (dcurves)
## ---------------------------------------------------------------------------
risk_at_t <- function(fit, t)
  1 - as.numeric(summary(survfit(fit, newdata = acrg), times = t)$surv)
acrg$p_clin <- risk_at_t(m_clin, t3)
acrg$p_comb <- risk_at_t(m_comb, t3)

d <- dca(Surv(time, event) ~ p_clin + p_comb, data = acrg,
         time = t3, thresholds = seq(0, 0.5, 0.01),
         label = list(p_clin = "Clinical", p_comb = "Combined"))
dsum <- as.data.frame(d$dca)
write.csv(dsum, file.path(outdir, "DCA_netbenefit_external.csv"), row.names = FALSE)
png(file.path(outdir, "DCA_external.png"), width = 1600, height = 1200, res = 200)
print(plot(d, smooth = TRUE) +
        ggplot2::labs(title = "External DCA @3y (ACRG/GSE62254)"))
dev.off()

# Does Combined dominate Clinical out-of-sample? Compare net benefit over a
# clinically relevant threshold range (5-40%).
nb <- dsum[dsum$threshold >= 0.05 & dsum$threshold <= 0.40, ]
nb_w <- reshape(nb[, c("threshold", "label", "net_benefit")],
                idvar = "threshold", timevar = "label", direction = "wide")
colnames(nb_w) <- sub("net_benefit\\.", "", colnames(nb_w))
combined_col <- nb_w[["Combined"]]; clinical_col <- nb_w[["Clinical"]]
frac_comb_ge_clin <- mean(combined_col >= clinical_col - 1e-6)
mean_gain <- mean(combined_col - clinical_col)
cat(sprintf(paste0("\nDCA (thr 5-40%%): Combined >= Clinical net benefit at ",
                   "%.0f%% of thresholds; mean NB gain=%.4f\n"),
            100 * frac_comb_ge_clin, mean_gain))
dca_conclusion <- if (frac_comb_ge_clin >= 0.9)
  "Combined dominates Clinical out-of-sample across the clinical threshold range" else if (mean_gain > 0)
  "Combined has higher mean net benefit but does not strictly dominate Clinical" else
  "Combined does NOT improve net benefit over Clinical out-of-sample"

## ---------------------------------------------------------------------------
## SUMMARY
## ---------------------------------------------------------------------------
sink(file.path(outdir, "SUMMARY.txt"))
cat("EXTERNAL clinical-utility validation in ACRG/GSE62254 (n=",
    nrow(acrg), ", events=", sum(acrg$event), ")\n\n", sep = "")
cat("Signature: 25 genes + coefs FIXED from TCGA; risk z-scored within ACRG.\n")
cat("Cox models fit in ACRG; ACRG has no tumour grade (Grade excluded).\n\n")
cat("== C-index (optimism-corrected Harrell, B=300) ==\n"); print(cindex_tab)
cat("\n== Added value: Combined vs Clinical ==\n"); print(added)
cat("\n== Time-dependent AUC (1/3/5y) ==\n"); print(auc_tab)
cat(sprintf("\n== DCA @3y conclusion ==\n%s\n", dca_conclusion))
cat(sprintf("(Combined >= Clinical at %.0f%% of thresholds 5-40%%, mean gain=%.4f)\n",
            100 * frac_comb_ge_clin, mean_gain))
sink()

cat("\nAll external-utility outputs written to", outdir, "\n")
cat("DCA conclusion:", dca_conclusion, "\n")

#!/usr/bin/env Rscript
# 13_combined_nomogram_DCA.R
# Combined clinical + genomic prognostic model + clinical-utility analysis.
# TCGA-STAD tumor complete cases. External validation on ACRG (GSE62254).
# HARD RULE: real data only. Complete-case analysis. No simulation/imputation.
#
# Models for OS:
#   (a) Clinical  = Age + Stage + Grade
#   (b) Signature = risk score alone
#   (c) Combined  = Age + Stage + Grade + risk score
# Deliverables: optimism-corrected C-index (Harrell bootstrap B=300),
#   added value (dC, LRT, NRI/IDI @3y), tdAUC @1/3/5y, DCA @3y,
#   combined nomogram + calibration, ACRG combined-model C-index.

suppressPackageStartupMessages({
  library(survival)
  library(rms)
})
set.seed(42)

outdir <- "results/nomogram_combined"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
USE_timeROC    <- has_pkg("timeROC")
USE_survIDINRI <- has_pkg("survIDINRI")
USE_dcurves    <- has_pkg("dcurves")
USE_survROC    <- has_pkg("survivalROC")
cat(sprintf("Optional pkgs: timeROC=%s survIDINRI=%s dcurves=%s survivalROC=%s\n",
            USE_timeROC, USE_survIDINRI, USE_dcurves, USE_survROC))

zscore_rows <- function(m) t(scale(t(m)))            # z-score each gene (row)
collapse_stage <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  out[grepl("Stage IV", x)]              <- "IV"
  out[grepl("Stage III", x)]             <- "III"
  out[grepl("Stage II($| |A|B|C)", x)]   <- "II"
  out[grepl("Stage I($| |A|B|C)", x)]    <- "I"
  factor(out, levels = c("I", "II", "III", "IV"))
}
# Harrell C for a fixed risk score (higher lp = higher risk = shorter time)
cidx <- function(time, event, lp) {
  survival::concordance(Surv(time, event) ~ lp, reverse = TRUE)$concordance
}

## ---------------------------------------------------------------------------
## 1. TCGA tumor complete-case cohort + risk score
## ---------------------------------------------------------------------------
load("results/rdata/tcga_processed.RData")          # col_data, tcga_vst, res
sig <- read.csv("results/validation/signature_coefficients.csv",
                stringsAsFactors = FALSE)
coefs <- setNames(sig$coefficient, sig$gene)
stopifnot(all(names(coefs) %in% rownames(tcga_vst)))

is_tum <- col_data$status == "Tumor"
cd  <- col_data[is_tum, ]
ex  <- tcga_vst[names(coefs), is_tum, drop = FALSE]  # signature genes x tumors

OS_time  <- ifelse(cd$vital_status == "Dead",
                   cd$days_to_death, cd$days_to_last_follow_up)
OS_event <- as.integer(cd$vital_status == "Dead")
Age      <- cd$age_at_diagnosis / 365.25
Stage    <- collapse_stage(cd$ajcc_pathologic_stage)
Grade    <- factor(ifelse(cd$tumor_grade %in% c("G1","G2","G3"),
                          cd$tumor_grade, NA), levels = c("G1","G2","G3"))

# Risk score: z-score each gene within this tumor cohort, then coef %*% z
Z    <- zscore_rows(ex)
risk_full <- as.numeric(coefs[rownames(Z)] %*% Z)

df <- data.frame(time = OS_time, event = OS_event, Age = Age,
                 Stage = Stage, Grade = Grade, risk = risk_full)
cc <- complete.cases(df) & df$time > 0
df <- df[cc, ]
cat(sprintf("TCGA tumor complete cases: n=%d, events=%d\n",
            nrow(df), sum(df$event)))

## ---------------------------------------------------------------------------
## 2. Three Cox models + optimism-corrected Harrell C-index (bootstrap B=300)
## ---------------------------------------------------------------------------
f_clin <- Surv(time, event) ~ Age + Stage + Grade
f_sig  <- Surv(time, event) ~ risk
f_comb <- Surv(time, event) ~ Age + Stage + Grade + risk

m_clin <- coxph(f_clin, data = df, x = TRUE, y = TRUE)
m_sig  <- coxph(f_sig,  data = df, x = TRUE, y = TRUE)
m_comb <- coxph(f_comb, data = df, x = TRUE, y = TRUE)

boot_c <- function(formula, data, B = 300) {
  lp0     <- predict(coxph(formula, data = data), type = "lp")
  apparent <- cidx(data$time, data$event, lp0)
  opt <- numeric(0)
  for (b in seq_len(B)) {
    idx <- sample(nrow(data), replace = TRUE)
    bd  <- data[idx, ]
    fit <- tryCatch(coxph(formula, data = bd), error = function(e) NULL)
    if (is.null(fit) || any(is.na(coef(fit)))) next
    cb <- tryCatch(cidx(bd$time, bd$event,
                        predict(fit, newdata = bd, type = "lp")),
                   error = function(e) NA)
    co <- tryCatch(cidx(data$time, data$event,
                        predict(fit, newdata = data, type = "lp")),
                   error = function(e) NA)
    if (is.na(cb) || is.na(co)) next
    opt <- c(opt, cb - co)
  }
  optimism  <- mean(opt)
  corrected <- apparent - optimism
  # bootstrap CI of corrected via SE of the per-model concordance
  se <- as.numeric(summary(coxph(formula, data = data))$concordance[2])
  c(apparent = apparent, optimism = optimism, corrected = corrected,
    lo = corrected - 1.96 * se, hi = corrected + 1.96 * se, nboot = length(opt))
}

set.seed(42); bc_clin <- boot_c(f_clin, df)
set.seed(42); bc_sig  <- boot_c(f_sig,  df)
set.seed(42); bc_comb <- boot_c(f_comb, df)

cindex_tab <- data.frame(
  model     = c("Clinical (Age+Stage+Grade)", "Signature (risk)",
                "Combined (Clinical+risk)"),
  cohort    = "TCGA-STAD",
  n         = nrow(df), events = sum(df$event),
  C_apparent  = c(bc_clin["apparent"],  bc_sig["apparent"],  bc_comb["apparent"]),
  C_corrected = c(bc_clin["corrected"], bc_sig["corrected"], bc_comb["corrected"]),
  C_lo        = c(bc_clin["lo"], bc_sig["lo"], bc_comb["lo"]),
  C_hi        = c(bc_clin["hi"], bc_sig["hi"], bc_comb["hi"]),
  row.names = NULL)
cat("\n== C-index (optimism-corrected, B=300) ==\n"); print(cindex_tab)

## ---------------------------------------------------------------------------
## 3. Added value of signature over clinical
## ---------------------------------------------------------------------------
dC <- bc_comb["corrected"] - bc_clin["corrected"]
lrt <- anova(m_clin, m_comb, test = "LRT")   # nested Cox LRT
lrt_p <- lrt$`Pr(>|Chi|)`[2]
lrt_chisq <- lrt$Chisq[2]
cat(sprintf("\ndeltaC (combined-clinical, corrected) = %.4f\n", dC))
cat(sprintf("LRT chisq=%.2f df=%d p=%.3e\n", lrt_chisq, lrt$Df[2], lrt_p))

t3 <- 3 * 365.25
nri_idi <- data.frame(metric = character(), estimate = numeric(),
                      lo = numeric(), hi = numeric(), method = character())
if (USE_survIDINRI) {
  library(survIDINRI)
  # design matrices (drop intercept); baseline = clinical, new = +risk
  X0 <- model.matrix(~ Age + Stage + Grade, data = df)[, -1, drop = FALSE]
  X1 <- model.matrix(~ Age + Stage + Grade + risk, data = df)[, -1, drop = FALSE]
  outc <- as.matrix(df[, c("time", "event")])
  set.seed(42)
  ii <- IDI.INF(outc, X0, X1, t0 = t3, npert = 300)
  o  <- IDI.INF.OUT(ii)   # rows: M1=IDI, M2=continuous NRI, M3=median improvement
  o  <- as.matrix(o)      # cols: Est., Lower, Upper, p-value (index by position)
  nri_idi <- data.frame(
    metric   = c("IDI@3y", "continuousNRI@3y", "medianDiff@3y"),
    estimate = o[, 1], lo = o[, 2], hi = o[, 3],
    method   = "survIDINRI::IDI.INF", row.names = NULL)
} else {
  # manual category-free NRI @ t3 from predicted event risk at t3
  risk_at <- function(fit) {
    sf <- survfit(fit, newdata = df)
    st <- summary(sf, times = t3)$surv
    1 - as.numeric(st)
  }
  p0 <- risk_at(m_clin); p1 <- risk_at(m_comb)
  # among events by t3 vs event-free by t3 (complete-case: exclude censored<t3)
  ev  <- df$event == 1 & df$time <= t3
  nev <- df$time > t3
  use <- ev | nev
  d   <- (p1 - p0)
  nri_ev  <- mean(sign(d[ev])  > 0) - mean(sign(d[ev])  < 0)
  nri_nev <- mean(sign(d[nev]) < 0) - mean(sign(d[nev]) > 0)
  cNRI <- nri_ev + nri_nev
  idi_manual <- mean(p1[ev]) - mean(p0[ev]) - (mean(p1[nev]) - mean(p0[nev]))
  nri_idi <- data.frame(
    metric   = c("continuousNRI@3y", "IDI@3y"),
    estimate = c(cNRI, idi_manual), lo = NA, hi = NA,
    method   = "manual (survIDINRI unavailable)", row.names = NULL)
}
cat("\n== Added value (NRI/IDI) ==\n"); print(nri_idi)

## ---------------------------------------------------------------------------
## 4. Time-dependent AUC @ 1/3/5 years
## ---------------------------------------------------------------------------
tvec <- c(1, 3, 5) * 365.25
lp_clin <- predict(m_clin, type = "lp")
lp_sig  <- predict(m_sig,  type = "lp")
lp_comb <- predict(m_comb, type = "lp")
auc_tab <- data.frame()
if (USE_timeROC) {
  library(timeROC)
  roc_of <- function(lp, lab) {
    r <- timeROC(T = df$time, delta = df$event, marker = lp,
                 cause = 1, times = tvec, iid = FALSE)
    data.frame(model = lab, t_years = c(1,3,5), AUC = as.numeric(r$AUC))
  }
  auc_tab <- rbind(roc_of(lp_clin, "Clinical"),
                   roc_of(lp_sig,  "Signature"),
                   roc_of(lp_comb, "Combined"))
} else if (USE_survROC) {
  library(survivalROC)
  roc_of <- function(lp, lab) {
    aucs <- sapply(tvec, function(tt)
      survivalROC(Stime = df$time, status = df$event, marker = lp,
                  predict.time = tt, method = "KM")$AUC)
    data.frame(model = lab, t_years = c(1,3,5), AUC = aucs)
  }
  auc_tab <- rbind(roc_of(lp_clin, "Clinical"),
                   roc_of(lp_sig,  "Signature"),
                   roc_of(lp_comb, "Combined"))
}
cat("\n== Time-dependent AUC ==\n"); print(auc_tab)

## ---------------------------------------------------------------------------
## 5. Decision Curve Analysis @ 3 years
## ---------------------------------------------------------------------------
risk_at_t <- function(fit, t) {
  sf <- survfit(fit, newdata = df)
  1 - as.numeric(summary(sf, times = t)$surv)
}
df$p_clin <- risk_at_t(m_clin, t3)
df$p_comb <- risk_at_t(m_comb, t3)
dca_done <- FALSE
if (USE_dcurves) {
  library(dcurves)
  d <- dcurves::dca(Surv(time, event) ~ p_clin + p_comb, data = df,
                    time = t3, thresholds = seq(0, 0.5, 0.01))
  saveRDS(d, file.path(outdir, "dca_object.rds"))
  dsum <- as.data.frame(d$dca)
  write.csv(dsum, file.path(outdir, "dca_netbenefit.csv"), row.names = FALSE)
  png(file.path(outdir, "DCA_3yr.png"), width = 1600, height = 1200, res = 200)
  print(plot(d, smooth = TRUE))
  dev.off()
  dca_done <- TRUE
} else {
  # manual net benefit: NB = TP/n - FP/n * (pt/(1-pt)); event by t3 (KM-free
  # complete-case: use observed status, censored-before-t3 excluded from n)
  thr <- seq(0.01, 0.5, 0.01)
  ev  <- df$event == 1 & df$time <= t3
  nev <- df$time > t3
  use <- ev | nev
  n <- sum(use); prev <- mean(ev[use])
  nb <- function(p) sapply(thr, function(pt) {
    pos <- p[use] >= pt
    tp <- mean(pos & ev[use]); fp <- mean(pos & nev[use])
    tp - fp * (pt / (1 - pt))
  })
  nb_all <- sapply(thr, function(pt) prev - (1 - prev) * (pt / (1 - pt)))
  man <- data.frame(threshold = thr, none = 0, all = nb_all,
                    clinical = nb(df$p_clin), combined = nb(df$p_comb))
  write.csv(man, file.path(outdir, "dca_netbenefit.csv"), row.names = FALSE)
  png(file.path(outdir, "DCA_3yr.png"), width = 1600, height = 1200, res = 200)
  matplot(man$threshold, man[, c("none","all","clinical","combined")],
          type = "l", lty = 1, lwd = 2,
          col = c("grey50","black","#1f77b4","#d62728"),
          xlab = "Threshold probability", ylab = "Net benefit",
          main = "Decision Curve Analysis @ 3 years (TCGA-STAD)")
  abline(h = 0, col = "grey70", lty = 3)
  legend("topright", c("Treat none","Treat all","Clinical","Combined"),
         col = c("grey50","black","#1f77b4","#d62728"), lwd = 2, bty = "n")
  dev.off()
  dca_done <- TRUE
}
cat(sprintf("\nDCA saved (dcurves=%s)\n", USE_dcurves))

## ---------------------------------------------------------------------------
## 6. Combined nomogram (rms) + calibration @ 1/3/5 years
## ---------------------------------------------------------------------------
dd <- datadist(df); options(datadist = "dd")
units(df$time) <- "Day"
cph_comb <- cph(Surv(time, event) ~ Age + Stage + Grade + risk, data = df,
                x = TRUE, y = TRUE, surv = TRUE, time.inc = t3)
surv_fun <- Survival(cph_comb)
nom <- nomogram(cph_comb,
                fun = list(function(x) surv_fun(365.25, x),
                           function(x) surv_fun(t3, x),
                           function(x) surv_fun(5 * 365.25, x)),
                funlabel = c("1-yr OS", "3-yr OS", "5-yr OS"),
                fun.at = seq(0.1, 0.9, 0.1))
png(file.path(outdir, "combined_nomogram.png"), width = 2400, height = 1500, res = 200)
plot(nom, xfrac = 0.35); dev.off()
pdf(file.path(outdir, "combined_nomogram.pdf"), width = 12, height = 7.5)
plot(nom, xfrac = 0.35); dev.off()

png(file.path(outdir, "calibration_combined.png"), width = 1800, height = 1500, res = 200)
plot(0, 0, type = "n", xlim = c(0,1), ylim = c(0,1),
     xlab = "Nomogram-predicted OS", ylab = "Observed OS",
     main = "Combined model calibration (TCGA-STAD)")
abline(0, 1, col = "grey60", lty = 2)
cols <- c("#1f77b4","#d62728","#2ca02c"); tt <- c(365.25, t3, 5*365.25)
for (i in seq_along(tt)) {
  fitc <- cph(Surv(time, event) ~ Age + Stage + Grade + risk, data = df,
              x = TRUE, y = TRUE, surv = TRUE, time.inc = tt[i])
  cal <- tryCatch(calibrate(fitc, cmethod = "KM", u = tt[i], m = 60, B = 200),
                  error = function(e) NULL)
  if (!is.null(cal)) {
    ok <- !is.na(cal[, "KM"])
    lines(cal[ok, "mean.predicted"], cal[ok, "KM"], col = cols[i], lwd = 2, type = "b")
  }
}
legend("topleft", c("1-yr","3-yr","5-yr"), col = cols, lwd = 2, bty = "n")
dev.off()
options(datadist = NULL)

## ---------------------------------------------------------------------------
## 7. ACRG external validation of transportable combined model
##    (Age + Stage + risk; ACRG lacks Grade -> excluded, noted)
## ---------------------------------------------------------------------------
load("data/geo/GSE62254.rda")                   # GSE62254.expr, GSE62254.subtype
st <- GSE62254.subtype
exA <- GSE62254.expr[names(coefs), , drop = FALSE]
ZA  <- zscore_rows(exA)
risk_acrg <- as.numeric(coefs[rownames(ZA)] %*% ZA)

acrg <- data.frame(
  time  = st$OS.m * 30.4375,                     # months -> days
  event = as.integer(st$Death),
  Age   = as.numeric(st$age),
  Stage = factor(st$Stage, levels = c("I","II","III","IV")),
  risk  = risk_acrg)
ccA <- complete.cases(acrg) & acrg$time > 0
acrg <- acrg[ccA, ]
cat(sprintf("\nACRG complete cases: n=%d, events=%d\n", nrow(acrg), sum(acrg$event)))

# Transport TCGA reduced combined coefs (no grade) -> ACRG linear predictor
m_comb_ng <- coxph(Surv(time, event) ~ Age + Stage + risk, data = df)
lp_acrg   <- predict(m_comb_ng, newdata = acrg, type = "lp")
c_acrg_comb <- cidx(acrg$time, acrg$event, lp_acrg)
c_acrg_sig  <- cidx(acrg$time, acrg$event, acrg$risk)
# refit-in-ACRG combined (upper bound) for reference
m_acrg_ref  <- coxph(Surv(time, event) ~ Age + Stage + risk, data = acrg)
c_acrg_refit <- as.numeric(summary(m_acrg_ref)$concordance[1])
cat(sprintf("ACRG C-index: signature=%.3f  combined(transported)=%.3f  combined(refit)=%.3f\n",
            c_acrg_sig, c_acrg_comb, c_acrg_refit))

## ---------------------------------------------------------------------------
## Write output tables
## ---------------------------------------------------------------------------
# External validation: out-of-sample C-index, no optimism correction needed
acrg_rows <- data.frame(
  model = c("Signature (risk)", "Combined (Age+Stage+risk, transported)",
            "Combined (Age+Stage+risk, refit)"),
  cohort = "ACRG/GSE62254", n = nrow(acrg), events = sum(acrg$event),
  C_apparent  = c(c_acrg_sig, c_acrg_comb, c_acrg_refit),
  C_corrected = c(c_acrg_sig, c_acrg_comb, c_acrg_refit),
  C_lo = NA, C_hi = NA)
cindex_all <- rbind(cindex_tab, acrg_rows)
write.csv(cindex_all, file.path(outdir, "cindex_comparison.csv"), row.names = FALSE)

added <- data.frame(
  statistic = c("deltaC_corrected(combined-clinical)", "LRT_chisq", "LRT_df",
                "LRT_p"),
  value = c(dC, lrt_chisq, lrt$Df[2], lrt_p))
write.csv(added, file.path(outdir, "added_value_stats.csv"), row.names = FALSE)
write.csv(nri_idi, file.path(outdir, "nri_idi.csv"), row.names = FALSE)
if (nrow(auc_tab)) write.csv(auc_tab, file.path(outdir, "time_dependent_AUC.csv"),
                             row.names = FALSE)

cat("\nAll outputs written to", outdir, "\n")
sink(file.path(outdir, "SUMMARY.txt"))
cat("== C-index comparison ==\n"); print(cindex_all)
cat("\n== Added value ==\n"); print(added); print(nri_idi)
cat("\n== Time-dependent AUC ==\n"); print(auc_tab)
cat(sprintf("\ntimeROC=%s survIDINRI=%s dcurves=%s survivalROC=%s\n",
            USE_timeROC, USE_survIDINRI, USE_dcurves, USE_survROC))
sink()

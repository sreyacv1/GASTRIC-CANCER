#!/usr/bin/env Rscript
# 32_nested_cv_signature.R
# UPGRADE 1: Leakage-free NESTED cross-validation of the 25-gene signature.
#
# The published "optimism-corrected" C only bootstrapped PRECOMPUTED risk
# scores -- gene screening + LASSO sat OUTSIDE resampling, so discrimination
# was optimistic. Here the ENTIRE signature-building procedure is re-run from
# scratch inside every outer training fold, and performance is measured only
# on untouched outer-fold predictions.
#
# Outer : 20 repeats of event-stratified 5-fold CV on TCGA tumours.
# Inside each outer TRAINING fold (nothing from the held-out fold is used):
#   (a) z-score mean/SD estimated on TRAINING only
#   (b) univariable Cox screen (p<0.05) on TRAINING only
#   (c) inner 5-fold cv.glmnet Cox LASSO (lambda.min AND lambda.1se)
#   (d) fit; (e) predict held-out linear predictors + survival.
# Reported from pooled out-of-fold predictions: Harrell C + Uno C (bootstrap
# 95% CI), time-dependent AUC @1/3/5y, integrated Brier score, model-size
# distribution, per-gene selection frequency, lambda.1se sensitivity.
#
# HARD RULE: real data only. No simulation / imputation.

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(survivalROC)
  library(pec)
})
set.seed(1105)

outdir <- "results/nested_cv"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

N_REPEAT <- 20
K_OUTER  <- 5
K_INNER  <- 5
P_SCREEN <- 0.05
DAY_YR   <- 365.25
tpts     <- c(1, 3, 5) * DAY_YR                 # 1/3/5 yr in days

## ---------------------------------------------------------------------------
## 1. TCGA tumours with usable OS + transferable candidate genes
##    (candidate set mirrors analysis/07: intersect TCGA & ACRG gene *names*
##     only -- uses no ACRG outcomes/values, hence no leakage).
## ---------------------------------------------------------------------------
load("results/rdata/tcga_processed.RData")       # col_data, tcga_vst
stopifnot(identical(colnames(tcga_vst), rownames(col_data)))
is_tumor <- col_data$status == "Tumor"
cd   <- col_data[is_tumor, ]
expr <- tcga_vst[, is_tumor, drop = FALSE]
OS_time  <- ifelse(cd$vital_status == "Dead",
                   cd$days_to_death, cd$days_to_last_follow_up)
OS_event <- as.integer(cd$vital_status == "Dead")
keep <- !is.na(OS_time) & OS_time > 0 & !is.na(OS_event)
expr <- expr[, keep, drop = FALSE]
OS_time <- OS_time[keep]; OS_event <- OS_event[keep]

load("data/geo/GSE62254.rda")                    # GSE62254.expr (names only)
candidates <- intersect(rownames(expr), rownames(GSE62254.expr))
expr <- expr[candidates, , drop = FALSE]
n <- ncol(expr)
cat(sprintf("TCGA tumours usable OS: n=%d events=%d; candidate genes=%d\n",
            n, sum(OS_event), length(candidates)))

surv_all <- Surv(OS_time, OS_event)

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
make_folds <- function(event, K) {              # event-stratified folds
  fold <- integer(length(event))
  for (cl in unique(event)) {
    idx <- sample(which(event == cl))
    fold[idx] <- rep_len(seq_len(K), length(idx))
  }
  fold
}

# Build signature from scratch on a training set; returns z-params + coef sets
build_model <- function(tr, ytime, yevent) {
  mu  <- rowMeans(tr)
  sdv <- apply(tr, 1, sd)
  g0  <- sdv > 0 & is.finite(sdv)
  tr  <- tr[g0, , drop = FALSE]; mu <- mu[g0]; sdv <- sdv[g0]
  Ztr <- (tr - mu) / sdv                         # genes x samples, z on train
  sv  <- Surv(ytime, yevent)
  pv  <- apply(Ztr, 1, function(g) {
    f <- tryCatch(coxph(sv ~ g), error = function(e) NULL)
    if (is.null(f)) NA_real_ else summary(f)$coefficients[1, "Pr(>|z|)"]
  })
  sig <- names(pv)[!is.na(pv) & pv < P_SCREEN]
  if (length(sig) < 2) return(NULL)
  X   <- t(Ztr[sig, , drop = FALSE])
  cvf <- tryCatch(cv.glmnet(X, sv, family = "cox", alpha = 1, nfolds = K_INNER,
                            cox.ties = "breslow"),
                  error = function(e) NULL)
  if (is.null(cvf)) return(NULL)
  grab <- function(lc) {
    cf <- as.matrix(coef(cvf, s = cvf[[lc]]))
    cf <- cf[cf[, 1] != 0, , drop = FALSE]
    setNames(cf[, 1], rownames(cf))
  }
  list(mu = mu, sdv = sdv,
       lambda.min = grab("lambda.min"),
       lambda.1se = grab("lambda.1se"))
}

# Linear predictor for new samples using TRAINING z-params
lp_new <- function(model, coefs, newexpr) {
  g <- names(coefs)
  if (length(g) == 0) return(rep(0, ncol(newexpr)))
  Z <- (newexpr[g, , drop = FALSE] - model$mu[g]) / model$sdv[g]
  as.numeric(coefs[g] %*% Z)
}

## ---------------------------------------------------------------------------
## 2. Nested CV loop
## ---------------------------------------------------------------------------
gtimes <- seq(90, min(5 * DAY_YR, quantile(OS_time[OS_event == 1], 0.95)),
              length.out = 40)                   # IBS integration grid (days)

# per-repeat storage
oof_lp_min  <- matrix(NA_real_, N_REPEAT, n)     # z-scored later
oof_lp_1se  <- matrix(NA_real_, N_REPEAT, n)
rep_C_min   <- rep(NA_real_, N_REPEAT)           # Harrell, per repeat (min)
rep_Cu_min  <- rep(NA_real_, N_REPEAT)           # Uno, per repeat (min)
rep_C_1se   <- rep(NA_real_, N_REPEAT)
rep_auc_min <- matrix(NA_real_, N_REPEAT, length(tpts))
rep_ibs_min <- rep(NA_real_, N_REPEAT)
rep_brier   <- matrix(NA_real_, N_REPEAT, length(tpts))
model_sizes <- integer(0)                        # sizes of lambda.min models
sel_counts  <- setNames(numeric(length(candidates)), candidates)
n_models    <- 0L

z <- function(v) {                                # robust standardise
  v <- v[is.finite(v)]
  s <- stats::sd(v)
  if (length(v) < 2 || !is.finite(s) || s == 0) return(v - mean(v))
  (v - mean(v)) / s
}

for (r in seq_len(N_REPEAT)) {
  fold <- make_folds(OS_event, K_OUTER)
  lp_min <- rep(NA_real_, n); lp_1se <- rep(NA_real_, n)
  Smat <- matrix(NA_real_, n, length(gtimes))    # OOF survival (lambda.min)
  for (k in seq_len(K_OUTER)) {
    tr_i <- which(fold != k); te_i <- which(fold == k)
    m <- build_model(expr[, tr_i, drop = FALSE], OS_time[tr_i], OS_event[tr_i])
    if (is.null(m)) next
    cmin <- m$lambda.min; c1se <- m$lambda.1se
    if (length(cmin) >= 1) {
      lp_min[te_i] <- lp_new(m, cmin, expr[, te_i, drop = FALSE])
      n_models <- n_models + 1L
      model_sizes <- c(model_sizes, length(cmin))
      sel_counts[names(cmin)] <- sel_counts[names(cmin)] + 1
      # baseline hazard from offset Cox on TRAINING -> OOF survival
      lp_tr <- lp_new(m, cmin, expr[, tr_i, drop = FALSE])
      cph <- tryCatch(coxph(Surv(OS_time[tr_i], OS_event[tr_i]) ~ offset(lp_tr)),
                      error = function(e) NULL)
      if (!is.null(cph)) {
        bh <- basehaz(cph, centered = FALSE)      # H0(t) at lp=0
        H0 <- stepfun(bh$time, c(0, bh$hazard))
        lp_te <- lp_min[te_i]
        Smat[te_i, ] <- exp(-outer(exp(lp_te), H0(gtimes)))
      }
    }
    if (length(c1se) >= 1)
      lp_1se[te_i] <- lp_new(m, c1se, expr[, te_i, drop = FALSE])
  }
  # per-repeat performance on complete OOF predictions (lambda.min)
  ok <- is.finite(lp_min)
  if (sum(ok) > 10 && length(unique(lp_min[ok])) > 1) {
    dd <- data.frame(t = OS_time[ok], e = OS_event[ok], lp = lp_min[ok])
    rep_C_min[r]  <- concordance(Surv(t, e) ~ lp, dd, reverse = TRUE)$concordance
    rep_Cu_min[r] <- tryCatch(
      concordance(Surv(t, e) ~ lp, dd, reverse = TRUE,
                  timewt = "n/G2")$concordance, error = function(e) NA_real_)
    rep_auc_min[r, ] <- sapply(tpts, function(tt)
      tryCatch(survivalROC(dd$t, dd$e, dd$lp, predict.time = tt,
                           method = "KM")$AUC, error = function(e) NA_real_))
    okS <- ok & stats::complete.cases(Smat)
    if (sum(okS) > 10) {
      pd <- data.frame(time = OS_time[okS], event = OS_event[okS])
      # pec prepends the time origin: supply S(0)=1 column + t=0
      Sp <- cbind(1, Smat[okS, , drop = FALSE]); g0 <- c(0, gtimes)
      pf <- tryCatch(pec(list(nested = Sp),
                         formula = Surv(time, event) ~ 1, data = pd,
                         times = g0, exact = FALSE, cens.model = "marginal",
                         verbose = FALSE), error = function(e) NULL)
      if (!is.null(pf)) {
        rep_ibs_min[r] <- crps(pf, times = max(g0))["nested", 1]
        # Brier at 1/3/5y via nearest grid time
        rep_brier[r, ] <- sapply(tpts, function(tt) {
          j <- which.min(abs(pf$time - tt)); pf$AppErr$nested[j] })
      }
    }
  }
  oko <- is.finite(lp_1se)
  if (sum(oko) > 10 && length(unique(lp_1se[oko])) > 1)
    rep_C_1se[r] <- concordance(Surv(OS_time[oko], OS_event[oko]) ~ lp_1se[oko],
                                reverse = TRUE)$concordance
  # store z-scored OOF LP for ensemble
  oof_lp_min[r, ok]  <- z(lp_min[ok])
  oof_lp_1se[r, oko] <- z(lp_1se[oko])
  cat(sprintf("repeat %2d/%d  C(min)=%.3f Uno=%.3f  C(1se)=%.3f  IBS=%.4f\n",
              r, N_REPEAT, rep_C_min[r], rep_Cu_min[r], rep_C_1se[r],
              rep_ibs_min[r]))
}

## ---------------------------------------------------------------------------
## 3. Ensemble OOF risk score + bootstrap CIs (Harrell & Uno)
## ---------------------------------------------------------------------------
ens_min <- colMeans(oof_lp_min, na.rm = TRUE)    # mean z-scored OOF LP / subj
ens_1se <- colMeans(oof_lp_1se, na.rm = TRUE)
okE <- is.finite(ens_min)
de  <- data.frame(t = OS_time[okE], e = OS_event[okE], lp = ens_min[okE])
Cens_h <- concordance(Surv(t, e) ~ lp, de, reverse = TRUE)$concordance
Cens_u <- concordance(Surv(t, e) ~ lp, de, reverse = TRUE,
                      timewt = "n/G2")$concordance

boot_ci <- function(df, uno = FALSE, B = 2000) {
  vals <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample(nrow(df), replace = TRUE)
    d <- df[idx, ]
    vals[b] <- tryCatch(
      if (uno) concordance(Surv(t, e) ~ lp, d, reverse = TRUE,
                           timewt = "n/G2")$concordance
      else concordance(Surv(t, e) ~ lp, d, reverse = TRUE)$concordance,
      error = function(e) NA_real_)
  }
  quantile(vals, c(.025, .975), na.rm = TRUE)
}
set.seed(1105); ci_h <- boot_ci(de, uno = FALSE)
set.seed(1105); ci_u <- boot_ci(de, uno = TRUE)
okE1 <- is.finite(ens_1se)
# lambda.1se can shrink to the NULL model (0 genes) in every fold -> no OOF
# predictions. Guard so this legitimate fragility result never crashes output.
lam1se_empty <- sum(okE1) < 10 || length(unique(ens_1se[okE1])) < 2
Cens_1se <- if (lam1se_empty) NA_real_ else
  concordance(Surv(OS_time[okE1], OS_event[okE1]) ~ ens_1se[okE1],
              reverse = TRUE)$concordance

## ---------------------------------------------------------------------------
## 4. Write outputs
## ---------------------------------------------------------------------------
pctl <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
perf <- rbind(
  data.frame(metric = "Harrell_C_ensemble", estimate = Cens_h,
             ci_low = ci_h[1], ci_high = ci_h[2],
             method = "ensemble OOF LP, 2000x subject bootstrap"),
  data.frame(metric = "Uno_C_ensemble", estimate = Cens_u,
             ci_low = ci_u[1], ci_high = ci_u[2],
             method = "ensemble OOF LP, IPCW (timewt=n/G2), 2000x bootstrap"),
  data.frame(metric = "Harrell_C_perRepeat", estimate = mean(rep_C_min, na.rm = TRUE),
             ci_low = pctl(rep_C_min)[1], ci_high = pctl(rep_C_min)[2],
             method = "mean +/- 2.5-97.5 pctile over 20 repeats"),
  data.frame(metric = "Uno_C_perRepeat", estimate = mean(rep_Cu_min, na.rm = TRUE),
             ci_low = pctl(rep_Cu_min)[1], ci_high = pctl(rep_Cu_min)[2],
             method = "mean +/- pctile over 20 repeats"),
  data.frame(metric = "IntegratedBrierScore", estimate = mean(rep_ibs_min, na.rm = TRUE),
             ci_low = pctl(rep_ibs_min)[1], ci_high = pctl(rep_ibs_min)[2],
             method = "IPCW IBS to 5y, mean +/- pctile over repeats"),
  data.frame(metric = "Harrell_C_lambda1se_ensemble", estimate = Cens_1se,
             ci_low = NA, ci_high = NA,
             method = "lambda.1se sensitivity, ensemble OOF LP"),
  data.frame(metric = "Harrell_C_lambda1se_perRepeat", estimate = mean(rep_C_1se, na.rm = TRUE),
             ci_low = pctl(rep_C_1se)[1], ci_high = pctl(rep_C_1se)[2],
             method = "lambda.1se sensitivity, mean +/- pctile"))
row.names(perf) <- NULL
write.csv(perf, file.path(outdir, "performance.csv"), row.names = FALSE)

auc_tab <- data.frame(
  t_years = c(1, 3, 5),
  AUC_mean = colMeans(rep_auc_min, na.rm = TRUE),
  AUC_lo = apply(rep_auc_min, 2, function(v) pctl(v)[1]),
  AUC_hi = apply(rep_auc_min, 2, function(v) pctl(v)[2]),
  Brier_mean = colMeans(rep_brier, na.rm = TRUE),
  method = "survivalROC KM AUC / IPCW Brier, mean over 20 repeats")
write.csv(auc_tab, file.path(outdir, "timeAUC.csv"), row.names = FALSE)

sel_freq <- data.frame(gene = names(sel_counts),
                       selection_frequency = as.numeric(sel_counts) / n_models,
                       times_selected = as.integer(sel_counts),
                       n_models = n_models)
sel_freq <- sel_freq[order(-sel_freq$selection_frequency), ]
sel_freq <- sel_freq[sel_freq$times_selected > 0, ]
write.csv(sel_freq, file.path(outdir, "gene_selection_frequency.csv"),
          row.names = FALSE)

sink(file.path(outdir, "SUMMARY.txt"))
cat("NESTED CROSS-VALIDATION of the 25-gene prognostic signature\n")
cat("TCGA-STAD tumours: n=", n, " events=", sum(OS_event),
    "; candidate genes=", length(candidates), "\n", sep = "")
cat("Design: 20 repeats x 5-fold outer CV; ALL steps (z-scoring, univariable\n")
cat("Cox screen p<0.05, inner 5-fold cv.glmnet Cox LASSO) re-run inside each\n")
cat("outer training fold; performance from untouched out-of-fold predictions.\n\n")
cat("== Honest (nested) discrimination ==\n")
cat(sprintf("Harrell C = %.3f (95%% CI %.3f-%.3f)\n", Cens_h, ci_h[1], ci_h[2]))
cat(sprintf("Uno C     = %.3f (95%% CI %.3f-%.3f)\n", Cens_u, ci_u[1], ci_u[2]))
cat(sprintf("Per-repeat Harrell C: mean=%.3f (range %.3f-%.3f over 20 repeats)\n",
            mean(rep_C_min, na.rm = TRUE), min(rep_C_min, na.rm = TRUE),
            max(rep_C_min, na.rm = TRUE)))
cat(sprintf("Integrated Brier Score (to 5y) = %.4f\n", mean(rep_ibs_min, na.rm = TRUE)))
cat("\n== Time-dependent AUC (mean over repeats) ==\n"); print(auc_tab[, 1:5])
cat("\n== lambda.1se sensitivity ==\n")
if (lam1se_empty) {
  cat("lambda.1se shrank to the NULL model (0 genes selected) in every outer\n")
  cat("fold -> no out-of-fold discrimination (C undefined, ~0.5). This is an\n")
  cat("HONEST fragility signal: under the more conservative 1SE rule the LASSO\n")
  cat("retains no genes, so the signature's apparent value rests entirely on\n")
  cat("the less-penalised lambda.min.\n")
} else {
  cat(sprintf("Harrell C (1se) = %.3f (per-repeat mean %.3f)\n",
              Cens_1se, mean(rep_C_1se, na.rm = TRUE)))
}
cat("\n== Model size distribution (lambda.min, ", n_models, " fold-models) ==\n", sep = "")
print(summary(model_sizes))
cat(sprintf("median genes selected = %g (IQR %g-%g)\n", median(model_sizes),
            quantile(model_sizes, .25), quantile(model_sizes, .75)))
cat("\n== Top 15 genes by selection frequency ==\n")
print(head(sel_freq, 15), row.names = FALSE)
cat("\nCONTEXT: the published apparent/optimism-corrected C was ~0.72. The\n")
cat("nested C above is the honest, un-inflated discrimination.\n")
sink()

cat("\n=== DONE. Outputs in", outdir, "===\n")
cat(sprintf("Nested Harrell C=%.3f (%.3f-%.3f); Uno C=%.3f (%.3f-%.3f)\n",
            Cens_h, ci_h[1], ci_h[2], Cens_u, ci_u[1], ci_u[2]))

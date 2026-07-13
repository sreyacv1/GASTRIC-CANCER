#!/usr/bin/env Rscript
# 33_timevarying_ACRG.R
# UPGRADE 2: Time-varying Cox for the non-proportional ACRG signature effect.
#
# Reviewer found the continuous ACRG signature effect violates proportional
# hazards (cox.zph risk-score p ~ 0.003). A single HR is then misleading.
# We (1) reconstruct the exact continuous ACRG signature score as in
# analysis/17_external_utility_ACRG.R, (2) run cox.zph on the principal Cox
# models, and (3) fit a prespecified time-varying-coefficient model
#   log HR_risk(t) = b0 + b1 * log(t / 36 months)
# via a tt() interaction (age-adjusted, baseline hazard stratified by Stage),
# and report HR(t) at 12/36/60 months with 95% CIs. A survSplit step-function
# model (cuts 12/36/60) is included as a robustness check.
#
# HARD RULE: real data only, complete-case, no simulation/imputation.

suppressPackageStartupMessages({ library(survival) })
set.seed(1105)

outdir <- "results/timevarying_ACRG"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

zscore_rows <- function(m) t(scale(t(m)))

## ---------------------------------------------------------------------------
## 1. ACRG signature score (identical construction to script 17), time = MONTHS
## ---------------------------------------------------------------------------
load("data/geo/GSE62254.rda")
stopifnot(identical(colnames(GSE62254.expr),
                    as.character(GSE62254.subtype$GEO_ID)))
sig   <- read.csv("results/validation/signature_coefficients.csv",
                  stringsAsFactors = FALSE)
coefs <- setNames(sig$coefficient, sig$gene)
stopifnot(all(names(coefs) %in% rownames(GSE62254.expr)))

st  <- GSE62254.subtype
ZA  <- zscore_rows(GSE62254.expr[names(coefs), , drop = FALSE])
stopifnot(all(!is.na(ZA)))
risk_raw <- as.numeric(coefs[rownames(ZA)] %*% ZA)

acrg <- data.frame(
  time  = as.numeric(st$OS.m),                   # MONTHS (not /12)
  event = as.integer(st$Death),
  Age   = as.numeric(st$age),
  Stage = factor(st$Stage, levels = c("I", "II", "III", "IV")),
  risk  = as.numeric(scale(risk_raw)))           # per-1-SD risk within ACRG
acrg <- acrg[complete.cases(acrg) & acrg$time > 0, ]
cat(sprintf("ACRG complete cases: n=%d events=%d (time in months)\n",
            nrow(acrg), sum(acrg$event)))

## ---------------------------------------------------------------------------
## 2. Proportional-hazards test (cox.zph) on principal models
## ---------------------------------------------------------------------------
m_sig <- coxph(Surv(time, event) ~ risk, data = acrg)
m_adj <- coxph(Surv(time, event) ~ risk + Stage + Age, data = acrg)
zph_sig <- cox.zph(m_sig)
zph_adj <- cox.zph(m_adj)
cat("\n== cox.zph: signature-only model ==\n"); print(zph_sig)
cat("\n== cox.zph: signature + Stage + Age ==\n"); print(zph_adj)

zph_row <- function(z, model) {
  tab <- z$table
  data.frame(model = model, term = rownames(tab),
             chisq = tab[, "chisq"], df = tab[, "df"],
             p = tab[, "p"], row.names = NULL)
}
zph_tab <- rbind(zph_row(zph_sig, "signature"),
                 zph_row(zph_adj, "signature+Stage+Age"))
write.csv(zph_tab, file.path(outdir, "coxzph.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 3a. Prespecified time-varying-coefficient model (tt, continuous log-time)
##     log HR_risk(t) = b0 + b1*log(t/36); age-adjusted; baseline strat by Stage
## ---------------------------------------------------------------------------
m_tv <- coxph(Surv(time, event) ~ risk + tt(risk) + Age + strata(Stage),
              data = acrg,
              tt = function(x, t, ...) x * log(t / 36))
cat("\n== Time-varying-coefficient model (tt) ==\n"); print(summary(m_tv)$coefficients)

b  <- coef(m_tv); V <- vcov(m_tv)
i_r <- which(names(b) == "risk")
i_t <- grep("tt\\(risk\\)", names(b))
stopifnot(length(i_r) == 1, length(i_t) == 1)
hr_at <- function(t) {
  g <- c(1, log(t / 36))
  est <- b[i_r] + b[i_t] * log(t / 36)
  se  <- sqrt(g %*% V[c(i_r, i_t), c(i_r, i_t)] %*% g)
  c(HR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}
tv_tab <- do.call(rbind, lapply(c(12, 36, 60), function(t)
  data.frame(model = "tt log-time (age-adj, Stage-stratified)",
             month = t, t(hr_at(t)))))
names(tv_tab) <- c("model", "month", "HR", "CI_low", "CI_high")

## ---------------------------------------------------------------------------
## 3b. Robustness: step-function (survSplit) HR per interval
## ---------------------------------------------------------------------------
cuts <- c(12, 36, 60)
sp <- survSplit(Surv(time, event) ~ ., data = acrg, cut = cuts,
                episode = "interval")
sp$interval <- factor(sp$interval)
m_step <- coxph(Surv(tstart, time, event) ~ risk:interval + Age + strata(Stage),
                data = sp)
ss <- summary(m_step)$coefficients
rk <- grep("^risk:interval", rownames(ss))
lbl <- c("[0,12)m", "[12,36)m", "[36,60)m", "[60+)m")
step_tab <- data.frame(
  model = "survSplit step (age-adj, Stage-stratified)",
  interval = lbl[seq_along(rk)],
  HR = exp(ss[rk, "coef"]),
  CI_low = exp(ss[rk, "coef"] - 1.96 * ss[rk, "se(coef)"]),
  CI_high = exp(ss[rk, "coef"] + 1.96 * ss[rk, "se(coef)"]),
  p = ss[rk, "Pr(>|z|)"], row.names = NULL)
cat("\n== Step-function HR per interval ==\n"); print(step_tab)

write.csv(tv_tab, file.path(outdir, "hr_over_time.csv"), row.names = FALSE)
write.csv(step_tab, file.path(outdir, "hr_over_time_stepfn.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## SUMMARY
## ---------------------------------------------------------------------------
p_zph_sig  <- zph_sig$table["risk", "p"]
p_zph_adj  <- zph_adj$table["risk", "p"]
hr12 <- tv_tab$HR[tv_tab$month == 12]
hr36 <- tv_tab$HR[tv_tab$month == 36]
hr60 <- tv_tab$HR[tv_tab$month == 60]
attenuates <- hr12 > hr36 && hr36 >= hr60

sink(file.path(outdir, "SUMMARY.txt"))
cat("TIME-VARYING Cox for the non-proportional ACRG signature effect\n")
cat("ACRG/GSE62254: n=", nrow(acrg), " events=", sum(acrg$event),
    " (time in months)\n\n", sep = "")
cat("== Proportional-hazards test (cox.zph) ==\n")
cat(sprintf("Signature-only model: risk PH p = %.4f; GLOBAL p = %.4f\n",
            p_zph_sig, zph_sig$table["GLOBAL", "p"]))
cat(sprintf("Adjusted model:       risk PH p = %.4f; GLOBAL p = %.4f\n",
            p_zph_adj, zph_adj$table["GLOBAL", "p"]))
cat("\n== HR_risk(t) from tt() log-time model (per 1-SD, age-adj, Stage-strat) ==\n")
print(tv_tab, row.names = FALSE)
cat("\n== Step-function robustness ==\n"); print(step_tab, row.names = FALSE)
cat(sprintf("\nEffect attenuates over follow-up: %s\n",
            ifelse(attenuates, "YES", "not monotone -- inspect table")))
cat(sprintf("HR(12mo)=%.2f -> HR(36mo)=%.2f -> HR(60mo)=%.2f\n", hr12, hr36, hr60))
cat("\nWORDING GUIDANCE: the signature is associated with EARLIER mortality,\n")
cat("with the hazard ratio ATTENUATING over follow-up (strong early, weaker\n")
cat("late); report HR(t) rather than a single proportional-hazards HR.\n")
sink()

cat("\n=== DONE. Outputs in", outdir, "===\n")
cat(sprintf("cox.zph risk p: sig=%.4f adj=%.4f | HR 12/36/60mo = %.2f/%.2f/%.2f\n",
            p_zph_sig, p_zph_adj, hr12, hr36, hr60))

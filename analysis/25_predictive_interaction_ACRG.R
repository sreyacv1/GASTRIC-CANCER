#!/usr/bin/env Rscript
# 25_predictive_interaction_ACRG.R
# PROGNOSTIC -> PREDICTIVE pivot for the 25-gene stromal/CAF signature.
#
# The signature adds no prognostic value over stage (shown elsewhere).
# Question here: is it PREDICTIVE, i.e. does its effect on OS differ by
#   (a) ACRG molecular subtype, (b) stage (I-II vs III-IV),
#   (c) adjuvant chemotherapy (if recoverable)?
# A significant interaction = genuinely predictive/actionable. A null = honest
# negative. We report whichever it is.
#
# Signature scoring is IDENTICAL to 07_external_validation.R /
# 17_external_utility_ACRG.R: fixed TCGA LASSO-Cox genes+coefs, genes z-scored
# WITHIN ACRG, linear predictor = coefs %*% Z. Risk is then z-scored across the
# whole cohort so per-SD HRs are comparable across subgroups.
#
# HARD RULE: real data only, complete-case, no simulation/imputation.

suppressPackageStartupMessages({ library(survival) })
set.seed(42)  # deterministic anyway (Cox LRT); set for reproducibility

outdir <- "results/predictive"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

zscore_rows <- function(m) t(scale(t(m)))   # z-score each gene (row)

## ---------------------------------------------------------------------------
## 1. ACRG cohort + fixed TCGA signature -> risk score (as in 07/17)
## ---------------------------------------------------------------------------
load("data/geo/GSE62254.rda")               # GSE62254.expr, GSE62254.subtype
stopifnot(identical(colnames(GSE62254.expr),
                    as.character(GSE62254.subtype$GEO_ID)))
sig   <- read.csv("results/validation/signature_coefficients.csv",
                  stringsAsFactors = FALSE)
coefs <- setNames(sig$coefficient, sig$gene)
stopifnot(all(names(coefs) %in% rownames(GSE62254.expr)))

st <- GSE62254.subtype
ZA <- zscore_rows(GSE62254.expr[names(coefs), , drop = FALSE])  # z within ACRG
stopifnot(all(!is.na(ZA)))                  # all signature genes non-constant
risk_raw <- as.numeric(coefs[rownames(ZA)] %*% ZA)

# ACRG molecular subtype: labels in .rda are EMT / MSI / TP53neg / TP53positive
# = the four ACRG subtypes MSS/EMT, MSI, MSS/TP53-, MSS/TP53+ (Cristescu 2015).
sub_map <- c(EMT = "MSS/EMT", MSI = "MSI",
             TP53neg = "MSS/TP53-", TP53positive = "MSS/TP53+")
subtype <- factor(sub_map[st$ACRG.sub],
                  levels = c("MSI", "MSS/TP53+", "MSS/TP53-", "MSS/EMT"))

acrg <- data.frame(
  time    = as.numeric(st$OS.m),            # months (HR invariant to time unit)
  event   = as.integer(st$Death),
  risk    = as.numeric(scale(risk_raw)),    # per-SD, scaled across whole cohort
  subtype = subtype,
  stage4  = factor(st$Stage, levels = c("I", "II", "III", "IV")),
  stringsAsFactors = FALSE)
acrg$stage_bin <- factor(ifelse(acrg$stage4 %in% c("I", "II"), "I-II", "III-IV"),
                         levels = c("I-II", "III-IV"))

cc <- !is.na(acrg$time) & acrg$time > 0 & !is.na(acrg$event) &
      !is.na(acrg$risk) & !is.na(acrg$subtype) & !is.na(acrg$stage_bin)
acrg <- acrg[cc, ]
cat(sprintf("ACRG complete cases: n=%d, events=%d\n",
            nrow(acrg), sum(acrg$event)))
cat("Subtype n:  "); print(table(acrg$subtype))
cat("Stage bin n:"); print(table(acrg$stage_bin))

## ---------------------------------------------------------------------------
## 2. Interaction tests (likelihood-ratio test of the interaction term)
##    (a) risk x subtype   (b) risk x stage_bin   (c) risk x adjuvant chemo
## ---------------------------------------------------------------------------
lrt_interaction <- function(data, xvar) {
  # main-effects vs main-effects + risk:xvar, LRT on the added interaction df
  f0 <- as.formula(sprintf("Surv(time, event) ~ risk + %s", xvar))
  f1 <- as.formula(sprintf("Surv(time, event) ~ risk * %s", xvar))
  m0 <- coxph(f0, data = data); m1 <- coxph(f1, data = data)
  a  <- anova(m0, m1, test = "LRT")
  list(chisq = a$Chisq[2], df = a$Df[2], p = a$`Pr(>|Chi|)`[2],
       m0 = m0, m1 = m1)
}

int_sub   <- lrt_interaction(acrg, "subtype")
int_stage <- lrt_interaction(acrg, "stage_bin")

# (c) adjuvant chemotherapy: NOT deposited in GSE62254 (GEO series-matrix
# characteristics carry only 'tissue' and 'patient'; getGEO pData is parsed
# from the same series matrix, so the field is genuinely unavailable). Skip.
chemo_available <- FALSE

int_tab <- data.frame(
  interaction = c("risk x molecular_subtype (4 ACRG subtypes)",
                  "risk x stage (I-II vs III-IV)",
                  "risk x adjuvant_chemotherapy"),
  test        = c("LRT (Cox)", "LRT (Cox)", "not tested"),
  LRT_chisq   = c(int_sub$chisq, int_stage$chisq, NA_real_),
  df          = c(int_sub$df,    int_stage$df,    NA_real_),
  p_value     = c(int_sub$p,     int_stage$p,     NA_real_),
  note        = c("", "",
    "adjuvant-chemo status not deposited in GSE62254 GEO characteristics; unavailable"),
  row.names = NULL)
cat("\n== Interaction tests ==\n"); print(int_tab)
write.csv(int_tab, file.path(outdir, "interaction_tests.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. Per-subgroup signature HRs (Cox of risk within each subgroup)
##    risk is on the same whole-cohort SD scale -> HRs directly comparable.
## ---------------------------------------------------------------------------
subgroup_HR <- function(data, label, group) {
  d <- data[data[[group]] == label, ]
  n <- nrow(d); ev <- sum(d$event)
  fit <- tryCatch(coxph(Surv(time, event) ~ risk, data = d),
                  error = function(e) NULL)
  if (is.null(fit) || any(is.na(coef(fit)))) {
    return(data.frame(grouping = group, subgroup = label, n = n, events = ev,
                      HR_per_SD = NA, CI_low = NA, CI_high = NA, p_value = NA))
  }
  s <- summary(fit)
  data.frame(grouping = group, subgroup = label, n = n, events = ev,
             HR_per_SD = s$conf.int[1, "exp(coef)"],
             CI_low    = s$conf.int[1, "lower .95"],
             CI_high   = s$conf.int[1, "upper .95"],
             p_value   = s$coefficients[1, "Pr(>|z|)"],
             row.names = NULL)
}

sub_rows <- do.call(rbind, lapply(levels(acrg$subtype),
                                  subgroup_HR, data = acrg, group = "subtype"))
stg_rows <- do.call(rbind, lapply(levels(acrg$stage_bin),
                                  subgroup_HR, data = acrg, group = "stage_bin"))
overall  <- {
  fit <- coxph(Surv(time, event) ~ risk, data = acrg); s <- summary(fit)
  data.frame(grouping = "overall", subgroup = "ALL",
             n = nrow(acrg), events = sum(acrg$event),
             HR_per_SD = s$conf.int[1, "exp(coef)"],
             CI_low = s$conf.int[1, "lower .95"],
             CI_high = s$conf.int[1, "upper .95"],
             p_value = s$coefficients[1, "Pr(>|z|)"], row.names = NULL)
}
subgroup_tab <- rbind(overall, sub_rows, stg_rows)
cat("\n== Per-subgroup signature HRs (per SD of risk) ==\n"); print(subgroup_tab)
write.csv(subgroup_tab, file.path(outdir, "subgroup_HRs.csv"), row.names = FALSE)

emt <- subgroup_tab[subgroup_tab$subgroup == "MSS/EMT", ]

## ---------------------------------------------------------------------------
## 4. SUMMARY / verdict
## ---------------------------------------------------------------------------
alpha <- 0.05
sig_sub   <- !is.na(int_sub$p)   && int_sub$p   < alpha
sig_stage <- !is.na(int_stage$p) && int_stage$p < alpha
any_sig   <- sig_sub || sig_stage

sink(file.path(outdir, "SUMMARY.txt"))
cat("PROGNOSTIC -> PREDICTIVE test of the 25-gene stromal/CAF signature\n")
cat("Cohort: ACRG / GSE62254 (independent). n=", nrow(acrg),
    ", events=", sum(acrg$event), "\n", sep = "")
cat("Signature: genes+coefs FIXED from TCGA; genes z-scored within ACRG;\n")
cat("risk z-scored across whole cohort (HR reported per 1 SD of risk).\n\n")

cat("== Interaction tests (LRT of interaction term) ==\n")
print(int_tab); cat("\n")

cat("== Per-subgroup signature HRs (per SD) ==\n")
print(subgroup_tab); cat("\n")

cat("== VERDICT ==\n")
if (any_sig) {
  cat("PREDICTIVE signal detected: at least one interaction is significant.\n")
  if (sig_sub)   cat(sprintf(" - risk x subtype: LRT p=%.3g (SIGNIFICANT)\n", int_sub$p))
  if (sig_stage) cat(sprintf(" - risk x stage:   LRT p=%.3g (SIGNIFICANT)\n", int_stage$p))
} else {
  cat("ALL NULL: no significant interaction. The signature is neither\n")
  cat("prognostic-above-stage nor predictive (effect does not differ by\n")
  cat("subtype or stage). Honest negative.\n")
}
cat(sprintf("\nHypothesis check (CAF signature strongest in MSS/EMT):\n"))
cat(sprintf(" MSS/EMT subgroup: n=%d, events=%d, HR/SD=%.2f (95%% CI %.2f-%.2f), p=%.3g\n",
            emt$n, emt$events, emt$HR_per_SD, emt$CI_low, emt$CI_high, emt$p_value))
cat(sprintf(" risk x subtype interaction p=%.3g -> %s\n", int_sub$p,
            ifelse(sig_sub, "subtype-dependent effect", "no evidence effect differs by subtype")))
sink()

cat("\nInteraction p (subtype)=", signif(int_sub$p, 3),
    "  interaction p (stage)=", signif(int_stage$p, 3), "\n")
cat("MSS/EMT HR/SD=", round(emt$HR_per_SD, 3), " p=", signif(emt$p_value, 3), "\n")
cat("Outputs in", outdir, "\n")

#!/usr/bin/env Rscript
# =============================================================================
#  NOMOGRAM FOR OVERALL SURVIVAL — TCGA-STAD (HONEST / REAL-COVARIATE REBUILD)
#  ---------------------------------------------------------------------------
#  This replaces nomogram_OS.R for publication. Differences vs the original:
#    * NO fabricated immune covariates. The original invented CD8+/Macrophage
#      fractions with rnorm() (orig L115-120); those are removed entirely.
#    * NO random imputation. The original filled missing Lauren/Stage/Grade and
#      median-imputed TMB via sample()/median (orig L113, L136-148); removed.
#      Missing covariates are handled by honest complete-case analysis.
#    * Model is built ONLY on real TCGA columns: Age, TNM Stage, Grade, TMB,
#      and real HAT1 expression (top WGCNA hub). Lauren is excluded from the
#      survival model (36% missing) but retained in the DEG analysis elsewhere.
#  Everything else (stepAIC selection, 500-bootstrap optimism-corrected
#  C-index, rms/ggplot nomogram, calibration) is retained unchanged.
# =============================================================================

cat("=================================================================\n")
cat("  Gastric Cancer OS Nomogram — REAL-COVARIATE REBUILD\n")
cat("  Started:", format(Sys.time()), "\n")
cat("=================================================================\n\n")

LIB  <- file.path(getwd(), "r_env/lib/R/library")
REPO <- "https://cloud.r-project.org"
auto_install <- function(pkg, lib = LIB, repos = REPO) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  Installing: ", pkg)
    install.packages(pkg, lib = lib, repos = repos, quiet = TRUE, dependencies = TRUE)
  }
}
for (p in c("survival", "MASS", "dplyr", "ggplot2", "gridExtra", "scales")) auto_install(p)
suppressPackageStartupMessages({
  library(survival); library(MASS); library(dplyr); library(ggplot2)
  library(gridExtra); library(scales)
})
RMS_OK <- requireNamespace("rms", quietly = TRUE)
if (RMS_OK) {
  suppressPackageStartupMessages(library(rms))
  options(contrasts = c("contr.treatment", "contr.treatment"))
  cat("[INFO] rms available — using rms::nomogram()\n\n")
} else {
  cat("[INFO] rms not available — using custom ggplot2 nomogram\n\n")
}

OUT <- "results/nomogram"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 1: LOAD REAL TCGA-STAD DATA
# ─────────────────────────────────────────────────────────────────────────────
cat("[Step 1] Loading TCGA-STAD data from local processed RData...\n")
LOCAL_RDATA <- "results/rdata/tcga_processed.RData"
if (!file.exists(LOCAL_RDATA)) stop("Local processed RData not found at: ", LOCAL_RDATA)
load(LOCAL_RDATA)  # col_data, tcga_vst, res
cat("  Loaded. col_data:", paste(dim(col_data), collapse="x"),
    "| tcga_vst:", paste(dim(tcga_vst), collapse="x"), "\n")

is_tumor <- col_data$status == "Tumor"
cd <- col_data[is_tumor, ]
expr_hat1 <- tcga_vst["HAT1", is_tumor]

# Survival (real)
os_days   <- suppressWarnings(as.numeric(as.character(
  ifelse(cd$vital_status == "Dead", cd$days_to_death, cd$days_to_last_follow_up))))
os_status <- ifelse(cd$vital_status == "Dead", 1, 0)

# Age (real, days -> years)
age <- suppressWarnings(as.numeric(as.character(cd$age_at_diagnosis)) / 365.25)

# TNM Stage (real; NA stays NA — no imputation)
stage_raw <- cd$ajcc_pathologic_stage
stage <- dplyr::case_when(
  grepl("Stage IV",  stage_raw, ignore.case=TRUE) ~ "IV",
  grepl("Stage III", stage_raw, ignore.case=TRUE) ~ "III",
  grepl("Stage II",  stage_raw, ignore.case=TRUE) ~ "II",
  grepl("Stage I",   stage_raw, ignore.case=TRUE) ~ "I",
  TRUE ~ NA_character_)
stage <- factor(stage, levels = c("I","II","III","IV"))

# Grade (real; NA stays NA)
grade <- factor(ifelse(cd$tumor_grade %in% c("G1","G2","G3"),
                       as.character(cd$tumor_grade), NA), levels = c("G1","G2","G3"))

# TMB (real; NA stays NA — no median imputation)
tmb <- suppressWarnings(as.numeric(as.character(cd$paper_Total.Mutation.Rate)))

df <- data.frame(
  OS_days = os_days, OS_status = os_status,
  Age = age, Stage = stage, Grade = grade, TMB = tmb, HAT1_expr = expr_hat1,
  stringsAsFactors = FALSE)
df$OS_years <- df$OS_days / 365.25

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 2: HONEST COMPLETE-CASE ANALYTIC COHORT (no imputation)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 2] Building complete-case analytic cohort (no imputation)...\n")
VARS <- c("OS_years","OS_status","Age","Stage","Grade","TMB","HAT1_expr")
# Winsorise TMB at 99th pct (outlier control only — not imputation)
df$TMB <- pmin(df$TMB, quantile(df$TMB, 0.99, na.rm = TRUE))
df_c <- df[complete.cases(df[, VARS]), ]
df_c <- df_c[df_c$OS_years > 0 & !is.na(df_c$OS_years), ]
cat(sprintf("  Tumor samples: %d | Complete-case analytic N = %d | Events = %d (%.1f%%)\n",
            sum(is_tumor), nrow(df_c), sum(df_c$OS_status), 100*mean(df_c$OS_status)))
if (nrow(df_c) < 50) stop("Too few complete cases.")

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 3: MULTIVARIATE COX + BACKWARD AIC (real covariates only)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 3] Multivariate Cox + backward AIC (Age+Stage+Grade+TMB+HAT1)...\n")
full_cox <- coxph(Surv(OS_years, OS_status) ~ Age + Stage + Grade + TMB + HAT1_expr,
                  data = df_c, ties = "efron")
cat("\n--- Full model ---\n"); print(summary(full_cox))
step_cox <- stepAIC(full_cox, direction = "backward", trace = FALSE)
cat("\n--- Stepwise-selected model ---\n"); print(summary(step_cox))

sel_vars <- all.vars(formula(step_cox))
sel_vars <- sel_vars[!sel_vars %in% c("OS_years","OS_status","Surv")]
cat(sprintf("Retained (%d): %s\n", length(sel_vars), paste(sel_vars, collapse=", ")))

sc <- summary(step_cox)$coefficients
cox_tbl <- data.frame(
  Variable = rownames(sc),
  coef = round(sc[,"coef"],4), HR = round(exp(sc[,"coef"]),4),
  HR_95_lo = round(exp(confint(step_cox)[,1]),4),
  HR_95_hi = round(exp(confint(step_cox)[,2]),4),
  se_coef = round(sc[,"se(coef)"],4), z = round(sc[,"z"],3),
  p_value = round(sc[,"Pr(>|z|)"],5))
write.csv(cox_tbl, file.path(OUT, "Cox_model_summary_real.csv"), row.names = FALSE)
cat(sprintf("  Cox table -> %s\n", file.path(OUT, "Cox_model_summary_real.csv")))

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 4: BOOTSTRAP OPTIMISM-CORRECTED C-INDEX (B = 500)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 4] Bootstrap optimism-corrected C-index (Harrell, B = 500)...\n")
# Correct Harrell/Efron optimism bootstrap:
#   for each bootstrap b: fit model on bootstrap sample, then
#   optimism_b = C(boot model on boot data) - C(boot model on ORIGINAL data).
#   corrected = apparent - mean(optimism_b).  This guarantees corrected <= apparent.
set.seed(42)
B <- 500; N <- nrow(df_c)
c_apparent <- concordance(step_cox)$concordance
boot_formula <- formula(step_cox)
opt_vec <- rep(NA_real_, B)
pb_interval <- max(1L, B %/% 10L)
for (b in seq_len(B)) {
  idx <- sample(N, N, replace = TRUE); dboot <- df_c[idx, ]
  fit_b <- tryCatch(coxph(boot_formula, data = dboot, ties = "efron",
                          init = coef(step_cox), control = coxph.control(iter.max = 20)),
                    error = function(e) NULL)
  if (is.null(fit_b)) next
  c_app_b <- tryCatch(concordance(fit_b)$concordance, error = function(e) NA_real_)
  lp_orig <- tryCatch(predict(fit_b, newdata = df_c, type = "lp"), error = function(e) NULL)
  if (is.null(lp_orig) || is.na(c_app_b)) next
  c_orig_b <- concordance(Surv(df_c$OS_years, df_c$OS_status) ~ lp_orig, reverse = TRUE)$concordance
  opt_vec[b] <- c_app_b - c_orig_b
  if (b %% pb_interval == 0) cat(sprintf("    bootstrap %d/%d...\n", b, B))
}
opt_vec <- opt_vec[!is.na(opt_vec)]
optimism <- mean(opt_vec)
c_corrected <- c_apparent - optimism
# Optimism-corrected 95% CI = percentiles of (apparent - optimism_b)
corr_dist <- c_apparent - opt_vec
c_lo <- quantile(corr_dist, 0.025); c_hi <- quantile(corr_dist, 0.975)
cat("\n=====================================\n")
cat(sprintf("  Apparent C-index       : %.4f\n", c_apparent))
cat(sprintf("  Optimism               : %.4f\n", optimism))
cat(sprintf("  Bias-corrected C-index : %.4f\n", c_corrected))
cat(sprintf("  Bootstrap 95%% CI       : [%.4f, %.4f]\n", c_lo, c_hi))
cat("=====================================\n\n")
writeLines(c(
  sprintf("analytic_N,%d", nrow(df_c)),
  sprintf("events,%d", sum(df_c$OS_status)),
  sprintf("retained_vars,%s", paste(sel_vars, collapse="|")),
  sprintf("c_apparent,%.4f", c_apparent),
  sprintf("c_corrected,%.4f", c_corrected),
  sprintf("c_lo,%.4f", c_lo), sprintf("c_hi,%.4f", c_hi)),
  file.path(OUT, "Cindex_summary_real.csv"))

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 5: NOMOGRAM (rms if available, else ggplot/base fallback)
# ─────────────────────────────────────────────────────────────────────────────
cat("[Step 5] Building nomogram...\n")
nom_png <- file.path(OUT, "Nomogram_GC_OS_real.png")
nom_pdf <- file.path(OUT, "Nomogram_GC_OS_real.pdf")

if (RMS_OK) {
  dd <- rms::datadist(df_c); options(datadist = "dd")
  frm <- as.formula(paste("Surv(OS_years, OS_status) ~", paste(sel_vars, collapse=" + ")))
  fit_rms <- rms::cph(frm, data = df_c, x=TRUE, y=TRUE, surv=TRUE, time.inc=1)
  surv_fun <- Survival(fit_rms)
  atv <- c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
  nom <- rms::nomogram(fit_rms,
    fun = list(function(lp) surv_fun(1,lp), function(lp) surv_fun(3,lp), function(lp) surv_fun(5,lp)),
    fun.at = list(atv, atv, atv),
    funlabel = c("1-Year Survival","3-Year Survival","5-Year Survival"),
    maxscale = 100, lp = FALSE)
  for (ext in c("png","pdf")) {
    fpath <- if (ext=="png") nom_png else nom_pdf
    if (ext=="png") png(fpath, width=14, height=10, units="in", res=300, bg="white")
    else            pdf(fpath, width=14, height=10)
    par(mar=c(2,1,3,1))
    plot(nom, xfrac=0.35, cex.axis=0.78, cex.var=0.88, label.every=1)
    title(main=sprintf("Nomogram — Overall Survival, Gastric Cancer (TCGA-STAD)\nN=%d, %d events | Bias-corrected C-index=%.3f (95%% CI %.3f–%.3f)",
                       nrow(df_c), sum(df_c$OS_status), c_corrected, c_lo, c_hi),
          cex.main=1.05, font.main=2)
    dev.off()
  }
} else {
  # Minimal base-graphics fallback: coefficient/points summary plot.
  lp_vals <- predict(step_cox, type="lp")
  png(nom_png, width=12, height=8, units="in", res=300, bg="white")
  dotchart(sort(coef(step_cox)), main=sprintf(
    "Cox coefficients (log HR) — N=%d, %d events, C=%.3f", nrow(df_c), sum(df_c$OS_status), c_corrected),
    xlab="log hazard ratio"); abline(v=0, lty=2, col="red"); dev.off()
  pdf(nom_pdf, width=12, height=8)
  dotchart(sort(coef(step_cox)), xlab="log hazard ratio"); abline(v=0, lty=2, col="red"); dev.off()
}
cat(sprintf("  Nomogram -> %s ; %s\n", nom_png, nom_pdf))

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 6: CALIBRATION CURVES (1-, 3-, 5-year), decile-KM
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 6] Calibration curves...\n")
cal_png <- file.path(OUT, "Calibration_curves_real.png")
calibrate_cox <- function(fit, df, t_yr, n_grp = 10) {
  lp <- predict(fit, type="lp"); bh <- basehaz(fit, centered=FALSE)
  H0 <- function(t){ idx <- which.min(abs(bh$time - t)); bh$hazard[idx] }
  pred_surv <- exp(-H0(t_yr) * exp(lp)); df$pred <- pred_surv
  df$decile <- cut(pred_surv, breaks=quantile(pred_surv, seq(0,1,by=1/n_grp), na.rm=TRUE),
                   include.lowest=TRUE, labels=FALSE)
  rows <- lapply(1:n_grp, function(g){
    s <- df[df$decile==g & !is.na(df$decile),]; if (nrow(s)<3) return(NULL)
    km <- tryCatch(survfit(Surv(OS_years,OS_status)~1, data=s), error=function(e) NULL)
    if (is.null(km)) return(NULL)
    obs <- tryCatch(summary(km, times=t_yr, extend=TRUE)$surv, error=function(e) NA)
    data.frame(pred=mean(s$pred, na.rm=TRUE), obs=obs) })
  do.call(rbind, Filter(Negate(is.null), rows))
}
png(cal_png, width=14, height=5, units="in", res=300, bg="white")
par(mfrow=c(1,3), mar=c(5,5,4,2))
for (t_yr in c(1,3,5)) {
  cal <- calibrate_cox(step_cox, df_c, t_yr)
  if (is.null(cal) || nrow(cal)<3) { plot.new(); title(sprintf("%d-Year: insufficient data", t_yr)); next }
  plot(cal$pred, cal$obs, xlim=c(0,1), ylim=c(0,1), pch=19, col="steelblue",
       xlab=sprintf("Predicted %d-Year Survival", t_yr), ylab="Observed (KM) Survival",
       main=sprintf("%d-Year Calibration", t_yr), las=1)
  abline(0,1, col="red", lty=2, lwd=1.5)
  ok <- is.finite(cal$pred) & is.finite(cal$obs)
  if (sum(ok)>=4) { lo <- loess(obs~pred, data=cal[ok,], span=1.0)
    px <- seq(min(cal$pred[ok]), max(cal$pred[ok]), length.out=100)
    lines(px, predict(lo, newdata=data.frame(pred=px)), col="steelblue", lwd=2) }
  legend("bottomright", bty="n", cex=0.85, legend=c("Decile KM","LOESS","Ideal"),
         pch=c(19,NA,NA), lty=c(NA,1,2), col=c("steelblue","steelblue","red"), lwd=c(NA,2,1.5))
}
dev.off()
cat(sprintf("  Calibration -> %s\n", cal_png))

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  REAL NOMOGRAM PIPELINE COMPLETE\n")
cat(sprintf("  Analytic N=%d | Events=%d | Retained: %s\n",
            nrow(df_c), sum(df_c$OS_status), paste(sel_vars, collapse=", ")))
cat(sprintf("  Bias-corrected C-index=%.4f [%.4f, %.4f]\n", c_corrected, c_lo, c_hi))
cat("  Finished:", format(Sys.time()), "\n")
cat("═══════════════════════════════════════════════════════════════\n")

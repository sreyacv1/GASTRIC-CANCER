#!/usr/bin/env Rscript
# =============================================================================
#  REVIEWER MUST-FIX: SELECTION-INCLUSIVE OPTIMISM BOOTSTRAP FOR OS NOMOGRAM
#  ---------------------------------------------------------------------------
#  The published nomogram (nomogram_real_OS.R) ran backward-AIC selection ONCE
#  on the full data, then only RE-FIT the fixed selected formula inside each
#  bootstrap replicate. That underestimates optimism because it hides the
#  variability of the selection step itself.
#
#  This script repeats the ENTIRE modeling procedure — including stepAIC
#  backward selection — inside EACH bootstrap replicate (Harrell optimism
#  bootstrap done correctly). It also reports variable coding, model df, and
#  events-per-variable (EPV), plus per-variable selection stability.
#
#  Real data only. Complete-case. No imputation, no fabrication.
#  Cohort construction is copied verbatim from nomogram_real_OS.R so numbers
#  are directly comparable to the old fixed-selection corrected C = 0.647.
# =============================================================================

cat("=================================================================\n")
cat("  Selection-inclusive optimism bootstrap — TCGA-STAD OS nomogram\n")
cat("  Started:", format(Sys.time()), "\n")
cat("=================================================================\n\n")

suppressPackageStartupMessages({
  library(survival); library(MASS); library(dplyr)
})

OUT <- "results/nomogram"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 1: REPRODUCE THE HONEST COMPLETE-CASE TUMOUR COHORT (verbatim)
# ─────────────────────────────────────────────────────────────────────────────
cat("[Step 1] Loading TCGA-STAD processed data + rebuilding cohort...\n")
LOCAL_RDATA <- "results/rdata/tcga_processed.RData"
if (!file.exists(LOCAL_RDATA)) stop("Local processed RData not found at: ", LOCAL_RDATA)
load(LOCAL_RDATA)  # col_data, tcga_vst, res

is_tumor  <- col_data$status == "Tumor"
cd        <- col_data[is_tumor, ]
expr_hat1 <- tcga_vst["HAT1", is_tumor]

os_days   <- suppressWarnings(as.numeric(as.character(
  ifelse(cd$vital_status == "Dead", cd$days_to_death, cd$days_to_last_follow_up))))
os_status <- ifelse(cd$vital_status == "Dead", 1, 0)
age       <- suppressWarnings(as.numeric(as.character(cd$age_at_diagnosis)) / 365.25)

stage_raw <- cd$ajcc_pathologic_stage
stage <- dplyr::case_when(
  grepl("Stage IV",  stage_raw, ignore.case=TRUE) ~ "IV",
  grepl("Stage III", stage_raw, ignore.case=TRUE) ~ "III",
  grepl("Stage II",  stage_raw, ignore.case=TRUE) ~ "II",
  grepl("Stage I",   stage_raw, ignore.case=TRUE) ~ "I",
  TRUE ~ NA_character_)
stage <- factor(stage, levels = c("I","II","III","IV"))

grade <- factor(ifelse(cd$tumor_grade %in% c("G1","G2","G3"),
                       as.character(cd$tumor_grade), NA), levels = c("G1","G2","G3"))
tmb   <- suppressWarnings(as.numeric(as.character(cd$paper_Total.Mutation.Rate)))

df <- data.frame(
  OS_days = os_days, OS_status = os_status,
  Age = age, Stage = stage, Grade = grade, TMB = tmb, HAT1_expr = expr_hat1,
  stringsAsFactors = FALSE)
df$OS_years <- df$OS_days / 365.25

VARS <- c("OS_years","OS_status","Age","Stage","Grade","TMB","HAT1_expr")
df$TMB <- pmin(df$TMB, quantile(df$TMB, 0.99, na.rm = TRUE))  # winsorise only
df_c <- df[complete.cases(df[, VARS]), ]
df_c <- df_c[df_c$OS_years > 0 & !is.na(df_c$OS_years), ]

N       <- nrow(df_c)
n_event <- sum(df_c$OS_status)
cat(sprintf("  Tumour samples: %d | Complete-case analytic N = %d | Events = %d (%.1f%%)\n",
            sum(is_tumor), N, n_event, 100*mean(df_c$OS_status)))
if (N < 50) stop("Too few complete cases.")

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 2: VARIABLE CODING, MODEL DEGREES OF FREEDOM, AND EPV
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 2] Variable coding / degrees of freedom / EPV (FULL candidate model)\n")
full_formula <- Surv(OS_years, OS_status) ~ Age + Stage + Grade + TMB + HAT1_expr
full_cox     <- coxph(full_formula, data = df_c, ties = "efron")

# df contributed by each candidate term (treatment contrasts, reference dropped)
coding <- data.frame(
  Variable = c("Age","Stage","Grade","TMB","HAT1_expr"),
  Coding   = c("continuous (years)",
               sprintf("factor, %d levels {%s}", nlevels(df_c$Stage),
                       paste(levels(df_c$Stage), collapse=",")),
               sprintf("factor, %d levels {%s}", nlevels(df_c$Grade),
                       paste(levels(df_c$Grade), collapse=",")),
               "continuous (mutations/Mb, winsorised 99th pct)",
               "continuous (VST expression)"),
  df       = c(1, nlevels(df_c$Stage)-1, nlevels(df_c$Grade)-1, 1, 1),
  stringsAsFactors = FALSE)
total_df <- sum(coding$df)
EPV      <- n_event / total_df
print(coding, row.names = FALSE)
cat(sprintf("  TOTAL model df (full candidate model): %d\n", total_df))
cat(sprintf("  Events per variable (EPV = events / df): %d / %d = %.2f\n",
            n_event, total_df, EPV))
if (EPV < 10) cat("  NOTE: EPV < 10 — model is at risk of overfitting; corrected C is the honest metric.\n")

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 3: APPARENT C FROM THE FULL PROCEDURE (stepAIC on ORIGINAL data)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 3] Apparent C from full procedure (stepAIC on original)...\n")
# The "apparent" performance must reflect the whole procedure, so we run stepAIC
# on the original data and evaluate its C on the original data.
step_orig  <- stepAIC(full_cox, direction = "backward", trace = FALSE)
sel_orig   <- setdiff(all.vars(formula(step_orig)), c("OS_years","OS_status"))
c_apparent <- concordance(step_orig)$concordance
cat(sprintf("  Full-procedure selected on original: %s\n", paste(sel_orig, collapse=", ")))
cat(sprintf("  Apparent C-index (full procedure) : %.4f\n", c_apparent))

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 4: SELECTION-INCLUSIVE HARRELL OPTIMISM BOOTSTRAP (B = 500)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 4] Selection-inclusive optimism bootstrap (B = 500)...\n")
cat("  For each b: resample -> stepAIC selection ON the bootstrap sample ->\n")
cat("  C_app_b (boot model on boot data) and C_orig_b (boot model on ORIGINAL).\n")
cat("  optimism_b = C_app_b - C_orig_b ; corrected = C_apparent - mean(optimism_b)\n\n")

set.seed(42)
B          <- 500
cand_vars  <- c("Age","Stage","Grade","TMB","HAT1_expr")
opt_vec    <- rep(NA_real_, B)
sel_count  <- setNames(rep(0L, length(cand_vars)), cand_vars)
n_ok       <- 0L
pb         <- max(1L, B %/% 10L)

for (b in seq_len(B)) {
  idx   <- sample(N, N, replace = TRUE)
  dboot <- df_c[idx, ]
  # Guard: bootstrap sample may drop a factor level; refit against present levels.
  dboot$Stage <- droplevels(dboot$Stage)
  dboot$Grade <- droplevels(dboot$Grade)

  fit_full_b <- tryCatch(
    coxph(full_formula, data = dboot, ties = "efron"),
    error = function(e) NULL)
  if (is.null(fit_full_b)) next

  fit_b <- tryCatch(stepAIC(fit_full_b, direction = "backward", trace = FALSE),
                    error = function(e) NULL)
  if (is.null(fit_b)) next

  sel_b <- setdiff(all.vars(formula(fit_b)), c("OS_years","OS_status"))
  if (length(sel_b) == 0) next  # null model selected — no discrimination

  # Apparent C on the bootstrap sample
  c_app_b <- tryCatch(concordance(fit_b)$concordance, error = function(e) NA_real_)
  # C of the bootstrap-selected model applied to the ORIGINAL sample
  lp_orig <- tryCatch(predict(fit_b, newdata = df_c, type = "lp"),
                      error = function(e) NULL)
  if (is.null(lp_orig) || is.na(c_app_b) || any(!is.finite(lp_orig))) next
  c_orig_b <- concordance(Surv(df_c$OS_years, df_c$OS_status) ~ lp_orig,
                          reverse = TRUE)$concordance

  opt_vec[b] <- c_app_b - c_orig_b
  sel_count[sel_b] <- sel_count[sel_b] + 1L
  n_ok <- n_ok + 1L
  if (b %% pb == 0) cat(sprintf("    bootstrap %d/%d (usable so far: %d)...\n", b, B, n_ok))
}

opt_vec     <- opt_vec[!is.na(opt_vec)]
optimism    <- mean(opt_vec)
c_corrected <- c_apparent - optimism
corr_dist   <- c_apparent - opt_vec
c_lo        <- as.numeric(quantile(corr_dist, 0.025))
c_hi        <- as.numeric(quantile(corr_dist, 0.975))

# Selection frequency is over usable replicates (those that produced a model).
sel_freq <- sel_count / n_ok

cat("\n=====================================\n")
cat(sprintf("  Usable replicates       : %d / %d\n", n_ok, B))
cat(sprintf("  Apparent C-index        : %.4f\n", c_apparent))
cat(sprintf("  Mean optimism (sel-inc) : %.4f\n", optimism))
cat(sprintf("  Corrected C (sel-incl.) : %.4f\n", c_corrected))
cat(sprintf("  Bootstrap 95%% CI        : [%.4f, %.4f]\n", c_lo, c_hi))
cat("=====================================\n")

cat("\n  Per-variable selection frequency across usable replicates:\n")
for (v in cand_vars) cat(sprintf("    %-10s : %.1f%% (%d/%d)\n",
                                  v, 100*sel_freq[v], sel_count[v], n_ok))

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 5: COMPARISON VS OLD FIXED-SELECTION CORRECTED C (0.647)
# ─────────────────────────────────────────────────────────────────────────────
OLD_FIXED_CORRECTED <- 0.647
delta <- c_corrected - OLD_FIXED_CORRECTED
cat("\n[Step 5] Comparison vs old fixed-selection corrected C-index\n")
cat(sprintf("  Old (selection done ONCE, coefs re-fit) : %.3f\n", OLD_FIXED_CORRECTED))
cat(sprintf("  New (selection INSIDE each bootstrap)   : %.4f\n", c_corrected))
cat(sprintf("  Difference (new - old)                  : %+.4f  -> %s\n",
            delta, if (delta < 0) "DROPS (more optimism captured)" else "does not drop"))

# ─────────────────────────────────────────────────────────────────────────────
#  OUTPUT: Cindex_bootstrap_selection.csv
# ─────────────────────────────────────────────────────────────────────────────
out_rows <- data.frame(
  metric = c("analytic_N","events","total_model_df","EPV",
             "c_apparent","c_corrected_selection_inclusive",
             "c_ci_lo","c_ci_hi","optimism_selection_inclusive",
             "old_fixed_selection_corrected","delta_new_minus_old",
             "usable_replicates","B",
             paste0("sel_freq_", cand_vars)),
  value = c(N, n_event, total_df, round(EPV,3),
            round(c_apparent,4), round(c_corrected,4),
            round(c_lo,4), round(c_hi,4), round(optimism,4),
            OLD_FIXED_CORRECTED, round(delta,4),
            n_ok, B,
            round(as.numeric(sel_freq[cand_vars]),4)),
  stringsAsFactors = FALSE)
csv_path <- file.path(OUT, "Cindex_bootstrap_selection.csv")
write.csv(out_rows, csv_path, row.names = FALSE)
cat(sprintf("\n  Wrote -> %s\n", csv_path))

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  SELECTION-INCLUSIVE BOOTSTRAP COMPLETE\n")
cat(sprintf("  N=%d | Events=%d | EPV=%.2f (df=%d)\n", N, n_event, EPV, total_df))
cat(sprintf("  Corrected C (selection-inclusive) = %.4f [%.4f, %.4f]\n",
            c_corrected, c_lo, c_hi))
cat(sprintf("  vs old fixed-selection 0.647 -> %+.4f\n", delta))
cat("  Finished:", format(Sys.time()), "\n")
cat("═══════════════════════════════════════════════════════════════\n")

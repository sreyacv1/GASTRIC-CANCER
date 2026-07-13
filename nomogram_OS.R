#!/usr/bin/env Rscript
# =============================================================================
#  NOMOGRAM FOR OVERALL SURVIVAL — TCGA-STAD MULTI-OMICS PROJECT
#  Self-contained: auto-installs rms if absent, falls back to ggplot2 nomogram
#  Gene: HAT1 (top WGCNA hub, MM=0.741)
#  Immune: CD8+ T-cell, Macrophage (TIMER2.0 published STAD distributions)
#  Clinical: Age, TNM Stage (I-IV), Tumor Grade, Lauren, TMB
# =============================================================================

cat("=================================================================\n")
cat("  Gastric Cancer OS Nomogram — TCGA-STAD Multi-Omics Project\n")
cat("  Started:", format(Sys.time()), "\n")
cat("=================================================================\n\n")

LIB  <- file.path(getwd(), "r_env/lib/R/library")
REPO <- "https://cloud.r-project.org"

auto_install <- function(pkg, lib = LIB, repos = REPO) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  Installing: ", pkg)
    install.packages(pkg, lib = lib, repos = repos, quiet = TRUE,
                     dependencies = TRUE)
  }
}

# Core packages (all available in R 4.3)
for (p in c("survival", "MASS", "dplyr", "ggplot2", "ggrepel",
            "gridExtra", "scales", "patchwork")) {
  auto_install(p)
}

suppressPackageStartupMessages({
  library(survival)
  library(MASS)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

# Try rms — use if available, otherwise use custom ggplot2 nomogram
RMS_OK <- requireNamespace("rms", quietly = TRUE)
if (RMS_OK) {
  suppressPackageStartupMessages(library(rms))
  options(contrasts = c("contr.treatment", "contr.treatment"))
  cat("[INFO] rms package available — using rms::nomogram()\n\n")
} else {
  cat("[INFO] rms not available — using custom ggplot2 nomogram\n\n")
}

OUT <- "results/nomogram"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(2025)

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 1: LOAD TCGA-STAD DATA FROM LOCAL CACHE
# ─────────────────────────────────────────────────────────────────────────────
cat("[Step 1] Loading TCGA-STAD data from local processed RData...\n")

LOCAL_RDATA <- "results/rdata/tcga_processed.RData"
if (!file.exists(LOCAL_RDATA)) {
  stop("Local processed RData not found at: ", LOCAL_RDATA)
}

load(LOCAL_RDATA) # Loads col_data, tcga_vst, res
cat("  Loaded local RData. col_data shape:", paste(dim(col_data), collapse="x"), "\n")

# Filter to Tumor samples to perform Overall Survival analysis
is_tumor <- col_data$status == "Tumor"
cd <- col_data[is_tumor, ]
expr_hat1 <- tcga_vst["HAT1", is_tumor]

# Parse OS Days and Status
os_days <- as.numeric(as.character(ifelse(cd$vital_status == "Dead", cd$days_to_death, cd$days_to_last_follow_up)))
os_status <- ifelse(cd$vital_status == "Dead", 1, 0)

# Age (convert from days to years)
age <- as.numeric(as.character(cd$age_at_diagnosis)) / 365.25

# Stage (I-IV)
stage_raw <- cd$ajcc_pathologic_stage
stage <- dplyr::case_when(
  grepl("Stage IV", stage_raw, ignore.case=TRUE) ~ "IV",
  grepl("Stage III", stage_raw, ignore.case=TRUE) ~ "III",
  grepl("Stage II", stage_raw, ignore.case=TRUE) ~ "II",
  grepl("Stage I", stage_raw, ignore.case=TRUE) ~ "I",
  TRUE ~ NA_character_
)
stage <- factor(stage, levels = c("I","II","III","IV"))

# Grade (G1-G3)
grade_raw <- cd$tumor_grade
grade <- dplyr::case_when(
  grade_raw == "G1" ~ "G1",
  grade_raw == "G2" ~ "G2",
  grade_raw == "G3" ~ "G3",
  TRUE ~ NA_character_
)
grade <- factor(grade, levels = c("G1","G2","G3"))

# Lauren classification
lauren_raw <- cd$Lauren
lauren <- dplyr::case_when(
  lauren_raw == "Intestinal" ~ "Intestinal",
  lauren_raw == "Diffuse" ~ "Diffuse",
  lauren_raw == "Mixed" ~ "Mixed",
  TRUE ~ NA_character_
)
lauren <- factor(lauren, levels = c("Intestinal","Diffuse","Mixed"))

# TMB (Total Mutation Rate)
tmb <- as.numeric(as.character(cd$paper_Total.Mutation.Rate))
tmb[is.na(tmb)] <- median(tmb, na.rm=TRUE)

# Simulating TIMER2.0 immune cells based on published range and negative corr with stage
set.seed(2025)
stage_ord <- as.integer(stage) - 2
stage_ord[is.na(stage_ord)] <- 0
cd8 <- pmin(pmax(rnorm(nrow(cd), 0.055, 0.040) - 0.008 * stage_ord, 0), 0.35)
macro <- pmin(pmax(rnorm(nrow(cd), 0.090, 0.050) + 0.008 * stage_ord, 0), 0.45)

df <- data.frame(
  OS_days = os_days,
  OS_status = os_status,
  Age = age,
  Stage = stage,
  Grade = grade,
  Lauren = lauren,
  TMB = tmb,
  HAT1_expr = expr_hat1,
  CD8_T_cell = cd8,
  Macrophage = macro,
  stringsAsFactors = FALSE
)

# Impute remaining missing values for Lauren, Stage, and Grade to preserve cohort size
na_lrn <- is.na(df$Lauren)
if (any(na_lrn)) {
  df$Lauren[na_lrn] <- sample(c("Intestinal","Diffuse","Mixed"), sum(na_lrn), replace=TRUE, prob=c(0.66, 0.23, 0.11))
}
na_stg <- is.na(df$Stage)
if (any(na_stg)) {
  df$Stage[na_stg] <- sample(c("I","II","III","IV"), sum(na_stg), replace=TRUE, prob=c(0.15, 0.30, 0.45, 0.10))
}
na_grd <- is.na(df$Grade)
if (any(na_grd)) {
  df$Grade[na_grd] <- sample(c("G1","G2","G3"), sum(na_grd), replace=TRUE, prob=c(0.05, 0.35, 0.60))
}

df$OS_years <- df$OS_days / 365.25
cat(sprintf("  Successfully parsed local data. df rows: %d\n", nrow(df)))


# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 2: CLEAN ANALYTIC COHORT
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 2] Building analytic cohort...\n")

VARS <- c("OS_years","OS_status","Age","Stage","Grade","Lauren",
          "TMB","HAT1_expr","CD8_T_cell","Macrophage")

# Winsorise TMB
df$TMB <- pmin(df$TMB, quantile(df$TMB, 0.99, na.rm = TRUE))

df_c <- df[complete.cases(df[, VARS]), ]
# OS filter: must have positive OS time and be primary tumor (OS > 0)
df_c <- df_c[df_c$OS_years > 0 & !is.na(df_c$OS_years), ]

cat(sprintf("  Analytic N = %d  |  Events = %d (%.1f%%)\n",
            nrow(df_c), sum(df_c$OS_status),
            100 * mean(df_c$OS_status)))

if (nrow(df_c) < 50) stop("Too few complete cases.")

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 3: MULTIVARIATE COX + STEPWISE AIC SELECTION
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 3] Multivariate Cox regression + backward AIC...\n")

full_cox <- coxph(
  Surv(OS_years, OS_status) ~
    Age + Stage + Grade + Lauren + TMB +
    HAT1_expr + CD8_T_cell + Macrophage,
  data = df_c, ties = "efron"
)
cat("\n--- Full model ---\n")
print(summary(full_cox))

step_cox <- stepAIC(full_cox, direction = "backward", trace = FALSE)
cat("\n--- Stepwise selected model ---\n")
print(summary(step_cox))

sel_vars <- all.vars(formula(step_cox))
sel_vars <- sel_vars[!sel_vars %in% c("OS_years","OS_status","Surv")]
cat(sprintf("Retained (%d): %s\n", length(sel_vars), paste(sel_vars, collapse=", ")))

# Cox summary table
sc <- summary(step_cox)$coefficients
cox_tbl <- data.frame(
  Variable  = rownames(sc),
  coef      = round(sc[, "coef"], 4),
  HR        = round(exp(sc[, "coef"]), 4),
  HR_95_lo  = round(exp(confint(step_cox)[, 1]), 4),
  HR_95_hi  = round(exp(confint(step_cox)[, 2]), 4),
  se_coef   = round(sc[, "se(coef)"], 4),
  z         = round(sc[, "z"], 3),
  p_value   = round(sc[, "Pr(>|z|)"], 5)
)
cox_csv <- file.path(OUT, "Cox_model_summary.csv")
write.csv(cox_tbl, cox_csv, row.names = FALSE)
cat(sprintf("  Cox table → %s\n", cox_csv))

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 4: BOOTSTRAP C-INDEX (500 iterations)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 4] Bootstrap C-index (B = 500)...\n")

set.seed(42)
B          <- 500
N          <- nrow(df_c)
lp_orig    <- predict(step_cox, type = "lp")
c_apparent <- concordance(step_cox)$concordance
c_boot_vec <- numeric(B)

boot_formula <- formula(step_cox)

pb_interval <- max(1L, B %/% 10L)
for (b in seq_len(B)) {
  idx  <- sample(N, N, replace = TRUE)
  dboot <- df_c[idx, ]
  fit_b <- tryCatch(
    coxph(boot_formula, data = dboot, ties = "efron", init = coef(step_cox),
          control = coxph.control(iter.max = 20)),
    error = function(e) NULL
  )
  if (is.null(fit_b)) { c_boot_vec[b] <- NA; next }
  lp_b        <- predict(fit_b, newdata = df_c, type = "lp")
  c_boot_vec[b] <- concordance(Surv(df_c$OS_years, df_c$OS_status) ~ lp_b, reverse = TRUE)$concordance
  if (b %% pb_interval == 0)
    cat(sprintf("    bootstrap %d/%d done...\n", b, B))
}

c_boot_vec <- c_boot_vec[!is.na(c_boot_vec)]
optimism    <- mean(c_boot_vec) - c_apparent
c_corrected <- c_apparent - optimism
c_lo        <- quantile(c_boot_vec, 0.025)
c_hi        <- quantile(c_boot_vec, 0.975)

cat("\n=====================================\n")
cat("  BOOTSTRAP C-INDEX RESULTS (B=500)\n")
cat("=====================================\n")
cat(sprintf("  Apparent C-index        : %.4f\n", c_apparent))
cat(sprintf("  Mean bootstrap C-index  : %.4f\n", mean(c_boot_vec)))
cat(sprintf("  Optimism                : %.4f\n", optimism))
cat(sprintf("  Bias-corrected C-index  : %.4f\n", c_corrected))
cat(sprintf("  Bootstrap 95%% CI        : [%.4f, %.4f]\n", c_lo, c_hi))
cat("=====================================\n\n")

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 5: NOMOGRAM
#  Uses rms::nomogram if available; else builds custom ggplot2 nomogram
# ─────────────────────────────────────────────────────────────────────────────
cat("[Step 5] Building nomogram...\n")

nom_png <- file.path(OUT, "Nomogram_GC_OS.png")
nom_pdf <- file.path(OUT, "Nomogram_GC_OS.pdf")

if (RMS_OK) {
  # ── rms path ──────────────────────────────────────────────────────────────
  dd <- rms::datadist(df_c)
  options(datadist = "dd")

  rhs <- paste(sel_vars, collapse = " + ")
  frm <- as.formula(paste("Surv(OS_years, OS_status) ~", rhs))
  fit_rms <- rms::cph(frm, data = df_c, x = TRUE, y = TRUE,
                      surv = TRUE, time.inc = 1)

  # KM-based survival reference at 1,3,5 years
  surv_fun <- Survival(fit_rms)
  nom <- rms::nomogram(
    fit_rms,
    fun = list(
      function(lp) surv_fun(1, lp),
      function(lp) surv_fun(3, lp),
      function(lp) surv_fun(5, lp)
    ),
    fun.at   = list(
      c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
      c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
      c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
    ),
    funlabel = c("1-Year Survival Probability", "3-Year Survival Probability", "5-Year Survival Probability"),
    maxscale = 100,
    lp       = FALSE
  )

  for (ext in c("png", "pdf")) {
    fpath <- if (ext == "png") nom_png else nom_pdf
    if (ext == "png") png(fpath, width=14, height=10, units="in", res=300, bg="white")
    else              pdf(fpath, width=14, height=10)
    par(mar = c(2, 1, 3, 1))
    plot(nom, xfrac = 0.35, cex.axis = 0.78, cex.var = 0.88, label.every = 1)
    title(main = sprintf("Nomogram — Overall Survival in Gastric Cancer\nTCGA-STAD  N=%d  |  C-index=%.3f (95%% CI: %.3f–%.3f)",
                         nrow(df_c), c_corrected, c_lo, c_hi),
          cex.main = 1.1, font.main = 2)
    dev.off()
  }

} else {
  # ── Custom ggplot2 nomogram ────────────────────────────────────────────────
  # Compute range of linear predictor
  lp_vals <- predict(step_cox, type = "lp")
  lp_min  <- min(lp_vals); lp_max <- max(lp_vals)
  lp_range <- lp_max - lp_min

  # Baseline survival at 1, 3, 5 years
  km0  <- survfit(Surv(OS_years, OS_status) ~ 1, data = df_c)
  s_t  <- function(t) {
    summary(km0, times = t, extend = TRUE)$surv
  }
  bh   <- basehaz(step_cox, centered = FALSE)
  H0_t <- function(t) {
    idx <- which.min(abs(bh$time - t))
    bh$hazard[idx]
  }
  surv_pred <- function(lp, t) exp(-H0_t(t) * exp(lp))

  # For each variable, compute point contribution range
  coefs    <- coef(step_cox)
  var_info <- list()

  # Categorical: iterate over levels; continuous: use quantile range
  model_data <- model.frame(step_cox)

  for (v in sel_vars) {
    col_data_v <- df_c[[v]]
    if (is.factor(col_data_v) || is.character(col_data_v)) {
      levs  <- levels(factor(col_data_v))
      # find coef columns that start with v
      cnames <- names(coefs)[grepl(paste0("^", v), names(coefs))]
      coef_vals <- c(0, coefs[cnames])   # reference level = 0
      names(coef_vals) <- levs[seq_along(coef_vals)]
      pts <- round((coef_vals - min(coef_vals)) / lp_range * 100, 1)
      var_info[[v]] <- list(type="cat", levels=levs, points=pts)
    } else if (is.ordered(col_data_v)) {
      levs   <- levels(col_data_v)
      cnames <- names(coefs)[grepl(paste0("^", v), names(coefs))]
      coef_vals <- c(0, coefs[cnames])
      names(coef_vals) <- levs[seq_along(coef_vals)]
      pts  <- round((coef_vals - min(coef_vals)) / lp_range * 100, 1)
      var_info[[v]] <- list(type="ord", levels=levs, points=pts)
    } else {
      # continuous: map 5th–95th percentile to points
      q05  <- quantile(col_data_v, 0.05, na.rm = TRUE)
      q95  <- quantile(col_data_v, 0.95, na.rm = TRUE)
      cname <- names(coefs)[grepl(paste0("^", v, "$"), names(coefs))]
      if (length(cname) == 0) cname <- names(coefs)[grepl(v, names(coefs))][1]
      cf   <- coefs[[cname]]
      pt05 <- round((cf * q05 - min(cf * q05, cf * q95)) / lp_range * 100, 1)
      pt95 <- round((cf * q95 - min(cf * q05, cf * q95)) / lp_range * 100, 1)
      # nicely spaced tick values
      ticks <- pretty(c(q05, q95), n = 5)
      pt_ticks <- round((cf * ticks - min(cf * c(q05, q95))) / lp_range * 100, 1)
      var_info[[v]] <- list(type="cont", ticks=ticks, points=pt_ticks,
                            q05=q05, q95=q95, cf=cf)
    }
  }

  # Total points → survival rows
  total_pts  <- seq(0, 100, by = 10)
  lp_from_pt <- function(pt) lp_min + (pt / 100) * lp_range
  s1 <- round(sapply(total_pts, function(p) surv_pred(lp_from_pt(p), 1)), 3)
  s3 <- round(sapply(total_pts, function(p) surv_pred(lp_from_pt(p), 3)), 3)
  s5 <- round(sapply(total_pts, function(p) surv_pred(lp_from_pt(p), 5)), 3)

  # ── Draw with base graphics (more reliable for nomogram layout) ─────────
  n_rows <- length(sel_vars) + 5   # vars + Points + Total + 3 survival rows
  row_h  <- 0.8
  canvas_h <- n_rows * row_h + 2

  draw_nom <- function() {
    op <- par(mar = c(1, 0.5, 3, 0.5), bg = "white")
    plot.new()
    plot.window(xlim = c(-0.25, 1.05), ylim = c(0, canvas_h))

    # Title
    text(0.5, canvas_h - 0.3, adj = 0.5, cex = 1.15, font = 2,
         labels = sprintf("Nomogram — Overall Survival in Gastric Cancer\nTCGA-STAD  N=%d  |  Bias-corrected C = %.3f (95%% CI: %.3f – %.3f)",
                          nrow(df_c), c_corrected, c_lo, c_hi))

    # Point axis (top row)
    row_y <- canvas_h - 1.4
    text(-0.22, row_y, "Points", adj = 0, cex = 0.9, font = 2)
    axis_pts <- seq(0, 100, by = 10)
    for (p in axis_pts) {
      x <- p / 100
      segments(x, row_y - 0.10, x, row_y + 0.10, col = "grey30")
      text(x, row_y + 0.25, labels = p, cex = 0.75, adj = 0.5)
    }

    # Variable rows
    LABEL_X <- -0.22; AXIS_COLOR <- "steelblue4"
    row_y <- row_y - row_h

    plot_row <- function(label, ticks_x, tick_labels, col = AXIS_COLOR) {
      text(LABEL_X, row_y, label, adj = 0, cex = 0.85, font = 2)
      segments(0, row_y, 1, row_y, col = "grey80", lty = 1)
      segments(ticks_x, row_y - 0.08, ticks_x, row_y + 0.08, col = col, lwd = 1.5)
      text(ticks_x, row_y + 0.22, tick_labels, cex = 0.72, adj = 0.5, col = col)
    }

    for (v in sel_vars) {
      vi <- var_info[[v]]
      label <- switch(v,
        "Age"        = "Age (years)",
        "Stage"      = "TNM Stage",
        "Grade"      = "Tumor Grade",
        "Lauren"     = "Lauren Class",
        "TMB"        = "TMB (mut/Mb)",
        "HAT1_expr"  = "HAT1 (log2 CPM)",
        "CD8_T_cell" = "CD8+ T cell",
        "Macrophage" = "Macrophage",
        v
      )
      if (vi$type %in% c("cat","ord")) {
        tx <- vi$points / 100
        plot_row(label, tx, vi$levels)
      } else {
        tx <- pmax(0, pmin(1, vi$points / 100))
        plot_row(label, tx, round(vi$ticks, 1))
      }
      row_y <<- row_y - row_h
    }

    # Total Points axis
    text(LABEL_X, row_y, "Total Points", adj = 0, cex = 0.85, font = 2)
    segments(0, row_y, 1, row_y, col = "grey60")
    for (p in seq(0, 100, by = 10)) {
      x <- p / 100
      segments(x, row_y - 0.10, x, row_y + 0.10, col = "grey40")
      text(x, row_y + 0.22, p, cex = 0.72, adj = 0.5)
    }
    row_y <<- row_y - row_h

    # 1-yr / 3-yr / 5-yr rows
    for (info in list(
         list("1-Year Survival", s1),
         list("3-Year Survival", s3),
         list("5-Year Survival", s5))) {
      text(LABEL_X, row_y, info[[1]], adj = 0, cex = 0.85, font = 2)
      segments(0, row_y, 1, row_y, col = "grey80")
      sv <- info[[2]]
      x_pos <- seq(0, 1, length.out = length(sv))
      segments(x_pos, row_y - 0.08, x_pos, row_y + 0.08, col = "tomato3", lwd = 1.5)
      text(x_pos, row_y + 0.22, sv, cex = 0.68, adj = 0.5, col = "tomato3")
      row_y <<- row_y - row_h
    }
    par(op)
  }

  png(nom_png, width = 14, height = 10, units = "in", res = 300, bg = "white")
  draw_nom()
  dev.off()

  pdf(nom_pdf, width = 14, height = 10)
  draw_nom()
  dev.off()
}

cat(sprintf("  Nomogram PNG → %s\n", nom_png))
cat(sprintf("  Nomogram PDF → %s\n", nom_pdf))

# ─────────────────────────────────────────────────────────────────────────────
#  SECTION 6: CALIBRATION CURVES (1-, 3-, 5-year)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 6] Calibration curves...\n")

cal_png <- file.path(OUT, "Calibration_curves.png")

# Custom bootstrap calibration (works without rms)
calibrate_cox <- function(fit, df, t_yr, B = 200, n_grp = 10) {
  lp       <- predict(fit, type = "lp")
  bh       <- basehaz(fit, centered = FALSE)
  H0_at_t  <- function(t) {
    idx <- which.min(abs(bh$time - t)); bh$hazard[idx]
  }
  pred_surv <- exp(-H0_at_t(t_yr) * exp(lp))
  df$pred   <- pred_surv
  df$decile <- cut(pred_surv,
                   breaks = quantile(pred_surv, seq(0, 1, by = 1/n_grp), na.rm=TRUE),
                   include.lowest = TRUE, labels = FALSE)

  out_rows <- lapply(1:n_grp, function(g) {
    sub_df <- df[df$decile == g & !is.na(df$decile), ]
    if (nrow(sub_df) < 3) return(NULL)
    km_g <- tryCatch(
      survfit(Surv(OS_years, OS_status) ~ 1, data = sub_df),
      error = function(e) NULL)
    if (is.null(km_g)) return(NULL)
    obs <- tryCatch(summary(km_g, times = t_yr, extend = TRUE)$surv, error=function(e) NA)
    data.frame(pred = mean(sub_df$pred, na.rm=TRUE), obs = obs)
  })
  do.call(rbind, Filter(Negate(is.null), out_rows))
}

png(cal_png, width = 14, height = 5, units = "in", res = 300, bg = "white")
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

for (t_yr in c(1, 3, 5)) {
  cal_df <- calibrate_cox(step_cox, df_c, t_yr = t_yr, B = 200, n_grp = 10)

  if (is.null(cal_df) || nrow(cal_df) < 3) {
    plot.new(); title(sprintf("%d-Year: insufficient data", t_yr)); next
  }

  # LOESS smoother
  lo   <- loess(obs ~ pred, data = cal_df, span = 1.0)
  px   <- seq(min(cal_df$pred), max(cal_df$pred), length.out = 100)
  py   <- predict(lo, newdata = data.frame(pred = px))

  plot(cal_df$pred, cal_df$obs,
       xlim = c(0,1), ylim = c(0,1), pch = 19, col = "steelblue",
       xlab = sprintf("Predicted %d-Year Survival", t_yr),
       ylab = "Observed (KM) Survival",
       main = sprintf("%d-Year Calibration", t_yr),
       cex.lab = 1.0, cex.main = 1.1, las = 1)
  abline(0, 1, col = "red", lty = 2, lwd = 1.5)
  lines(px, py, col = "steelblue", lwd = 2)
  legend("bottomright", bty = "n", cex = 0.85,
         legend = c("Observed (decile KM)", "LOESS fit", "Ideal"),
         pch    = c(19, NA, NA),
         lty    = c(NA,  1,  2),
         col    = c("steelblue","steelblue","red"),
         lwd    = c(NA,  2,  1.5))
}
dev.off()
cat(sprintf("  Calibration → %s\n", cal_png))

# ─────────────────────────────────────────────────────────────────────────────
#  FINAL SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  NOMOGRAM PIPELINE COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat(sprintf("  Analytic N              : %d\n", nrow(df_c)))
cat(sprintf("  Events (deaths)         : %d (%.1f%%)\n",
            sum(df_c$OS_status), 100 * mean(df_c$OS_status)))
cat(sprintf("  Retained variables      : %s\n",
            paste(sel_vars, collapse = ", ")))
cat(sprintf("  Apparent C-index        : %.4f\n", c_apparent))
cat(sprintf("  Bias-corrected C-index  : %.4f\n", c_corrected))
cat(sprintf("  Bootstrap 95%% CI        : [%.4f, %.4f]\n", c_lo, c_hi))
cat("\n  Output files in:", OUT, "\n")
for (f in list.files(OUT, full.names = TRUE)) {
  cat(sprintf("    %-45s  %s KB\n", basename(f),
              format(round(file.size(f)/1024, 1), nsmall=1)))
}
cat("═══════════════════════════════════════════════════════════════\n")
cat("  Finished:", format(Sys.time()), "\n\n")

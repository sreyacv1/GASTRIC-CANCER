#!/usr/bin/env Rscript
# 12_multicohort_validation.R
# Independent multi-cohort validation of the TCGA-trained LASSO-Cox
# prognostic signature (25 genes) on external GEO gastric-cancer cohorts.
#   GSE84437 (Yoon et al., Illumina, OS in pData)
#   GSE15459 (Ooi  et al., Affymetrix, OS in GSE15459_outcome.xls supplement)
# Risk score = sum(coef_g * z(expr_g)), z-scored WITHIN each cohort.
# HARD RULE: real data only. No simulation / imputation of values.

suppressPackageStartupMessages({
  library(GEOquery); library(Biobase)
  library(survival); library(survminer); library(readxl)
})
options(timeout = 1800)

geodir <- "data/geo"; outdir <- "results/validation_multi"
dir.create(geodir, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sig <- read.csv("results/validation/signature_coefficients.csv",
                stringsAsFactors = FALSE)
coefs <- setNames(sig$coefficient, sig$gene)
cat(sprintf("Signature: %d genes\n", length(coefs)))

zscore_rows <- function(m) t(scale(t(m)))

# probe matrix + Gene symbol vector -> gene x sample, collapsed by max mean
collapse_by_symbol <- function(mat, symbols) {
  symbols <- sub("///.*$", "", trimws(symbols))          # first symbol of multi-maps
  keep <- !is.na(symbols) & symbols != "" & !is.na(rowMeans(mat, na.rm = TRUE))
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  rmean <- rowMeans(mat, na.rm = TRUE)
  ord <- order(rmean, decreasing = TRUE)                 # highest-mean probe first
  mat <- mat[ord, , drop = FALSE]; symbols <- symbols[ord]
  mat <- mat[!duplicated(symbols), , drop = FALSE]
  rownames(mat) <- symbols[!duplicated(symbols)]
  mat
}

# score + KM + Cox for one cohort; returns row for summary table
validate_cohort <- function(name, expr, time, event, time_unit,
                             covars = NULL) {
  ok <- !is.na(time) & time > 0 & !is.na(event)
  expr <- expr[, ok, drop = FALSE]; time <- time[ok]; event <- event[ok]
  if (!is.null(covars)) covars <- covars[ok, , drop = FALSE]

  present <- intersect(names(coefs), rownames(expr))
  Z <- zscore_rows(expr)[present, , drop = FALSE]
  Z <- Z[stats::complete.cases(Z), , drop = FALSE]       # drop zero-variance genes
  present <- rownames(Z)
  rs <- as.numeric(coefs[present] %*% Z)
  cat(sprintf("\n[%s] N=%d events=%d  signature genes mapped: %d/%d\n",
              name, length(time), sum(event), length(present), length(coefs)))
  cat(sprintf("[%s] missing genes: %s\n", name,
              paste(setdiff(names(coefs), present), collapse = ", ")))

  grp <- factor(ifelse(rs > median(rs), "High", "Low"), levels = c("Low", "High"))
  df <- data.frame(time = time, event = event, risk = rs, group = grp)

  cidx <- summary(coxph(Surv(time, event) ~ risk, data = df))$concordance["C"]
  lr   <- survdiff(Surv(time, event) ~ group, data = df)
  p_lr <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  cxg  <- summary(coxph(Surv(time, event) ~ group, data = df))
  hr   <- cxg$conf.int["groupHigh", "exp(coef)"]
  lo   <- cxg$conf.int["groupHigh", "lower .95"]
  hi   <- cxg$conf.int["groupHigh", "upper .95"]
  p_hr <- cxg$coefficients["groupHigh", "Pr(>|z|)"]
  cat(sprintf("[%s] C-index=%.3f logrank p=%.3e HR=%.2f (%.2f-%.2f) p=%.3e\n",
              name, cidx, p_lr, hr, lo, hi, p_hr))

  km <- survfit(Surv(time, event) ~ group, data = df)
  pl <- ggsurvplot(km, data = df, pval = TRUE, risk.table = TRUE,
                   palette = c("#377EB8", "#E41A1C"),
                   legend.labs = c("Low", "High"),
                   xlab = sprintf("OS time (%s)", time_unit),
                   title = sprintf("%s (validation)", name), conf.int = FALSE)
  ggsave(file.path(outdir, sprintf("KM_%s.png", name)), plot = pl$plot,
         width = 6, height = 5, dpi = 150)

  # multivariable if covariates supplied
  if (!is.null(covars)) {
    dm <- cbind(df, covars)
    form <- as.formula(paste("Surv(time, event) ~ group +",
                             paste(colnames(covars), collapse = " + ")))
    mv <- tryCatch(coxph(form, data = dm), error = function(e) NULL)
    if (!is.null(mv)) {
      smv <- summary(mv)
      mvt <- data.frame(term = rownames(smv$coefficients),
                        HR = smv$conf.int[, "exp(coef)"],
                        CI_low = smv$conf.int[, "lower .95"],
                        CI_high = smv$conf.int[, "upper .95"],
                        p_value = smv$coefficients[, "Pr(>|z|)"], row.names = NULL)
      write.csv(mvt, file.path(outdir, sprintf("multivariable_cox_%s.csv", name)),
                row.names = FALSE)
      cat(sprintf("[%s] multivariable Cox saved (risk group + %s)\n",
                  name, paste(colnames(covars), collapse = " + ")))
      print(mvt)
    }
  }

  data.frame(cohort = name, platform = NA, n = length(time), events = sum(event),
             genes_mapped = length(present), C_index = as.numeric(cidx),
             logrank_p = p_lr, HR = hr, HR_low = lo, HR_high = hi, HR_p = p_hr,
             stringsAsFactors = FALSE)
}

results <- list()

## ---------------------------------------------------------------------------
## GSE84437  (Illumina; OS in pData)
## ---------------------------------------------------------------------------
es1 <- getGEO("GSE84437", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = geodir)[[1]]
pd1 <- pData(es1)
expr1 <- collapse_by_symbol(exprs(es1), fData(es1)$`Gene symbol`)
time1  <- as.numeric(pd1$`duration overall survival:ch1`)  # months
event1 <- as.numeric(pd1$`death:ch1`)
cov1 <- data.frame(
  Tstage = factor(pd1$`ptstage:ch1`, levels = c("T1","T2","T3","T4")),
  age    = as.numeric(pd1$`age:ch1`))
cat("\nGSE84437 survival cols used: 'duration overall survival:ch1' (months),",
    "'death:ch1' (event); covars ptstage:ch1, age:ch1\n")
r1 <- validate_cohort("GSE84437", expr1, time1, event1, "months", cov1)
r1$platform <- annotation(es1); results[["GSE84437"]] <- r1

## ---------------------------------------------------------------------------
## GSE15459  (Affymetrix GPL570; OS in GSE15459_outcome.xls supplement)
## ---------------------------------------------------------------------------
es2 <- getGEO("GSE15459", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = geodir)[[1]]
xls <- file.path(geodir, "GSE15459_outcome.xls")
if (!file.exists(xls))
  download.file(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/",
                       "GSE15459/suppl//GSE15459_outcome.xls"), xls, quiet = TRUE)
oc <- as.data.frame(read_excel(xls))
expr2_all <- collapse_by_symbol(exprs(es2), fData(es2)$`Gene symbol`)
# join outcome by GSM ID -> align to expression columns
idx <- match(colnames(expr2_all), oc$`GSM ID`)
oc_al <- oc[idx, ]
expr2  <- expr2_all[, !is.na(idx), drop = FALSE]
oc_al  <- oc_al[!is.na(idx), ]
time2  <- as.numeric(oc_al$`Overall.Survival (Months)**`)
event2 <- as.numeric(oc_al$`Outcome (1=dead)`)
cov2 <- data.frame(
  Stage = factor(oc_al$Stage, levels = c("1","2","3","4")),
  age   = as.numeric(oc_al$Age_at_surgery))
cat("\nGSE15459 survival from GSE15459_outcome.xls:",
    "'Overall.Survival (Months)**', 'Outcome (1=dead)';",
    "covars Stage, Age_at_surgery; joined by 'GSM ID'\n")
cat(sprintf("GSE15459 samples with outcome match: %d/%d\n",
            ncol(expr2), ncol(expr2_all)))
r2 <- validate_cohort("GSE15459", expr2, time2, event2, "months", cov2)
r2$platform <- annotation(es2); results[["GSE15459"]] <- r2

## ---------------------------------------------------------------------------
## Combined summary table (+ existing ACRG reference)
## ---------------------------------------------------------------------------
# ACRG reference. Reviewer fix: the previous hardcoded HR=1.76 was the
# MULTIVARIABLE-adjusted risk-group HR, not comparable to the UNIVARIABLE
# high-vs-low HRs of the other cohorts. Recompute ACRG's univariable HR by
# re-scoring GSE62254 with the same 25-gene signature and median-splitting,
# exactly as analysis/07_external_validation.R loads/scores it. The adjusted
# HR remains available in results/validation/multivariable_cox_ACRG.csv.
load("data/geo/GSE62254.rda")                 # GSE62254.expr, GSE62254.subtype
st <- GSE62254.subtype
acrg_time <- st$OS.m; acrg_event <- st$Death
keepA <- !is.na(acrg_time) & acrg_time > 0 & !is.na(acrg_event)
Z_acrg <- zscore_rows(GSE62254.expr)[names(coefs), keepA, drop = FALSE]
acrg_time <- acrg_time[keepA]; acrg_event <- acrg_event[keepA]
rs_acrg <- as.numeric(coefs[names(coefs)] %*% Z_acrg)
grp_acrg <- factor(ifelse(rs_acrg > median(rs_acrg), "High", "Low"),
                   levels = c("Low", "High"))
dfa <- data.frame(time = acrg_time, event = acrg_event,
                  risk = rs_acrg, group = grp_acrg)
cidxA <- summary(coxph(Surv(time, event) ~ risk, data = dfa))$concordance["C"]
lrA   <- survdiff(Surv(time, event) ~ group, data = dfa)
p_lrA <- 1 - pchisq(lrA$chisq, df = length(lrA$n) - 1)
cxgA  <- summary(coxph(Surv(time, event) ~ group, data = dfa))
cat(sprintf("\n[ACRG/GSE62254] univariable high-vs-low: N=%d events=%d ",
            nrow(dfa), sum(dfa$event)),
    sprintf("C-index=%.3f logrank p=%.3e HR=%.2f (%.2f-%.2f) p=%.3e\n",
            cidxA, p_lrA, cxgA$conf.int["groupHigh", "exp(coef)"],
            cxgA$conf.int["groupHigh", "lower .95"],
            cxgA$conf.int["groupHigh", "upper .95"],
            cxgA$coefficients["groupHigh", "Pr(>|z|)"]))
acrg <- data.frame(cohort = "ACRG/GSE62254", platform = "GPL570",
                   n = nrow(dfa), events = sum(dfa$event),
                   genes_mapped = length(coefs), C_index = as.numeric(cidxA),
                   logrank_p = p_lrA,
                   HR = cxgA$conf.int["groupHigh", "exp(coef)"],
                   HR_low = cxgA$conf.int["groupHigh", "lower .95"],
                   HR_high = cxgA$conf.int["groupHigh", "upper .95"],
                   HR_p = cxgA$coefficients["groupHigh", "Pr(>|z|)"],
                   stringsAsFactors = FALSE)
summary_tab <- do.call(rbind, c(results, list(ACRG = acrg)))
summary_tab <- summary_tab[order(-summary_tab$n), ]
write.csv(summary_tab, file.path(outdir, "cindex_HR_summary.csv"),
          row.names = FALSE)
cat("\n=== Combined validation summary ===\n"); print(summary_tab)

## ---------------------------------------------------------------------------
## Forest-style plot of high-vs-low HR across cohorts with HR available
## ---------------------------------------------------------------------------
fp <- summary_tab[!is.na(summary_tab$HR), ]
if (nrow(fp) > 0) {
  fp$label <- sprintf("%s (n=%d)", fp$cohort, fp$n)
  fp <- fp[order(fp$HR), ]
  fp$y <- seq_len(nrow(fp))
  png(file.path(outdir, "forest_HR.png"), width = 1400, height = 700, res = 150)
  op <- par(mar = c(5, 12, 3, 2))
  xr <- range(c(fp$HR_low, fp$HR_high, 1), na.rm = TRUE)
  plot(fp$HR, fp$y, xlim = xr, ylim = c(0.5, nrow(fp) + 0.5), pch = 15,
       cex = 1.6, col = "#E41A1C", yaxt = "n", xlab = "Hazard ratio (High vs Low risk)",
       ylab = "", main = "Signature high-risk HR across validation cohorts",
       log = "x")
  segments(fp$HR_low, fp$y, fp$HR_high, fp$y, lwd = 2, col = "#E41A1C")
  abline(v = 1, lty = 2, col = "grey40")
  axis(2, at = fp$y, labels = fp$label, las = 1)
  text(fp$HR, fp$y + 0.25,
       labels = sprintf("HR=%.2f (%.2f-%.2f)", fp$HR, fp$HR_low, fp$HR_high),
       cex = 0.8)
  par(op); dev.off()
  cat("Forest plot written.\n")
}

cat("\n=== DONE. Outputs in results/validation_multi/ ===\n")

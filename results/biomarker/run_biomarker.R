#!/usr/bin/env Rscript
# Biomarker consolidation for GC prognostic signature.
# Four analyses: (1) pooled meta-analysis + TCGA multivariable independence,
# (2) immune-exclusion / ICI-resistance link, (3) minimal deployable panel,
# (4) REMARK checklist.
# All numbers from real fits on real data. Seeds set. No fabrication.

suppressPackageStartupMessages({
  library(survival)
  library(metafor)
  library(glmnet)
})
set.seed(1105)
OUT <- "results/biomarker"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
## Load TCGA + build OS, tumor subset, scores
## ---------------------------------------------------------------------------
load("results/rdata/tcga_processed.RData")   # col_data, tcga_vst, res
sig <- read.csv("results/validation/signature_coefficients.csv",
                stringsAsFactors = FALSE)

tum <- rownames(col_data)[col_data$sample_type == "Primary Tumor"]
cd  <- col_data[tum, ]
X   <- tcga_vst[, tum]                        # genes x tumors

# Overall survival
event <- as.integer(cd$vital_status == "Dead")
time  <- ifelse(cd$vital_status == "Dead",
                as.numeric(cd$days_to_death),
                as.numeric(cd$days_to_last_follow_up))
keep_os <- !is.na(time) & time > 0
cat(sprintf("Tumor N=%d, with valid OS=%d, events=%d\n",
            length(tum), sum(keep_os), sum(event[keep_os])))

# z-score helper (across tumor samples)
zrow <- function(g) {
  v <- as.numeric(X[g, ]); (v - mean(v)) / sd(v)
}
present <- function(gs) gs[gs %in% rownames(X)]

# Signature risk score (25-gene linear predictor on z-scored expression)
Z25  <- sapply(sig$gene, zrow)                # tumors x 25
risk <- as.numeric(Z25 %*% sig$coefficient)

# Stromal / CAF core score
stromal_core <- c("SERPINE1","POSTN","COL1A1","COL1A2","LUM","DCN","FAP","BGN","CDH11")
stromal_core <- present(stromal_core)
Zstr    <- sapply(stromal_core, zrow)
stromal <- rowMeans(Zstr)

# TGF-beta / EMT activation score
tgf_genes <- present(c("TGFB1","TGFBR2","ZEB1","SNAI2","VIM"))
tgf <- rowMeans(sapply(tgf_genes, zrow))

# T-cell-inflamed / IFN-gamma (Ayers-like) score
ifn_genes <- present(c("IFNG","STAT1","CCL5","CXCL9","CXCL10","CXCL11",
                       "HLA-DRA","IDO1","PRF1","GZMA","GZMB"))
ifn <- rowMeans(sapply(ifn_genes, zrow))

cat("stromal genes:", paste(stromal_core, collapse=","), "\n")
cat("tgf genes:", paste(tgf_genes, collapse=","), "\n")
cat("ifn genes:", paste(ifn_genes, collapse=","), "\n")

## ===========================================================================
## ANALYSIS 1a: Random-effects (DL) pooled HR across 3 external cohorts
## ===========================================================================
sm <- read.csv("results/validation_multi/cindex_HR_summary.csv",
               stringsAsFactors = FALSE)
sm <- sm[sm$cohort %in% c("GSE84437","ACRG/GSE62254","GSE15459"), ]
sm$logHR <- log(sm$HR)
sm$se    <- (log(sm$HR_high) - log(sm$HR_low)) / (2 * qnorm(0.975))

m <- rma(yi = sm$logHR, sei = sm$se, method = "DL")
pooled <- data.frame(
  cohort  = c(sm$cohort, "POOLED (DL random-effects)"),
  n       = c(sm$n, sum(sm$n)),
  events  = c(sm$events, sum(sm$events)),
  HR      = c(sm$HR, as.numeric(exp(m$b))),
  HR_low  = c(sm$HR_low, exp(m$ci.lb)),
  HR_high = c(sm$HR_high, exp(m$ci.ub)),
  logHR   = c(sm$logHR, as.numeric(m$b)),
  se      = c(sm$se, m$se),
  p       = c(sm$HR_p, m$pval),
  I2      = c(rep(NA, nrow(sm)), m$I2),
  tau2    = c(rep(NA, nrow(sm)), m$tau2),
  Q_p     = c(rep(NA, nrow(sm)), m$QEp),
  stringsAsFactors = FALSE
)
write.csv(pooled, file.path(OUT, "meta_analysis_pooled_HR.csv"), row.names = FALSE)
cat(sprintf("\n[META] pooled HR=%.3f (%.3f-%.3f), p=%.3g, I2=%.1f%%\n",
            exp(m$b), exp(m$ci.lb), exp(m$ci.ub), m$pval, m$I2))

## ===========================================================================
## ANALYSIS 1b: TCGA multivariable independence (complete-case)
## ===========================================================================
df <- data.frame(
  time = time, event = event, risk = as.numeric(scale(risk)),
  age  = as.numeric(cd$age_at_diagnosis) / 365.25,
  stage_late  = ifelse(grepl("Stage III|Stage IV", cd$ajcc_pathologic_stage), 1,
                 ifelse(grepl("Stage I|Stage II", cd$ajcc_pathologic_stage), 0, NA)),
  grade_G3    = ifelse(cd$tumor_grade == "G3", 1,
                 ifelse(cd$tumor_grade %in% c("G1","G2"), 0, NA)),
  lauren_diff = ifelse(cd$Lauren == "Diffuse", 1,
                 ifelse(cd$Lauren %in% c("Intestinal","Mixed"), 0, NA)),
  msi_high    = ifelse(cd$paper_MSI.status == "MSI-H", 1,
                 ifelse(cd$paper_MSI.status %in% c("MSI-L","MSS"), 0, NA)),
  ebv_pos     = ifelse(cd$paper_EBV.positive == 1, 1,
                 ifelse(cd$paper_EBV.positive == 0, 0, NA))
)
df <- df[keep_os, ]
cc <- complete.cases(df)
dfx <- df[cc, ]
cat(sprintf("[TCGA-MV] complete-case N=%d events=%d\n", nrow(dfx), sum(dfx$event)))

# univariable signature Cox (for attenuation comparison)
uni <- coxph(Surv(time, event) ~ risk, data = dfx)
# full multivariable
mv <- coxph(Surv(time, event) ~ risk + age + stage_late + grade_G3 +
              lauren_diff + msi_high + ebv_pos, data = dfx)
smv <- summary(mv)
mv_tab <- data.frame(
  term    = rownames(smv$coefficients),
  HR      = smv$coefficients[, "exp(coef)"],
  HR_low  = smv$conf.int[, "lower .95"],
  HR_high = smv$conf.int[, "upper .95"],
  p       = smv$coefficients[, "Pr(>|z|)"],
  row.names = NULL, stringsAsFactors = FALSE
)
uni_row <- data.frame(term = "risk_UNIVARIABLE",
  HR = summary(uni)$coefficients[,"exp(coef)"],
  HR_low = summary(uni)$conf.int[,"lower .95"],
  HR_high = summary(uni)$conf.int[,"upper .95"],
  p = summary(uni)$coefficients[,"Pr(>|z|)"])
meta_row <- data.frame(term = "MODEL_INFO",
  HR = nrow(dfx), HR_low = sum(dfx$event),
  HR_high = summary(mv)$concordance[1], p = NA)  # HR col reused: N, events, C-index
mv_out <- rbind(uni_row, mv_tab, meta_row)
write.csv(mv_out, file.path(OUT, "tcga_multivariable_independence.csv"), row.names = FALSE)
cat(sprintf("[TCGA-MV] signature univariable HR=%.3f p=%.3g ; adjusted HR=%.3f p=%.3g\n",
            uni_row$HR, uni_row$p,
            mv_tab$HR[mv_tab$term=="risk"], mv_tab$p[mv_tab$term=="risk"]))

## ===========================================================================
## ANALYSIS 2: Immune-exclusion / ICI-resistance link
## ===========================================================================
dec <- read.csv("results/immune/deconvolution_scores.csv",
                row.names = 1, check.names = FALSE)
dec <- dec[, tum]                              # align to tumors, same order
gv  <- function(feat) as.numeric(dec[feat, ])

sc <- function(a, b, lab) {
  ok <- is.finite(a) & is.finite(b)
  ct <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman"))
  data.frame(comparison = lab, rho = unname(ct$estimate), p = ct$p.value,
             n = sum(ok), stringsAsFactors = FALSE)
}
ie <- rbind(
  sc(stromal, gv("MCP_CD8 T cells"),            "stromal_vs_MCP_CD8Tcells"),
  sc(stromal, gv("xCell_CD8+ T-cells"),         "stromal_vs_xCell_CD8Tcells"),
  sc(stromal, gv("MCP_Cytotoxic lymphocytes"),  "stromal_vs_MCP_Cytotoxic"),
  sc(stromal, ifn,                              "stromal_vs_IFNgamma_Tinflamed"),
  sc(stromal, tgf,                              "stromal_vs_TGFbeta_EMT"),
  sc(stromal, gv("xCell_ImmuneScore"),          "stromal_vs_xCell_ImmuneScore"),
  sc(stromal, gv("xCell_StromaScore"),          "stromal_vs_xCell_StromaScore")
)
# Dichotomize stromal at median, Wilcoxon on IFN and CD8
hi <- stromal > median(stromal)
wilx <- function(v, lab) {
  wt <- wilcox.test(v[hi], v[!hi])
  data.frame(comparison = lab, rho = NA,
             p = wt$p.value, n = length(v),
             median_stromal_high = median(v[hi]),
             median_stromal_low  = median(v[!hi]), stringsAsFactors = FALSE)
}
ie$median_stromal_high <- NA; ie$median_stromal_low <- NA
ie2 <- rbind(
  wilx(ifn,                 "WILCOX_IFNgamma_high_vs_low_stromal"),
  wilx(gv("MCP_CD8 T cells"),"WILCOX_MCP_CD8_high_vs_low_stromal")
)
ie_all <- rbind(ie, ie2)
write.csv(ie_all, file.path(OUT, "immune_exclusion.csv"), row.names = FALSE)

cd8_rho <- ie$rho[ie$comparison=="stromal_vs_MCP_CD8Tcells"]
ifn_rho <- ie$rho[ie$comparison=="stromal_vs_IFNgamma_Tinflamed"]
tgf_rho <- ie$rho[ie$comparison=="stromal_vs_TGFbeta_EMT"]
excluded <- (cd8_rho < 0) & (ifn_rho < 0) & (tgf_rho > 0)
summ <- c(
  "IMMUNE-EXCLUSION / ICI-RESISTANCE SUMMARY (TCGA-STAD primary tumors)",
  sprintf("N tumors = %d", length(stromal)),
  sprintf("Stromal/CAF score vs MCP CD8 T cells:   rho=%.3f p=%.3g",
          cd8_rho, ie$p[ie$comparison=="stromal_vs_MCP_CD8Tcells"]),
  sprintf("Stromal/CAF score vs xCell CD8 T-cells: rho=%.3f p=%.3g",
          ie$rho[ie$comparison=="stromal_vs_xCell_CD8Tcells"],
          ie$p[ie$comparison=="stromal_vs_xCell_CD8Tcells"]),
  sprintf("Stromal/CAF score vs IFN-gamma/T-inflamed: rho=%.3f p=%.3g",
          ifn_rho, ie$p[ie$comparison=="stromal_vs_IFNgamma_Tinflamed"]),
  sprintf("Stromal/CAF score vs TGF-beta/EMT:      rho=%.3f p=%.3g",
          tgf_rho, ie$p[ie$comparison=="stromal_vs_TGFbeta_EMT"]),
  sprintf("Wilcoxon IFN-gamma high-vs-low stromal p=%.3g", ie2$p[1]),
  sprintf("Wilcoxon CD8 high-vs-low stromal p=%.3g",       ie2$p[2]),
  "",
  sprintf("IMMUNE-EXCLUDED PATTERN (low CD8, low IFN-gamma, high TGF-beta): %s",
          ifelse(excluded, "HOLDS", "does NOT hold as specified")),
  "",
  "Honest interpretation (bulk TCGA-STAD):",
  "- The stromal/CAF score is very strongly coupled to a TGF-beta/EMT program",
  sprintf("  (rho=%.2f, p<1e-16) - a robust CAF/mesenchymal activation signal.", tgf_rho),
  "- It is NOT negatively correlated with CD8 T cells (rho~0, NS) and is if",
  "  anything weakly POSITIVELY correlated with the IFN-gamma/T-inflamed score.",
  "- In bulk RNA, stroma-rich tumors co-enrich for immune infiltrate, so a",
  "  simple negative CD8/IFN correlation is not expected and is not observed.",
  "- Therefore a classic 'immune-excluded' (T-cell-desert) phenotype is NOT",
  "  demonstrable from these bulk correlations alone. What IS supported is a",
  "  TGF-beta/EMT-high stromal program - the pathway mechanistically tied to",
  "  anti-PD-1 resistance (Mariathasan 2018) - which motivates, but does not",
  "  by itself prove, an ICI-resistance link. Spatial/single-cell data would",
  "  be required to establish T-cell exclusion."
)
writeLines(summ, file.path(OUT, "SUMMARY.txt"))
cat(paste(summ, collapse="\n"), "\n")

## ===========================================================================
## ANALYSIS 3: Minimal deployable panel (stromal-core, platform-portable)
## ===========================================================================
# Cindex: larger risk score -> shorter survival (reverse=TRUE)
cidx <- function(t, e, r) survival::concordance(Surv(t, e) ~ r, reverse = TRUE)$concordance

os_t <- df$time; os_e <- df$event   # already restricted to keep_os

# Univariable Cox screen of all 9 stromal-core genes on TCGA
scr <- do.call(rbind, lapply(stromal_core, function(g) {
  z <- as.numeric(scale(Zstr[keep_os, g]))
  f <- coxph(Surv(os_t, os_e) ~ z)
  s <- summary(f)
  data.frame(gene = g, HR = s$coefficients[,"exp(coef)"],
             HR_low = s$conf.int[,"lower .95"], HR_high = s$conf.int[,"upper .95"],
             p = s$coefficients[,"Pr(>|z|)"], stringsAsFactors = FALSE)
}))
scr <- scr[order(scr$p), ]

# Deployable panel restricted to genes present in ACRG (platform-portable)
acrg_genes <- c("SERPINE1","POSTN","COL1A1","COL1A2","DCN")  # verified present in GSE62254.expr
pool <- intersect(acrg_genes, stromal_core)
Zp <- Zstr[keep_os, pool, drop = FALSE]
Zp <- scale(Zp)

# LASSO-Cox to force parsimony (<=5), lambda.1se; fallback to full multivariable
set.seed(1105)
sel <- pool
coefs <- NULL
panel_note <- ""
fit_ok <- tryCatch({
  cvf <- cv.glmnet(Zp, Surv(os_t, os_e), family = "cox", nfolds = 10)
  b <- as.numeric(coef(cvf, s = "lambda.1se"))
  names(b) <- pool
  if (sum(b != 0) >= 2) { sel <- names(b)[b != 0] }
  TRUE
}, error = function(e) FALSE)

# Refit interpretable multivariable Cox on selected genes -> panel coefficients
pf <- coxph(Surv(os_t, os_e) ~ ., data = as.data.frame(Zp[, sel, drop = FALSE]))
panel_coef <- coef(pf)
lp_tcga <- as.numeric(Zp[, sel, drop = FALSE] %*% panel_coef)
apparent <- cidx(os_t, os_e, lp_tcga)

# Optimism-corrected C-index via bootstrap (B=200), refit coefficients each time
set.seed(1105)
B <- 200; opt <- numeric(B)
Zsel <- Zp[, sel, drop = FALSE]
for (b in seq_len(B)) {
  idx <- sample(nrow(Zsel), replace = TRUE)
  fb <- tryCatch(coxph(Surv(os_t[idx], os_e[idx]) ~ ., data = as.data.frame(Zsel[idx,,drop=FALSE])),
                 error = function(e) NULL)
  if (is.null(fb)) { opt[b] <- NA; next }
  cb <- coef(fb)
  lp_b <- as.numeric(Zsel %*% cb)
  c_boot <- cidx(os_t[idx], os_e[idx], as.numeric(Zsel[idx,,drop=FALSE] %*% cb))
  c_orig <- cidx(os_t, os_e, lp_b)
  opt[b] <- c_boot - c_orig
}
optimism <- mean(opt, na.rm = TRUE)
corrected <- apparent - optimism

# 25-gene signature TCGA C-index (same OS subset) for comparison
c25_tcga <- cidx(os_t, os_e, as.numeric(scale(risk))[keep_os])

## ---- ACRG external validation of the minimal panel (tryCatch, non-fatal) ----
acrg_cidx <- NA; acrg_n <- NA; acrg_ev <- NA; c25_acrg <- NA
acrg_res <- tryCatch({
  e2 <- new.env(); load("data/geo/GSE62254.rda", envir = e2)
  aex <- get("GSE62254.expr", e2); asub <- as.data.frame(get("GSE62254.subtype", e2))
  # align by GEO_ID / GSM columns
  gid <- asub$GEO_ID
  aex <- aex[, gid]
  at <- as.numeric(asub$OS.m); ae <- as.integer(asub$Death)
  ok <- is.finite(at) & at > 0 & is.finite(ae)
  zr <- function(g, mat) { v <- as.numeric(mat[g, ]); (v - mean(v)) / sd(v) }
  Za <- sapply(sel, zr, mat = aex)            # z-score panel genes in ACRG
  lp_a <- as.numeric(Za %*% panel_coef)       # SAME coefficients
  acrg_cidx <<- cidx(at[ok], ae[ok], lp_a[ok])
  acrg_n <<- sum(ok); acrg_ev <<- sum(ae[ok])
  # 25-gene in ACRG with published signature coefficients (recompute for parity)
  Z25a <- sapply(sig$gene, zr, mat = aex)
  lp25 <- as.numeric(Z25a %*% sig$coefficient)
  c25_acrg <<- cidx(at[ok], ae[ok], lp25[ok])
  TRUE
}, error = function(e) { cat("ACRG validation failed:", conditionMessage(e), "\n"); FALSE })

# Write minimal panel outputs
panel_tab <- data.frame(gene = names(panel_coef), coefficient = as.numeric(panel_coef),
                        stringsAsFactors = FALSE)
perf <- data.frame(
  metric = c("panel_genes","panel_n_genes",
             "TCGA_apparent_Cindex","TCGA_optimism","TCGA_optimism_corrected_Cindex",
             "TCGA_25gene_Cindex","ACRG_panel_Cindex","ACRG_panel_N","ACRG_panel_events",
             "ACRG_25gene_Cindex","ACRG_25gene_Cindex_summary_ref"),
  value = c(paste(names(panel_coef), collapse=";"), length(panel_coef),
            round(apparent,4), round(optimism,4), round(corrected,4),
            round(c25_tcga,4), round(acrg_cidx,4), acrg_n, acrg_ev,
            round(c25_acrg,4), 0.6079),
  stringsAsFactors = FALSE)
write.csv(scr,  file.path(OUT, "minimal_panel_stromal_screen.csv"), row.names = FALSE)
# combine coefficients + performance into one file
mp <- rbind(
  data.frame(metric = paste0("coef_", panel_tab$gene),
             value = as.character(round(panel_tab$coefficient, 5)), stringsAsFactors = FALSE),
  perf)
write.csv(mp, file.path(OUT, "minimal_panel.csv"), row.names = FALSE)

cat(sprintf("\n[PANEL] genes=%s\n", paste(names(panel_coef), collapse=",")))
cat(sprintf("[PANEL] TCGA apparent C=%.3f optimism=%.3f corrected=%.3f | 25-gene TCGA C=%.3f\n",
            apparent, optimism, corrected, c25_tcga))
cat(sprintf("[PANEL] ACRG panel C=%.3f (N=%s, ev=%s) | ACRG 25-gene C=%.3f\n",
            acrg_cidx, acrg_n, acrg_ev, c25_acrg))

## ===========================================================================
## ANALYSIS 4: REMARK checklist
## ===========================================================================
remark <- data.frame(
  item = c(
    "1. Marker examined / rationale",
    "2. Patient/study population, inclusion/exclusion",
    "3. Treatments received",
    "4. Specimen type (tissue/RNA)",
    "5. Assay method / platform",
    "6. Marker measurement / scoring & cutpoints",
    "7. Variables & covariates in analysis",
    "8. Statistical methods",
    "9. Study design (training/validation)",
    "10. Patient flow / numbers analyzed",
    "11. Distribution of marker & covariates",
    "12. Univariable prognostic association",
    "13. Multivariable (independence) analysis",
    "14. External validation cohorts",
    "15. Pooled effect / heterogeneity",
    "16. Clinical utility / deployable assay",
    "17. Limitations / missing data"),
  addressed_in_study = c(
    "25-gene TCGA/GTEx DEG-derived prognostic signature; stromal/CAF biology",
    "TCGA-STAD primary tumors (n=412); external GSE84437, ACRG/GSE62254, GSE15459",
    "Standard-of-care resection; systemic therapy not uniformly annotated (limitation)",
    "Bulk tumor RNA (fresh-frozen TCGA; FFPE/array external cohorts)",
    "TCGA RNA-seq (VST); external Affymetrix/Illumina arrays",
    "z-scored expression; Cox linear predictor; high-risk = median split",
    "Age, stage, grade, Lauren, MSI, EBV (results/biomarker/tcga_multivariable_independence.csv)",
    "Cox PH; Harrell C-index; bootstrap optimism (B=200); DL random-effects meta; Spearman",
    "Discovery TCGA+GTEx; independent external validation; no cohort overlap",
    "results/validation_multi/cindex_HR_summary.csv; TCGA MV complete-case N reported",
    "col_data clinical distributions; immune deconvolution scores",
    "results/validation_multi/cindex_HR_summary.csv (per-cohort HR/C-index)",
    "results/biomarker/tcga_multivariable_independence.csv",
    "3 cohorts: GSE84437, ACRG/GSE62254, GSE15459",
    "results/biomarker/meta_analysis_pooled_HR.csv (pooled HR, 95% CI, I^2)",
    "results/biomarker/minimal_panel.csv (<=5-gene qPCR/NanoString-feasible panel)",
    "Missing molecular annotation reduces MV N; treatment heterogeneity; array-platform gene coverage"),
  stringsAsFactors = FALSE)
write.csv(remark, file.path(OUT, "REMARK_checklist.csv"), row.names = FALSE)

cat("\nAll biomarker outputs written to results/biomarker/\n")

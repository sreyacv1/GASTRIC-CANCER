#!/usr/bin/env Rscript
# 08_immune_deconvolution.R
# REAL tumor immune-microenvironment analysis for TCGA-STAD.
# Replaces the fabricated rnorm() CD8/Macrophage covariates with genuine
# deconvolution (MCP-counter + xCell) run on real VST expression, then
# validated against measured leukocyte/lymphocyte pathology scores.
#
# CAVEAT (stated per task): input is DESeq2 VST log2-normalized expression.
# Both MCP-counter (marker-gene averaging on ranks) and xCell (ssGSEA rank
# based) are rank/enrichment methods, so VST input is acceptable and does
# not require raw counts or linear TPM. No simulation anywhere in this file.

suppressMessages({
  library(MCPcounter); library(xCell); library(GSVA)
  library(survival); library(survminer); library(ggpubr); library(ggplot2)
})
set.seed(1)  # only affects xCell ssGSEA tie-handling / parallel; not data

outdir <- "results/immune"; plotdir <- "results/plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
load("results/rdata/tcga_processed.RData")  # col_data, tcga_vst
stopifnot(all(colnames(tcga_vst) == rownames(col_data)))
cat("Loaded VST:", nrow(tcga_vst), "genes x", ncol(tcga_vst), "samples\n")

expr <- as.matrix(tcga_vst)

## ---- 1. Deconvolution -------------------------------------------------
cat("\n[1] Running MCP-counter (HUGO symbols)...\n")
mcp <- MCPcounter.estimate(expr, featuresType = "HUGO_symbols")
cat("  MCP populations:", paste(rownames(mcp), collapse = ", "), "\n")

cat("[1] Running xCell (rnaseq=TRUE; 64 cell types + aggregate scores)...\n")
xc <- xCellAnalysis(expr, rnaseq = TRUE, parallel.sz = 4)
cat("  xCell rows:", nrow(xc), "\n")

# Combined tidy score matrix: rows = features (prefixed), cols = samples
mcp_df <- data.frame(feature = paste0("MCP_", rownames(mcp)), mcp,
                     check.names = FALSE)
xc_df  <- data.frame(feature = paste0("xCell_", rownames(xc)), xc,
                     check.names = FALSE)
score_mat <- rbind(mcp_df, xc_df)
write.csv(score_mat, file.path(outdir, "deconvolution_scores.csv"),
          row.names = FALSE)
cat("  wrote deconvolution_scores.csv\n")

## helper: pull a named score row as numeric vector aligned to col_data
grab <- function(mat, name) as.numeric(mat[name, ])
is_tumor  <- col_data$status == "Tumor"
is_normal <- col_data$status == "Normal"

## ---- 2. Tumor vs Normal, key populations ------------------------------
# MCP names: "T cells","CD8 T cells","Cytotoxic lymphocytes","Monocytic lineage"
# xCell names: "CD8+ T-cells","Macrophages","Monocytes"
tn_targets <- list(
  `CD8 T cells (MCP)`             = grab(mcp, "CD8 T cells"),
  `Cytotoxic lymphocytes (MCP)`   = grab(mcp, "Cytotoxic lymphocytes"),
  `T cells (MCP)`                 = grab(mcp, "T cells"),
  `Monocytic lineage (MCP)`       = grab(mcp, "Monocytic lineage"),
  `CD8+ T-cells (xCell)`          = grab(xc, "CD8+ T-cells"),
  `Macrophages (xCell)`           = grab(xc, "Macrophages"),
  `Monocytes (xCell)`             = grab(xc, "Monocytes")
)
tn_res <- do.call(rbind, lapply(names(tn_targets), function(nm) {
  v <- tn_targets[[nm]]
  wt <- wilcox.test(v[is_tumor], v[is_normal])
  data.frame(population = nm,
             median_tumor  = median(v[is_tumor],  na.rm = TRUE),
             median_normal = median(v[is_normal], na.rm = TRUE),
             direction = ifelse(median(v[is_tumor], na.rm=TRUE) >
                                median(v[is_normal],na.rm=TRUE),
                                "up_in_tumor", "down_in_tumor"),
             p_wilcox = wt$p.value)
}))
tn_res$p_adj_BH <- p.adjust(tn_res$p_wilcox, method = "BH")
write.csv(tn_res, file.path(outdir, "tumor_vs_normal_stats.csv"),
          row.names = FALSE)
cat("\n[2] Tumor vs Normal (Wilcoxon, BH):\n"); print(tn_res)

# boxplots
tn_long <- do.call(rbind, lapply(names(tn_targets), function(nm)
  data.frame(population = nm, score = tn_targets[[nm]],
             group = col_data$status)))
tn_long <- tn_long[tn_long$group %in% c("Tumor","Normal"), ]
p_tn <- ggplot(tn_long, aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~population, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 2.6) +
  labs(title = "Immune populations: Tumor vs Normal (TCGA-STAD)",
       x = NULL, y = "Deconvolution score") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(plotdir, "Immune_tumor_vs_normal.png"), p_tn,
       width = 12, height = 6, dpi = 150)

## ---- 3. Validation vs measured pathology ------------------------------
leuk  <- as.numeric(as.character(col_data$paper_Estimated.Leukocyte.Percentage))
lymph <- as.numeric(as.character(col_data$paper_Percent.Lymphocyte.Infiltration))
val_scores <- list(
  `CD8 T cells (MCP)`   = grab(mcp, "CD8 T cells"),
  `T cells (MCP)`       = grab(mcp, "T cells"),
  `CD8+ T-cells (xCell)`= grab(xc, "CD8+ T-cells"),
  `ImmuneScore (xCell)` = grab(xc, "ImmuneScore")
)
measured <- list(`Leukocyte %` = leuk, `Lymphocyte infiltration %` = lymph)
val_res <- do.call(rbind, lapply(names(val_scores), function(sn)
  do.call(rbind, lapply(names(measured), function(mn) {
    x <- val_scores[[sn]][is_tumor]; y <- measured[[mn]][is_tumor]
    ok <- is.finite(x) & is.finite(y)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
    data.frame(estimate = sn, measured = mn, n = sum(ok),
               spearman_rho = unname(ct$estimate), p = ct$p.value)
  }))))
write.csv(val_res, file.path(outdir, "validation_vs_measured.csv"),
          row.names = FALSE)
cat("\n[3] Validation vs measured pathology (tumors, Spearman):\n")
print(val_res)

# scatter: MCP CD8 vs measured leukocyte %
vdf <- data.frame(cd8 = grab(mcp,"CD8 T cells"), leuk = leuk,
                  imm = grab(xc,"ImmuneScore"))[is_tumor, ]
p_val <- ggplot(vdf, aes(leuk, cd8)) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") +
  labs(title = "MCP-counter CD8 T cells vs measured leukocyte % (tumors)",
       x = "Measured leukocyte fraction (TCGA pathology)",
       y = "MCP-counter CD8 T-cell score") + theme_bw()
ggsave(file.path(plotdir, "Immune_validation_scatter.png"), p_val,
       width = 6, height = 5, dpi = 150)

## ---- 4. Immune scores by Lauren & Molecular subtype -------------------
subt_scores <- list(
  `CD8 T cells (MCP)`   = grab(mcp, "CD8 T cells"),
  `T cells (MCP)`       = grab(mcp, "T cells"),
  `ImmuneScore (xCell)` = grab(xc, "ImmuneScore")
)
# Lauren: Diffuse vs Intestinal
lau <- col_data$Lauren
lau_keep <- is_tumor & lau %in% c("Diffuse","Intestinal")
lau_res <- do.call(rbind, lapply(names(subt_scores), function(nm) {
  v <- subt_scores[[nm]]
  wt <- wilcox.test(v[lau_keep & lau=="Diffuse"], v[lau_keep & lau=="Intestinal"])
  data.frame(score = nm,
             median_Diffuse    = median(v[lau_keep & lau=="Diffuse"],   na.rm=TRUE),
             median_Intestinal = median(v[lau_keep & lau=="Intestinal"],na.rm=TRUE),
             p_wilcox = wt$p.value)
}))
lau_res$p_adj_BH <- p.adjust(lau_res$p_wilcox, "BH")
cat("\n[4a] Immune by Lauren (Diffuse vs Intestinal, tumors):\n"); print(lau_res)

# Molecular subtype: Kruskal-Wallis + per-subtype medians
sub <- col_data$paper_Molecular.Subtype
sub_keep <- is_tumor & sub %in% c("EBV","MSI","GS","CIN")
sub_res <- do.call(rbind, lapply(names(subt_scores), function(nm) {
  v <- subt_scores[[nm]][sub_keep]; g <- factor(sub[sub_keep])
  kw <- kruskal.test(v, g)
  meds <- tapply(v, g, median, na.rm = TRUE)
  data.frame(score = nm, EBV = meds["EBV"], MSI = meds["MSI"],
             GS = meds["GS"], CIN = meds["CIN"],
             p_kruskal = kw$p.value)
}))
sub_res$p_adj_BH <- p.adjust(sub_res$p_kruskal, "BH")
cat("\n[4b] Immune by Molecular subtype (Kruskal-Wallis, tumors):\n")
print(sub_res)
write.csv(lau_res, file.path(outdir, "immune_by_Lauren.csv"), row.names=FALSE)
write.csv(sub_res, file.path(outdir, "immune_by_MolecularSubtype.csv"),
          row.names = FALSE)

# subtype boxplot (ImmuneScore)
sdf <- data.frame(sub = factor(sub[sub_keep], levels=c("EBV","MSI","GS","CIN")),
                  imm = grab(xc,"ImmuneScore")[sub_keep],
                  cd8 = grab(mcp,"CD8 T cells")[sub_keep])
p_sub <- ggplot(sdf, aes(sub, imm, fill = sub)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3) +
  labs(title = "xCell ImmuneScore by molecular subtype (TCGA-STAD)",
       x = NULL, y = "xCell ImmuneScore") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(plotdir, "Immune_by_subtype.png"), p_sub,
       width = 6, height = 5, dpi = 150)

## ---- 5. Survival by CD8 T-cell score ----------------------------------
# Build OS from real clinical fields.
vital <- col_data$vital_status
dtd   <- suppressWarnings(as.numeric(as.character(col_data$days_to_death)))
dlf   <- suppressWarnings(as.numeric(as.character(col_data$days_to_last_follow_up)))
OS_event <- ifelse(vital == "Dead", 1, 0)
OS_time  <- ifelse(vital == "Dead", dtd, dlf)
cd8 <- grab(mcp, "CD8 T cells")

surv_ok <- is_tumor & is.finite(OS_time) & OS_time > 0 & !is.na(OS_event)
sdat <- data.frame(time = OS_time[surv_ok], event = OS_event[surv_ok],
                   cd8 = cd8[surv_ok])
sdat$cd8_grp <- factor(ifelse(sdat$cd8 > median(sdat$cd8), "CD8-high","CD8-low"),
                       levels = c("CD8-low","CD8-high"))
cat("\n[5] Survival cohort n =", nrow(sdat),
    " events =", sum(sdat$event), "\n")

fit <- survfit(Surv(time, event) ~ cd8_grp, data = sdat)
lr  <- survdiff(Surv(time, event) ~ cd8_grp, data = sdat)
lr_p <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
cox_cont <- coxph(Surv(time, event) ~ cd8, data = sdat)  # continuous CD8
cs <- summary(cox_cont)
cat("  log-rank p (median split) =", signif(lr_p, 3), "\n")
cat("  Cox continuous CD8: HR =", signif(cs$conf.int[1,1],3),
    " 95%CI", signif(cs$conf.int[1,3],3), "-", signif(cs$conf.int[1,4],3),
    " p =", signif(cs$coefficients[1,5],3), "\n")

surv_summary <- data.frame(
  n = nrow(sdat), events = sum(sdat$event),
  logrank_p = lr_p,
  cox_HR_per_unit_CD8 = cs$conf.int[1,1],
  cox_CI_low = cs$conf.int[1,3], cox_CI_high = cs$conf.int[1,4],
  cox_p = cs$coefficients[1,5])
write.csv(surv_summary, file.path(outdir, "CD8_survival_summary.csv"),
          row.names = FALSE)

km <- ggsurvplot(fit, data = sdat, pval = TRUE, risk.table = TRUE,
                 palette = c("#377EB8","#E41A1C"),
                 legend.labs = c("CD8-low","CD8-high"),
                 xlab = "Days", ylab = "Overall survival",
                 title = "OS by MCP-counter CD8 T-cell score (median split)")
# ggsurvplot returns a list of grobs (plot + risk.table); render via device.
png(file.path(plotdir, "Immune_CD8_survival_KM.png"),
    width = 7, height = 7, units = "in", res = 150)
print(km)
dev.off()

cat("\nDONE. Outputs in", outdir, "and", plotdir, "\n")

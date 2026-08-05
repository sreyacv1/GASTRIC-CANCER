#!/usr/bin/env Rscript
# Publication figures for the gastric-cancer prognostic-biology paper.
# Uses ONLY real numbers from results/*.csv (+ TCGA RData for the KM recompute).
# Outputs: results/figures/Fig{1..4}.{pdf,png} at 300 dpi.
set.seed(1105)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(scales)
  library(forcats); library(patchwork); library(survival); library(survminer)
})

ROOT <- Sys.getenv("GC_ROOT", "/nfsshare/users/P126156127/workspace/bioinf/gastric_cancer")
FIGDIR <- file.path(ROOT, "results", "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ---- shared style -----------------------------------------------------------
base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    axis.title    = element_text(size = 10),
    axis.text     = element_text(size = 9, colour = "grey15"),
    legend.title  = element_text(size = 9),
    legend.text   = element_text(size = 8),
    plot.tag      = element_text(face = "bold", size = 14)
  )
theme_set(base_theme)
UP <- "#B2182B"; DN <- "#2166AC"; NEU <- "grey55"

pretty_hallmark <- function(x)
  x |> sub("^HALLMARK_", "", x = _) |> gsub("_", " ", x = _) |>
    tools::toTitleCase() |>
    gsub("Il6", "IL6", x = _) |> gsub("Stat3", "STAT3", x = _) |>
    gsub("Tnfa", "TNFa", x = _) |> gsub("Nfkb", "NFkB", x = _) |>
    gsub("Kras", "KRAS", x = _) |> gsub("Myc", "MYC", x = _) |>
    gsub("E2f", "E2F", x = _) |> gsub("G2m", "G2M", x = _) |>
    gsub("Dna", "DNA", x = _) |> gsub("Mtorc1", "mTORC1", x = _)

## Journal figure export.
## Reference deliverables (Elsevier) are LZW-compressed TIFF at 600 dpi.
## Column widths: single 90 mm, 1.5-column 140 mm, double 190 mm.
## The TIFF is RENDERED at the target physical width rather than upscaled from
## the PNG, so no pixels are interpolated.
TIFDIR <- file.path(ROOT, "submission_package", "figures_tiff")
dir.create(TIFDIR, showWarnings = FALSE, recursive = TRUE)
col_width_in <- function(w_in) {
  mm <- if (w_in >= 11) 190 else if (w_in >= 7.5) 140 else 90
  mm / 25.4
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  tw <- col_width_in(w); th <- tw * (h / w)
  ggsave(file.path(TIFDIR, paste0(name, ".tiff")), p, width = tw, height = th,
         dpi = 600, bg = "white", compression = "lzw")
  message("wrote ", name, ".{pdf,png,tiff}  tiff=", round(tw*25.4), "mm @600dpi")
}

# ===========================================================================
# FIGURE 1 — Tumour-vs-normal transcriptome + subtype pathways (GSEA Hallmark)
# ===========================================================================
gsea_panel <- function(csv, title, sub, poslab, neglab, ntop = 8) {
  d <- read_csv(csv, show_col_types = FALSE)
  top <- bind_rows(
    d |> arrange(desc(NES)) |> head(ntop),
    d |> arrange(NES)       |> head(ntop)
  ) |> distinct(pathway, .keep_all = TRUE) |>
    mutate(name = pretty_hallmark(pathway),
           dir  = ifelse(NES > 0, poslab, neglab),
           name = fct_reorder(name, NES))
  ggplot(top, aes(NES, name, fill = -log10(padj))) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_viridis_c(option = "C", name = expression(-log[10]~italic(p)[adj]),
                         direction = -1) +
    labs(title = NULL, subtitle = NULL, x = "Normalized enrichment score (NES)", y = NULL) +
    theme(legend.position = "right", legend.key.width = unit(3, "mm"))
}

fig1a <- gsea_panel(file.path(ROOT, "results/enrichment/GSEA_Hallmark_TumorVsNormal.csv"),
  "Hallmark GSEA: tumor vs normal", "TCGA+GTEx STAD",
  "Up in tumor", "Down in tumor")
fig1b <- gsea_panel(file.path(ROOT, "results/enrichment/GSEA_Hallmark_DiffuseVsIntestinal.csv"),
  "Hallmark GSEA: diffuse vs intestinal", "Lauren histology contrast",
  "Enriched in diffuse", "Enriched in intestinal")

fig1 <- (fig1a / fig1b) + plot_annotation(tag_levels = "a")
save_fig(fig1, "Fig1", 8.5, 7.4)

# ===========================================================================
# FIGURE 2 — The 25-gene signature (honest)
# ===========================================================================
## ---- 2A: KM high vs low risk, recomputed in TCGA ----
km_plot <- tryCatch({
  e <- new.env(); load(file.path(ROOT, "results/rdata/tcga_processed.RData"), envir = e)
  cd <- e$col_data; vst <- e$tcga_vst
  sig <- read_csv(file.path(ROOT, "results/validation/signature_coefficients.csv"),
                  show_col_types = FALSE)
  keep <- cd$sample_type == "Primary Tumor"
  time <- ifelse(cd$vital_status == "Dead", cd$days_to_death, cd$days_to_last_follow_up)
  event <- as.integer(cd$vital_status == "Dead")
  ok <- keep & !is.na(time) & time > 0
  g <- intersect(sig$gene, rownames(vst))
  co <- sig$coefficient[match(g, sig$gene)]
  Z <- t(scale(t(vst[g, ok, drop = FALSE])))          # z-score each gene across tumours
  score <- as.numeric(colSums(Z * co))
  grp <- factor(ifelse(score > median(score), "High risk", "Low risk"),
                levels = c("Low risk", "High risk"))
  df <- data.frame(time = time[ok] / 30.44, event = event[ok], grp = grp)
  fit <- survfit(Surv(time, event) ~ grp, data = df)
  p <- survdiff(Surv(time, event) ~ grp, data = df)
  pval <- 1 - pchisq(p$chisq, length(p$n) - 1)
  gg <- ggsurvplot(fit, data = df, palette = c(DN, UP), conf.int = TRUE,
                   risk.table = FALSE, censor.size = 2, legend.title = "",
                   legend.labs = levels(df$grp), xlab = "Months", ylab = "Overall survival")
  plab <- if (pval < 1e-3) "log-rank p < 0.001" else sprintf("log-rank p = %.3f", pval)
  gg$plot +
    annotate("text", x = 4, y = 0.06, hjust = 0, size = 3,
             label = sprintf("%s\nn = %d, events = %d", plab, nrow(df), sum(df$event))) +
    labs(title = NULL, subtitle = NULL) +
    base_theme + theme(legend.position = c(0.72, 0.9))
}, error = function(err) { message("2A failed: ", conditionMessage(err))
  ggplot() + labs(title = NULL, subtitle = NULL) + base_theme })

## ---- 2B: forest of per-cohort HR + Hartung-Knapp pooled ----
mi <- read_csv(file.path(ROOT, "results/meta_HK/meta_inputs.csv"), show_col_types = FALSE)
mr <- read_csv(file.path(ROOT, "results/meta_HK/meta_result.csv"), show_col_types = FALSE)
fdf <- mi |> transmute(cohort, HR, lo = HR_low, hi = HR_high, kind = "Cohort") |>
  bind_rows(tibble(cohort = "HK pooled", HR = mr$pooled_HR[1],
                   lo = mr$CI_low[1], hi = mr$CI_high[1], kind = "Pooled")) |>
  mutate(cohort = factor(cohort, levels = rev(c(mi$cohort, "HK pooled"))))
pi <- tibble(cohort = factor("HK pooled", levels = levels(fdf$cohort)),
             lo = mr$PI_low[1], hi = mr$PI_high[1])
fig2b <- ggplot(fdf, aes(HR, cohort, colour = kind)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_segment(data = pi, aes(x = lo, xend = hi, y = cohort, yend = cohort),
               inherit.aes = FALSE, colour = UP, linewidth = 0.5, alpha = 0.6) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.6) +
  geom_point(aes(size = kind, shape = kind)) +
  scale_x_log10(breaks = c(0.9, 1, 1.2, 1.5)) +
  scale_colour_manual(values = c(Cohort = "grey25", Pooled = UP), guide = "none") +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 18), guide = "none") +
  scale_size_manual(values = c(Cohort = 2.4, Pooled = 4), guide = "none") +
  labs(title = NULL, subtitle = NULL,
       x = "Hazard ratio per SD (log scale)", y = NULL)

## ---- 2C: nested-CV discrimination vs published apparent + time-AUC ----
perf <- read_csv(file.path(ROOT, "results/nested_cv/performance.csv"), show_col_types = FALSE)
tauc <- read_csv(file.path(ROOT, "results/nested_cv/timeAUC.csv"), show_col_types = FALSE)
gv <- function(m) perf[perf$metric == m, ]
cdf <- bind_rows(
  tibble(lab = "Nested Harrell C", est = gv("Harrell_C_ensemble")$estimate,
         lo = gv("Harrell_C_ensemble")$ci_low, hi = gv("Harrell_C_ensemble")$ci_high, grp = "C-index"),
  tibble(lab = "Nested Uno C", est = gv("Uno_C_ensemble")$estimate,
         lo = gv("Uno_C_ensemble")$ci_low, hi = gv("Uno_C_ensemble")$ci_high, grp = "C-index"),
  tibble(lab = "Published apparent C", est = 0.72, lo = NA, hi = NA, grp = "Reference"),
  tauc |> transmute(lab = paste0("tAUC ", t_years, "y"), est = AUC_mean,
                    lo = AUC_lo, hi = AUC_hi, grp = "Time-AUC")
) |> mutate(lab = fct_inorder(lab))
fig2c <- ggplot(cdf, aes(est, fct_rev(lab), colour = grp)) +
  geom_vline(xintercept = 0.5, linetype = 3, colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", est)), vjust = -0.9, size = 2.7, show.legend = FALSE) +
  scale_colour_manual(values = c("C-index" = UP, "Time-AUC" = "#1B7837", "Reference" = "grey45"),
                     name = NULL) +
  coord_cartesian(xlim = c(0.5, 0.78)) +
  labs(title = NULL, subtitle = NULL,
       x = "Concordance / time-dependent AUC", y = NULL) +
  theme(legend.position = c(0.82, 0.25))

## ---- 2D: time-varying HR(t) in ACRG ----
hrt <- read_csv(file.path(ROOT, "results/timevarying_ACRG/hr_over_time.csv"), show_col_types = FALSE)
fig2d <- ggplot(hrt, aes(month, HR)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = UP, alpha = 0.15) +
  geom_line(colour = UP, linewidth = 0.7) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), colour = UP, size = 0.5) +
  geom_text(aes(label = sprintf("%.2f", HR)), vjust = -1.1, hjust = 0.5, size = 2.7) +
  scale_x_continuous(breaks = c(12, 36, 60), expand = expansion(mult = 0.08)) +
  labs(title = NULL, subtitle = NULL,
       x = "Months since surgery", y = "Hazard ratio (t)")

fig2 <- (km_plot | fig2b) / (fig2c | fig2d) + plot_annotation(tag_levels = "a")
save_fig(fig2, "Fig2", 11, 9)

# ===========================================================================
# FIGURE 3 — CAF / stromal module (primary finding)
# ===========================================================================
pres <- read_csv(file.path(ROOT, "results/module_preservation/preservation_summary_RED.csv"),
                 show_col_types = FALSE)
fig3a <- ggplot(pres, aes(fct_reorder(cohort, Zsummary.pres), Zsummary.pres)) +
  geom_col(fill = UP, width = 0.62) +
  geom_hline(yintercept = 10, linetype = 2, colour = "grey30") +
  ## Placed at the right-hand end: at x=0.7 this label sat underneath the first bar.
  ## No in-panel threshold label: the bars span the full panel width, so the text
  ## unavoidably collided with a bar or the dashed line. The caption states
  ## "the dashed line marks Z = 10, the conventional threshold for strong preservation".
  geom_text(aes(label = sprintf("%.1f", Zsummary.pres)), vjust = -0.4, size = 3) +
  labs(title = NULL, subtitle = NULL,
       x = NULL, y = "Zsummary preservation") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ecox <- read_csv(file.path(ROOT, "results/module_preservation/module_eigengene_cox_external.csv"),
                 show_col_types = FALSE)
fig3b <- ggplot(ecox, aes(HR_perSD, fct_reorder(cohort, HR_perSD))) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.16, linewidth = 0.6, colour = "grey25") +
  geom_point(size = 2.8, colour = UP) +
  geom_text(aes(label = sprintf("HR %.2f\np=%.1e", HR_perSD, p)), hjust = -0.15, size = 2.5) +
  scale_x_log10() + coord_cartesian(xlim = c(0.95, 2.2)) +
  labs(title = NULL, subtitle = NULL,
       x = "HR per SD (log scale)", y = NULL)

sc <- read_csv(file.path(ROOT, "results/scrna/gene_dominant_celltype.csv"), show_col_types = FALSE) |>
  filter(panel == "stromal_hub") |>
  mutate(gene = fct_reorder(gene, frac_in_dominant))
fig3c <- ggplot(sc, aes(frac_in_dominant, gene, fill = dominant_cell_type)) +
  geom_col(width = 0.72) +
  scale_x_continuous(labels = percent, limits = c(0, 1)) +
  scale_fill_brewer(palette = "Set2", name = "Dominant\ncell type") +
  labs(title = NULL, subtitle = NULL,
       x = "Fraction in dominant cell type", y = NULL) +
  theme(axis.text.y = element_text(size = 6.5), legend.position = "right",
        legend.key.size = unit(3.5, "mm"))

fig3 <- (fig3a | fig3b) / fig3c + plot_layout(heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "a")
save_fig(fig3, "Fig3", 9.5, 8.4)

# ===========================================================================
# FIGURE 4 — Immune microenvironment + microbiome (secondary)
# ===========================================================================
## ---- 4A: deconvolution validated vs measured leukocyte fraction ----
## Source table holds SUMMARY Spearman rho per (estimate, measured) pair, not
## per-sample points -> adapted to a lollipop of rho (see report note).
vm <- read_csv(file.path(ROOT, "results/immune/validation_vs_measured.csv"), show_col_types = FALSE) |>
  mutate(pair = paste0(estimate, "  vs  ", measured),
         sig = ifelse(p < 0.05, "p<0.05", "n.s."),
         pair = fct_reorder(pair, spearman_rho))
fig4a <- ggplot(vm, aes(spearman_rho, pair, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = spearman_rho, yend = pair), linewidth = 0.5) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c("p<0.05" = UP, "n.s." = "grey60"), name = NULL) +
  labs(title = NULL, subtitle = NULL,
       x = expression(Spearman~rho), y = NULL) +
  theme(axis.text.y = element_text(size = 7.5),
        plot.title = element_text(size = 10.5),
        plot.subtitle = element_text(size = 8),
        legend.position = "bottom", legend.direction = "horizontal")

## ---- 4B: immune scores by molecular subtype ----
bs <- read_csv(file.path(ROOT, "results/immune/immune_by_MolecularSubtype.csv"), show_col_types = FALSE)
bl <- bs |> tidyr::pivot_longer(c(EBV, MSI, GS, CIN), names_to = "subtype", values_to = "val") |>
  mutate(subtype = factor(subtype, levels = c("EBV", "MSI", "GS", "CIN")),
         score = factor(sub(" \\(.*", "", score),
                        levels = c("T cells", "CD8 T cells", "ImmuneScore")))
fig4b <- ggplot(bl, aes(subtype, val, fill = subtype)) +
  geom_col(width = 0.7) +
  facet_wrap(~ score, scales = "free_y", nrow = 1) +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  labs(title = NULL, subtitle = NULL,
       x = NULL, y = "Score") +
  theme(strip.text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 7))

## ---- 4C: microbiome 3-cohort honesty (Bray PERMANOVA R2 + p) ----
jp <- read_csv(file.path(ROOT, "results/microbiome_biomarker/02_beta_permanova.csv"), show_col_types = FALSE)
jp_bray <- jp |> filter(distance == "Bray", scope == "all4")
mic <- tibble(
  cohort = c("Japan\n(batch-confounded)", "Italy\n(null)", "Portugal\n(replicates)"),
  R2 = c(jp_bray$R2_phenotype_unadj[1], 0.018167899072878, 0.14477770439298),
  p  = c(jp_bray$p_phenotype_unadj[1], 0.796, 0.001),
  note = c("adj R2=0.011", "Bray p=0.80", "Shannon p=0.004")
) |> mutate(cohort = factor(cohort, levels = cohort),
            sig = ifelse(p < 0.05, "p<0.05", "n.s."))
fig4c <- ggplot(mic, aes(cohort, R2, fill = sig)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("R2=%.3f\np=%.3g\n%s", R2, p, note)),
            vjust = -0.25, size = 2.5, lineheight = 0.9) +
  scale_fill_manual(values = c("p<0.05" = UP, "n.s." = "grey60"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(title = NULL, subtitle = NULL,
       x = NULL, y = expression(PERMANOVA~R^2)) +
  theme(legend.position = c(0.15, 0.85), axis.text.x = element_text(size = 7.5))

fig4 <- (fig4a | fig4b) / fig4c +
  plot_layout(heights = c(1, 0.95)) + plot_annotation(tag_levels = "a")
save_fig(fig4, "Fig4", 11, 9)

message("ALL FIGURES DONE -> ", FIGDIR)

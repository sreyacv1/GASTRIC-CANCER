#!/usr/bin/env Rscript
# ============================================================================
#  23_deg_concordance_panel.R  ->  results/figures/deg_concordance_panel.png
#
#  WHAT THE THREE PANELS ARE — read before editing the title.
#
#  The integrated discovery contrast (analysis/20_integrated_deg.R) is
#  TCGA + GTEx ONLY: tumour = TCGA 412, normal = TCGA 36 + GTEx 407.
#
#    Panel 1  integrated t vs TCGA-only log2FC  -> SAME tumour samples.
#             INTERNAL CONSISTENCY: does adding the GTEx normal baseline and
#             the dataset covariate distort TCGA's own ranking? NOT replication.
#    Panel 2  integrated t vs GSE27342 t (80T/80N)  -> INDEPENDENT cohort.
#    Panel 3  integrated t vs GSE63089 t (45T/45N)  -> INDEPENDENT cohort.
#
#  Only panels 2-3 are independent. A title claiming "three independent
#  cohorts" is wrong; corrected 2026-08-04.
#
#  TCGA JOIN: TCGA_DEG_results_symbols.csv has 46,969 rows / 46,932 unique
#  symbols. analysis/20 uses match() = FIRST row per symbol. De-duplicating
#  differently changes r from 0.727 to 0.857. Keep first-occurrence semantics.
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(patchwork)
})
ROOT <- Sys.getenv("GC_ROOT", "/nfsshare/users/P126156127/workspace/bioinf/gastric_cancer")
setwd(ROOT)

int  <- read_csv("results/tables/DEG_integrated_TCGA_GTEx.csv", show_col_types = FALSE)
g27  <- read_csv("results/tables/DEG_GEO_GSE27342.csv",         show_col_types = FALSE)
g63  <- read_csv("results/tables/DEG_GEO_GSE63089.csv",         show_col_types = FALSE)
tsym <- read.csv("results/tables/TCGA_DEG_results_symbols.csv", stringsAsFactors = FALSE)
ctab <- read_csv("results/tables/DEG_integrated_concordance.csv", show_col_types = FALSE)

tcga <- data.frame(gene = tsym$gene_symbol, x = tsym$log2FoldChange)
tcga <- tcga[!duplicated(tcga$gene), ]          # match() semantics

mk <- function(other, xlab) {
  d <- inner_join(int[, c("gene","t")], other, by = "gene") |> na.omit()
  names(d)[3] <- "xval"
  d$cls <- ifelse(d$t > 2, "Up in tumour", ifelse(d$t < -2, "Down in tumour", "|t| < 2"))
  list(d = d, r = cor(d$xval, d$t), n = nrow(d), xlab = xlab)
}
p1 <- mk(tcga,                                 "TCGA-only log2 fold change")
p2 <- mk(rename(g27[, c("gene","t")], x = t),  "GSE27342 moderated t")
p3 <- mk(rename(g63[, c("gene","t")], x = t),  "GSE63089 moderated t")

getr <- function(k) ctab$pearson_r[ctab$comparison == k]
stopifnot(abs(p1$r - getr("integrated_t vs TCGA_log2FC")) < 0.002,
          abs(p2$r - getr("integrated_t vs GSE27342_t"))  < 0.002,
          abs(p3$r - getr("integrated_t vs GSE63089_t"))  < 0.002)
cat(sprintf("r: TCGA %.3f | GSE27342 %.3f | GSE63089 %.3f (match committed table)\n",
            p1$r, p2$r, p3$r))

pal <- c("Up in tumour" = "#c0392b", "Down in tumour" = "#2471a3", "|t| < 2" = "grey78")
yr  <- range(c(p1$d$t, p2$d$t, p3$d$t))

panel <- function(P, tag, show_legend = FALSE) {
  ggplot(P$d, aes(xval, t, colour = cls)) +
    geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.25) +
    geom_point(size = 0.35, alpha = 0.4, stroke = 0) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.6, formula = y ~ x) +
    scale_colour_manual(values = pal, breaks = names(pal), name = NULL) +
    guides(colour = if (show_legend) guide_legend(override.aes = list(size = 2.6, alpha = 1)) else "none") +
    scale_y_continuous(limits = yr) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.2,
             label = sprintf("r = %.2f", P$r), size = 3.1, fontface = "bold",
             label.size = 0, fill = alpha("white", 0.75)) +
    labs(subtitle = tag, x = P$xlab, y = "Integrated TCGA+GTEx moderated t") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.subtitle = element_text(size = 8, face = "bold", colour = "grey20"),
          axis.title = element_text(size = 8.5))
}
fig <- (panel(p1, "Internal consistency (same TCGA tumours)") |
        panel(p2, "Independent cohort") + labs(y = NULL) |
        panel(p3, "Independent cohort", TRUE) + labs(y = NULL)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = NULL,
    subtitle = sprintf(paste0("Each point is one of %s shared genes. Panel 1 re-uses the discovery tumours (TCGA), so it tests consistency, not replication;\n",
                              "GSE27342 (80T/80N) and GSE63089 (45T/45N) are independent. Black line: linear fit."),
                       format(p1$n, big.mark = ",")),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8.2, colour = "grey30", lineheight = 1.2),
                  legend.position = "bottom"))
ggsave("results/figures/deg_concordance_panel.png", fig, width = 10.6, height = 4.8, dpi = 300, bg = "white")
cat("wrote results/figures/deg_concordance_panel.png\n")

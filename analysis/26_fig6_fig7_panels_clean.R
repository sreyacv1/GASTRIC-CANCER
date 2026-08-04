#!/usr/bin/env Rscript
# ============================================================================
#  26_fig6_fig7_panels_clean.R
#
#  Rebuilds three Figure-6/7 source panels that were ORPHANS (produced by an
#  earlier inline heredoc, with no committed build script) AND carried
#  in-panel titles. The target journal's convention puts descriptive text in
#  the caption, so all titles/subtitles are omitted here.
#
#    (Fig6a) volcano            <- results/tables/TCGA_DEG_results.csv
#    (Fig6b) top-30 DEG heatmap <- results/rdata/tcga_processed.RData
#    (Fig7a) LASSO coefficients <- results/validation/signature_coefficients.csv
#
#  Values are NOT recomputed: each panel is drawn from the committed result
#  table, so the plotted numbers are identical to those already reported.
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ggrepel); library(ComplexHeatmap)
  library(circlize); library(grid)
})
ROOT <- Sys.getenv("GC_ROOT", "/nfsshare/users/P126156127/workspace/bioinf/gastric_cancer")
setwd(ROOT)
OUT <- "results/figures/clean"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_blank(), plot.subtitle = element_blank(),
        legend.position = "top", legend.title = element_blank())

## ---- (Fig6a) volcano, no title -------------------------------------------
res <- read.csv("results/tables/TCGA_DEG_results.csv")
stopifnot(nrow(res) == 21446)
nu <- sum(res$sig == "Up"); nd <- sum(res$sig == "Down")
stopifnot(nu == 2134, nd == 2362)          # assert against reported counts
res$lab <- factor(res$sig, levels = c("Up","Down","NS"),
                  labels = c(sprintf("Up in tumour (%s)", format(nu, big.mark=",")),
                             sprintf("Down in tumour (%s)", format(nd, big.mark=",")),
                             "Not significant"))
top <- res %>% filter(sig != "NS", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
  group_by(sig) %>% slice_min(padj, n = 10) %>% ungroup()
pv <- ggplot(res, aes(log2FoldChange, -log10(padj), colour = lab)) +
  geom_point(size = 0.35, alpha = 0.45, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey45", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey45", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = top, aes(label = hgnc_symbol), size = 2.5,
                           fontface = "italic", max.overlaps = Inf, box.padding = 0.32,
                           segment.size = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = setNames(c("#c0392b","#2471a3","grey78"), levels(res$lab))) +
  labs(x = expression(log[2]~"fold change (tumour/normal)"),
       y = expression(-log[10]~"(adjusted"~italic(P)*")")) + th
ggsave(file.path(OUT, "fig6a_volcano_clean.png"), pv, width = 6.2, height = 5.0,
       dpi = 300, bg = "white")

## ---- (Fig7a) LASSO-Cox coefficients, no title ----------------------------
co <- read.csv("results/validation/signature_coefficients.csv")
stopifnot(nrow(co) == 25)
co <- co[order(co$coefficient), ]
co$gene <- factor(co$gene, levels = co$gene)
co$dir  <- ifelse(co$coefficient > 0, "Higher risk", "Protective")
pc <- ggplot(co, aes(coefficient, gene, colour = dir)) +
  geom_segment(aes(x = 0, xend = coefficient, y = gene, yend = gene), linewidth = 0.5) +
  geom_point(size = 1.9) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  scale_colour_manual(values = c("Higher risk" = "#c0392b", "Protective" = "#2471a3")) +
  labs(x = "LASSO-Cox coefficient (log-hazard per SD)", y = NULL) +
  th + theme(axis.text.y = element_text(face = "italic", size = 7))
ggsave(file.path(OUT, "fig7a_coefficients_clean.png"), pc, width = 4.6, height = 5.0,
       dpi = 300, bg = "white")

## ---- (Fig6b) top-30 DEG heatmap, no title --------------------------------
e <- new.env(); load("results/rdata/tcga_processed.RData", envir = e)
vst <- e$tcga_vst; cd <- e$col_data
sel <- res %>% filter(sig != "NS", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
  group_by(sig) %>% slice_min(padj, n = 15) %>% ungroup()
## tcga_vst is keyed by HGNC symbol (not Ensembl id) -- match on symbol.
idx <- match(sel$hgnc_symbol, rownames(vst))
keep <- !is.na(idx)
mat <- vst[idx[keep], , drop = FALSE]
rownames(mat) <- sel$hgnc_symbol[keep]
grp <- factor(ifelse(cd$status == "Tumor", "Tumour", "Normal"), levels = c("Normal","Tumour"))
z   <- t(scale(t(mat)))
z[z >  2] <-  2; z[z < -2] <- -2
png(file.path(OUT, "fig6b_heatmap_clean.png"), width = 2100, height = 1750, res = 260)
ht <- Heatmap(z, name = "z-scored\nexpression",
  col = colorRamp2(c(-2,0,2), c("#2166AC","white","#B2182B")),
  column_split = grp, cluster_column_slices = FALSE,
  show_column_names = FALSE, row_names_gp = gpar(fontsize = 7, fontface = "italic"),
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_split = factor(ifelse(sel$sig[keep] == "Up", "tumour-up", "tumour-down"),
                     levels = c("tumour-up","tumour-down")),
  row_title_gp = gpar(fontsize = 8), heatmap_legend_param = list(labels_gp = gpar(fontsize = 7)))
draw(ht, merge_legend = TRUE)
dev.off()
message("wrote 3 clean panels -> ", OUT)

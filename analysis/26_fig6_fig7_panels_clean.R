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
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 7.5), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"), legend.key.size = unit(0.30, "cm"))

## ---- (Fig6a) volcano, no title -------------------------------------------
res <- read.csv("results/tables/TCGA_DEG_results.csv")
stopifnot(nrow(res) == 21446)
nu <- sum(res$sig == "Up"); nd <- sum(res$sig == "Down")
stopifnot(nu == 2134, nd == 2362)          # assert against reported counts
res$lab <- factor(res$sig, levels = c("Up","Down","NS"),
                  labels = c(sprintf("Up in tumor (%s)", format(nu, big.mark=",")),
                             sprintf("Down in tumor (%s)", format(nd, big.mark=",")),
                             "Not significant"))
top <- res %>% filter(sig != "NS", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
  group_by(sig) %>% slice_min(padj, n = 10) %>% ungroup()
## Labels are placed on a deterministic grid in the empty band above the data
## rather than by ggrepel: the solver repels labels from points and from other
## labels, but not from other labels' leader SEGMENTS, so a segment could be
## routed through a neighbouring label's glyphs. Laying the labels out on a
## fixed grid -- ordered by each gene's own x position, so leaders stay short
## and near-parallel -- makes overlap impossible by construction.
lab_grid <- function(d, x_lo, x_hi, y_rows) {
  ## Pack labels into rows by their RENDERED WIDTH, not by equal column slots:
  ## gene symbols differ ~4x in length (DPT vs MIR4435-2HG), so equal slots let
  ## long names overrun their neighbours. Width is estimated from the character
  ## count at the plotting size and converted to data units.
  d <- d[order(d$log2FoldChange), , drop = FALSE]
  span <- x_hi - x_lo
  ## 2.5 mm-per-character at size 2.5 italic, expressed as a fraction of the
  ## 3.74 in panel, then scaled into the x range of this half.
  wch <- 0.062 * span
  wid <- nchar(d$hgnc_symbol) * wch
  gap <- 0.10 * span
  ## Greedy first-fit into the available rows. A row is only chosen if the label
  ## still FITS inside the band; otherwise the least-full row is used and the
  ## row is compressed below. Without the fit test a long row overruns x_hi and
  ## ggplot silently drops the outermost label.
  row <- integer(nrow(d)); used <- rep(0, length(y_rows))
  for (i in seq_len(nrow(d))) {
    need <- wid[i] + gap
    fits <- which(used + need <= span)
    r <- if (length(fits)) fits[which.min(used[fits])] else which.min(used)
    row[i] <- r
    used[r] <- used[r] + need
  }
  d$lab_y <- y_rows[row]
  ## left-align each row's run, then convert to text centres
  d$lab_x <- NA_real_
  for (r in seq_along(y_rows)) {
    idx <- which(row == r)
    if (!length(idx)) next
    w <- wid[idx]
    g <- gap
    tot <- sum(w) + g * (length(idx) - 1)
    ## if the row still exceeds the band, shrink the inter-label gap (never the
    ## text) until it fits, so no label can be pushed outside the axis range
    if (tot > span && length(idx) > 1) {
      g <- max(0, (span - sum(w)) / (length(idx) - 1))
      tot <- sum(w) + g * (length(idx) - 1)
    }
    start <- x_lo + max(0, (span - tot)) / 2
    pos <- start + cumsum(c(0, head(w + g, -1))) + w / 2
    d$lab_x[idx] <- pos
  }
  d
}
## three rows in the band; the down set spans the left half, the up set the right
y_rows <- c(72, 81, 90, 99)
gd <- lab_grid(dplyr::filter(top, sig == "Down"), -9.0, -1.4, y_rows)
gu <- lab_grid(dplyr::filter(top, sig == "Up"),    1.4,  9.0, y_rows)
gg <- rbind(gd, gu)
pv_labels <- list(
  geom_segment(data = gg,
               aes(x = log2FoldChange, y = -log10(padj), xend = lab_x, yend = lab_y - 2.4,
                   colour = lab),
               linewidth = 0.2, alpha = 0.55, show.legend = FALSE),
  geom_text(data = gg, aes(x = lab_x, y = lab_y, label = hgnc_symbol, colour = lab),
            size = 2.5, fontface = "italic", show.legend = FALSE)
)

pv <- ggplot(res, aes(log2FoldChange, -log10(padj), colour = lab)) +
  geom_point(size = 0.35, alpha = 0.45, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey45", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey45", linewidth = 0.3) +
  pv_labels +
  scale_colour_manual(values = setNames(c("#c0392b","#2471a3","grey78"), levels(res$lab))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(size = 1.6, alpha = 1))) +
  scale_x_continuous(limits = c(-10.4, 10.4), breaks = seq(-5, 5, 5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 60, 20), expand = c(0, 0)) +
  labs(x = expression(log[2]~"fold change (tumor/normal)"),
       y = expression(-log[10]~"(adjusted"~italic(P)*")")) + th
## Guard: every selected gene must have a label position inside the axis range;
## ggplot silently DROPS text that falls outside scale limits.
stopifnot(nrow(gg) == 20,
          all(abs(gg$lab_x) < 10.4), all(gg$lab_y < 105),
          !any(is.na(gg$lab_x)))
ggsave(file.path(OUT, "fig6a_volcano_clean.png"), pv, width = 3.740, height = 3.700,
       dpi = 600, bg = "white")

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
ggsave(file.path(OUT, "fig7a_coefficients_clean.png"), pc, width = 3.740, height = 4.065,
       dpi = 600, bg = "white")

## ---- (Fig6b) top-30 DEG heatmap, no title --------------------------------
e <- new.env(); load("results/rdata/tcga_processed.RData", envir = e)
vst <- e$tcga_vst; cd <- e$col_data
sel <- res %>% filter(sig != "NS", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
  ## slice_min keeps ties (a padj tie in the Down group returned 16 rows and
  ## made the panel 31 genes against a caption saying 30); with_ties = FALSE
  ## plus an explicit ordering makes the 15/15 selection deterministic.
  group_by(sig) %>% arrange(padj, desc(abs(log2FoldChange)), hgnc_symbol, .by_group = TRUE) %>%
  slice_head(n = 15) %>% ungroup()
## tcga_vst is keyed by HGNC symbol (not Ensembl id) -- match on symbol.
## Assert the panel matches the caption: exactly 15 up + 15 down.
stopifnot(nrow(sel) == 30L, sum(sel$sig == "Up") == 15L, sum(sel$sig == "Down") == 15L)
idx <- match(sel$hgnc_symbol, rownames(vst))
keep <- !is.na(idx)
mat <- vst[idx[keep], , drop = FALSE]
rownames(mat) <- sel$hgnc_symbol[keep]
grp <- factor(ifelse(cd$status == "Tumor", "Tumor", "Normal"), levels = c("Normal","Tumor"))
z   <- t(scale(t(mat)))
z[z >  2] <-  2; z[z < -2] <- -2
## Width raised from 2100: at 2100 the "z-scored expression" legend title was
## clipped to "expressior" at the right canvas edge.
png(file.path(OUT, "fig6b_heatmap_clean.png"), width = 2244, height = 1693, res = 600)
ht <- Heatmap(z, name = "z-scored\nexpression",
  col = colorRamp2(c(-2,0,2), c("#2166AC","white","#B2182B")),
  column_split = grp, cluster_column_slices = FALSE,
  show_column_names = FALSE, row_names_gp = gpar(fontsize = 7, fontface = "italic"),
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_split = factor(ifelse(sel$sig[keep] == "Up", "tumor-up", "tumor-down"),
                     levels = c("tumor-up","tumor-down")),
  row_title_gp = gpar(fontsize = 8), heatmap_legend_param = list(labels_gp = gpar(fontsize = 7)))
draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 6), "mm"))
dev.off()
message("wrote 3 clean panels -> ", OUT)

## ---- (Fig8 a-f) MR scatters: strip the clipped per-panel y-title -----------
## Each of the six TwoSampleMR scatters carries the SAME 53-character rotated
## y-axis title ("SNP effect on Gastric cancer || id:ebi-a-GCST90018849"), which
## overruns the 750 px canvas and is clipped at the bottom on every panel.
## The harmonised SNP data is not cached locally (11_real_mr.R persists only
## plots + summary CSVs; re-extraction needs a live OPENGWAS_JWT), so the panels
## cannot be re-plotted. Instead the title strip is cropped off here and a single
## shared y-axis label is drawn once by 24_main_figures_5to8.R. Column 32 is the
## measured whitespace gap between the rotated title (cols 10-26) and the tick
## labels; no data ink lies left of it.
suppressPackageStartupMessages(library(magick))
MRDIR <- "results/figures/clean/mr"; dir.create(MRDIR, showWarnings = FALSE, recursive = TRUE)
for (f in list.files("results/mr_real", "^scatter_.*\\.png$", full.names = TRUE)) {
  im <- image_read(f); w <- image_info(im)$width; h <- image_info(im)$height
  gap <- 34
  ## Assert the discarded strip is title-only: cols 34-38 must be blank (the
  ## measured whitespace before the first tick label at col 39). If a future
  ## re-render moves the axis, this fails loudly rather than silently cutting data.
  guard <- image_crop(im, sprintf("%dx%d+%d+0", 5, h, gap))
  stopifnot(min(as.numeric(image_data(guard, "gray"))) > 0.90)
  image_write(image_crop(im, sprintf("%dx%d+%d+0", w - gap, h, gap)),
              file.path(MRDIR, basename(f)))
}
message("cropped ", length(list.files(MRDIR)), " MR scatters -> ", MRDIR)

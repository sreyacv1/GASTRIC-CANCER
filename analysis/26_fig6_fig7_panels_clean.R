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
## How many genes can be labelled without collisions is a geometric question:
## each label needs its own column at least as wide as the label itself. Solve
## for the largest n that fits the 7.6-unit band on BOTH sides, then take the
## top-n by adjusted p per direction.
PANEL_W_IN <- 6.90      # plotting area of the 3.74 in device, less axis furniture
X_SPAN     <- 20.8      # scale_x_continuous limits: -10.4 .. 10.4

## Rendered half-width of a label in DATA units, measured from the actual text
## grob (a per-character estimate overstates wide glyphs by ~40%).
lab_halfwidth <- function(sym) {
  ins <- vapply(sym, function(s)
    grid::convertWidth(grid::grobWidth(grid::textGrob(
      s, gp = grid::gpar(fontsize = 2.5 * ggplot2::.pt, fontface = "italic"))),
      "in", valueOnly = TRUE), numeric(1))
  unname(0.5 * ins * (X_SPAN / PANEL_W_IN))
}

## Lay one row out in x: start each label at its own point, separate to remove
## overlap, then shift the run (never individual labels) back inside the band.

## Largest label count that VERIFIES collision-free. Rather than predict
## capacity from a width formula (which mispredicts, because a wide label
## overhangs its column into neighbours), build the actual layout for each
## candidate N and run the same geometric test used in the final assertion.
lab_grid <- function(d, x_lo, x_hi, y_rows) {
  ## INTERLEAVED COLUMNS.
  ## Every label gets its own column across the band, and consecutive columns
  ## alternate between the available levels. A leader rises vertically in its
  ## own column: the only boxes it could meet are those in the SAME column, and
  ## each column holds exactly one label. Crossing is impossible by construction
  ## regardless of how many levels are used -- extra levels exist only so that a
  ## wide name may overhang its neighbours' columns without touching them.
  d <- d[order(d$log2FoldChange), , drop = FALSE]
  n <- nrow(d); nl <- length(y_rows)
  w <- lab_halfwidth(d$hgnc_symbol)
  d$lab_x <- x_lo + (x_hi - x_lo) * (seq_len(n) - 0.5) / n
  lev <- rep(seq_len(nl), length.out = n)      # interleave across levels
  d$lab_y <- y_rows[lev]
  d$lab_w <- w
  d$lab_lev <- lev
  d$corr_y <- min(y_rows) - 5.0
  ## constraint check: labels sharing a level must not overlap
  bad <- 0
  for (l in unique(lev)) {
    i <- which(lev == l); if (length(i) < 2) next
    o <- i[order(d$lab_x[i])]
    bad <- bad + sum(diff(d$lab_x[o]) - (head(w[o], -1) + tail(w[o], -1)) < 0)
  }
  attr(d, "overlaps") <- bad
  d
}

layout_clean <- function(d, x_lo, x_hi, y_rows, HH = 2.2) {
  g <- lab_grid(d, x_lo, x_hi, y_rows)
  if (attr(g, "overlaps") > 0) return(FALSE)
  for (a in seq_len(nrow(g))) for (b in seq_len(nrow(g))) {
    if (a == b) next
    xa <- g$lab_x[a]; y1 <- g$corr_y[a]; y2 <- g$lab_y[a] - 2.6
    if (xa > g$lab_x[b] - g$lab_w[b] && xa < g$lab_x[b] + g$lab_w[b] &&
        y2 > g$lab_y[b] - HH && y1 < g$lab_y[b] + HH) return(FALSE)
  }
  TRUE
}
n_fit <- function(res, y_rows) {
  for (N in 10:3) {
    d_up <- res %>% dplyr::filter(sig == "Up", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
      dplyr::slice_min(padj, n = N)
    d_dn <- res %>% dplyr::filter(sig == "Down", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
      dplyr::slice_min(padj, n = N)
    if (layout_clean(d_up, 1.3, 9.6, y_rows) && layout_clean(d_dn, -9.6, -1.3, y_rows))
      return(N)
  }
  3
}
Y_ROWS <- c(72, 81, 90, 99)
N_LAB <- n_fit(res, Y_ROWS)
message(sprintf("volcano: labelling top %d genes per direction (collision-free capacity)", N_LAB))
top <- res %>%
  dplyr::filter(sig != "NS", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
  dplyr::group_by(sig) %>% dplyr::slice_min(padj, n = N_LAB) %>% dplyr::ungroup()
## Labels are placed on a deterministic grid in the empty band above the data
## rather than by ggrepel: the solver repels labels from points and from other
## labels, but not from other labels' leader SEGMENTS, so a segment could be
## routed through a neighbouring label's glyphs. Laying the labels out on a
## fixed grid -- ordered by each gene's own x position, so leaders stay short
## and near-parallel -- makes overlap impossible by construction.
## Crossing-free label layout by CORRIDOR ROUTING.
##
## The hard constraint that defeats ggrepel -- and defeated a plain grid -- is
## that a leader travelling to an upper row can pass through a label sitting in
## a lower row. Rather than search for an assignment with no such pair (often
## infeasible: 10 labels per side do not leave enough clear x), each row is given
## its own horizontal corridor immediately below it. A leader runs obliquely to
## its row's corridor, along the corridor to its label's x, then straight up a
## short stub into the label. Corridors lie in the gaps BETWEEN rows, so a leader
## never enters a row it is not destined for, and within a row labels are laid
## out non-overlapping. Both constraints hold by construction.

row_layout <- function(x_pref, w, x_lo, x_hi) {
  o <- order(x_pref); x <- x_pref[o]; ww <- w[o]
  pad <- 0.04 * (x_hi - x_lo)
  for (i in seq_along(x)[-1]) {
    lim <- x[i - 1] + ww[i - 1] + ww[i] + pad
    if (x[i] < lim) x[i] <- lim
  }
  over <- (x[length(x)] + ww[length(ww)]) - x_hi
  if (over > 0) x <- x - over
  under <- x_lo - (x[1] - ww[1])
  if (under > 0) x <- x + under
  out <- numeric(length(x)); out[o] <- x; out
}


## four rows in the empty band; down set on the left half, up set on the right
y_rows <- Y_ROWS
gd <- lab_grid(dplyr::filter(top, sig == "Down"), -9.6, -1.3, y_rows)
gu <- lab_grid(dplyr::filter(top, sig == "Up"),    1.3,  9.6, y_rows)
gg <- rbind(gd, gu)

pv_labels <- list(
  ## oblique run from the data point to the foot of the label's own column
  geom_segment(data = gg, aes(x = log2FoldChange, y = -log10(padj),
                              xend = lab_x, yend = corr_y, colour = lab),
               linewidth = 0.2, alpha = 0.5, show.legend = FALSE),
  ## vertical rise inside that column -- no other label occupies this x
  geom_segment(data = gg, aes(x = lab_x, y = corr_y,
                              xend = lab_x, yend = lab_y - 2.6, colour = lab),
               linewidth = 0.2, alpha = 0.5, show.legend = FALSE),
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
stopifnot(attr(gd, "overlaps") == 0, attr(gu, "overlaps") == 0,
          nrow(gg) == 2 * N_LAB,
          all(abs(gg$lab_x) < 10.4), all(gg$lab_y < 105),
          !any(is.na(gg$lab_x)))
ggsave(file.path(OUT, "fig6a_volcano_clean.png"), pv, width = 7.400, height = 3.900,
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
## (b) shares the lower row with (c), so it is authored at half the figure
## width and taller: at the old 2244 x 1693 the 30 gene labels ran into the
## panel edge once the montage scaled it into the narrower slot.
png(file.path(OUT, "fig6b_heatmap_clean.png"), width = 2100, height = 2280, res = 600)
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

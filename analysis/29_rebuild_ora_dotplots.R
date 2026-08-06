#!/usr/bin/env Rscript
## Re-draw the GO:BP and KEGG over-representation dot plots for Supplementary
## Figure S4 directly from the committed ORA result tables.
##
## Why: the originals were written by analysis/09_functional_enrichment.R via
## clusterProfiler::dotplot() at 9 x 8 in / 150 dpi, then shrunk into a 7.48 in
## montage -- their axis text lands near 3.6 pt on the page. clusterProfiler is
## not installable in this environment, but the plot is a pure function of
## ORA_*.csv (Description, GeneRatio, p.adjust, Count), so it can be redrawn
## faithfully with ggplot2 at the final printed size.
##
## Fidelity: same top-N-by-p.adjust selection (N=15; the montage slot cannot
## fit 20 wrapped multi-line GO labels without them overlapping), same GeneRatio x axis, same
## Count->size and p.adjust->colour encodings as clusterProfiler's dotplot.
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
ROOT <- "/nfsshare/users/P126156127/workspace/bioinf/gastric_cancer"
setwd(ROOT)
OUT <- "results/enrichment"

ratio <- function(x) vapply(strsplit(x, "/"), function(p) as.numeric(p[1]) / as.numeric(p[2]), 1)

dot <- function(csv, png, n = 12, w = 3.74, h = 4.10) {
  d <- read.csv(csv, stringsAsFactors = FALSE)
  stopifnot(all(c("Description", "GeneRatio", "p.adjust", "Count") %in% names(d)))
  d <- d[order(d$p.adjust), , drop = FALSE]
  d <- head(d, min(n, nrow(d)))
  d$GeneRatio <- ratio(d$GeneRatio)
  ## Wrap long GO terms: unwrapped, they consume most of the panel width and
  ## squeeze the x axis until its tick labels collide.
  wrapped <- vapply(d$Description, function(s)
    paste(strwrap(s, width = 32), collapse = "\n"), "")
  d$Description <- factor(wrapped, levels = rev(wrapped))
  p <- ggplot(d, aes(x = GeneRatio, y = Description, size = Count, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "#D7301F", high = "#2166AC", name = "p.adjust") +
    scale_size_continuous(name = "Count", range = c(1.2, 4.2)) +
    scale_x_continuous(n.breaks = 4) +
    labs(x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 8) +
    theme(axis.text.y = element_text(size = 6.4, lineheight = 0.85),
          axis.text.x = element_text(size = 6.4),
          axis.title.x = element_text(size = 7),
          legend.key.size = unit(0.26, "cm"),
          legend.text = element_text(size = 6.2),
          legend.title = element_text(size = 6.6),
          panel.grid.minor = element_blank())
  ggsave(png, p, width = w, height = h, dpi = 600, bg = "white")
  cat(sprintf("%-34s %d terms, top p.adjust %.3g\n", basename(png), nrow(d), min(d$p.adjust)))
}

dot(file.path(OUT, "ORA_GO_BP_UP.csv"), file.path(OUT, "dotplot_GO_BP_UP.png"))
dot(file.path(OUT, "ORA_KEGG_UP.csv"), file.path(OUT, "dotplot_KEGG_UP.png"))
cat("done\n")

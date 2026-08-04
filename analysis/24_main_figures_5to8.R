#!/usr/bin/env Rscript
# ============================================================================
#  24_main_figures_5to8.R
#
#  Builds main Figures 5-8 as journal-spec TIFF + PNG.
#
#  These four were previously ORPHANS: montages of standalone PNGs with no
#  build script, so the repo could not regenerate them. Each is assembled here
#  from committed source images with explicit (a)/(b)/(c) panel tags matching
#  the lower-case convention used by the target journal.
#
#  Output: submission_package/figures_tiff/Fig{5..8}.tiff  (600 dpi, LZW)
#          results/figures/Fig{5..8}.png                   (300 dpi preview)
# ============================================================================
suppressPackageStartupMessages({
  library(magick); library(ggplot2); library(patchwork); library(cowplot)
})
ROOT <- Sys.getenv("GC_ROOT", "/nfsshare/users/P126156127/workspace/bioinf/gastric_cancer")
setwd(ROOT)
TIFDIR <- "submission_package/figures_tiff"; dir.create(TIFDIR, showWarnings = FALSE, recursive = TRUE)

panel_of <- function(path, tag) {
  stopifnot(file.exists(path))
  ggdraw() +
    draw_image(path) +
    draw_label(tag, x = 0.012, y = 0.985, hjust = 0, vjust = 1,
               fontface = "bold", size = 13)
}
# figure -> list(panels = c(file = tag), width_mm, aspect)
FIGS <- list(
  `5` = list(p = c("results/wgcna_real/wgcna_dendrogram.png"   = "(a)",
                   "results/wgcna_real/wgcna_module_trait.png"      = "(b)",
                   "results/wgcna_real/wgcna_power.png"         = "(c)"),
             ncol = 2, mm = 190),
  ## Fig6 uses an explicit layout: (a) and (b) share the top row, (c) is a wide
  ## strip beneath. A 1-column stack gives aspect 1.96, far outside the
  ## reference range (0.34-1.00) and unreadable at 190 mm.
  `6` = list(p = c("results/plots/transcriptome/deg_volcano.png"       = "(a)",
                   "results/plots/transcriptome/deg_heatmap_top30_clean.png" = "(b)",
                   "results/figures/deg_concordance_panel.png"         = "(c)"),
             ncol = 2, mm = 190, design = "AB\nCC", aspect = 0.62),
  `7` = list(p = c("results/validation/signature_coefficients_clean.png" = "(a)",
                   "results/validation_multi/forest_HR.png"               = "(b)",
                   "results/composite_figures/s16_km.png"           = "(c)"),
             ncol = 2, mm = 190),
  `8` = list(p = c("results/composite_figures/s18_mr_scatter_all.png" = ""),
             ncol = 1, mm = 190)
)
for (nm in names(FIGS)) {
  f <- FIGS[[nm]]
  present <- file.exists(names(f$p))
  if (!all(present)) {
    message("Fig", nm, ": MISSING -> ", paste(names(f$p)[!present], collapse = ", "))
    next
  }
  pl <- Map(panel_of, names(f$p), unname(f$p))
  comb <- if (!is.null(f$design)) wrap_plots(pl, design = f$design)
          else wrap_plots(pl, ncol = f$ncol)
  w_in <- f$mm / 25.4
  if (!is.null(f$aspect)) {
    h_in <- w_in * f$aspect
  } else {
    ars <- vapply(names(f$p), function(x) { d <- image_info(image_read(x)); d$height/d$width }, 1)
    nrow <- ceiling(length(pl) / f$ncol)
    h_in <- w_in * mean(ars) * nrow / f$ncol
  }
  ## guard: journal figures should not exceed a full page (aspect <= 1.35)
  if (h_in / w_in > 1.35) warning(sprintf("Fig%s aspect %.2f exceeds page proportion", nm, h_in/w_in))
  ggsave(file.path(TIFDIR, sprintf("Fig%s.tiff", nm)), comb, width = w_in, height = h_in,
         dpi = 600, bg = "white", compression = "lzw", limitsize = FALSE)
  ggsave(sprintf("results/figures/Fig%s.png", nm), comb, width = w_in, height = h_in,
         dpi = 300, bg = "white", limitsize = FALSE)
  message("wrote Fig", nm, " (", f$mm, "mm, ", length(pl), " panels)")
}

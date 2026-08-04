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

## Fig8's six panels shared one identical y-axis title that overran the source
## canvas and was clipped at the bottom. 26_fig6_fig7_panels_clean.R crops it
## off; it is drawn once here for the whole composed figure.
YLAB_FIG8 <- "SNP effect on gastric cancer (ebi-a-GCST90018849)"

TAG_PT <- 14   # measured to match plot.tag size = 14 of Figs 1-4 at 300 dpi
panel_of <- function(path, tag) {
  stopifnot(file.exists(path))
  ggdraw() +
    draw_image(path) +
    ## TAG_PT is tuned so the rendered tag height matches the patchwork
    ## plot.tag of Figs 1-4 exactly (both files render at 300 dpi).
    draw_label(tag, x = 0.012, y = 0.985, hjust = 0, vjust = 1,
               fontface = "bold", size = TAG_PT)
}
# figure -> list(panels = c(file = tag), width_mm, aspect)
FIGS <- list(
  ## Panels re-rendered title-free by analysis/25_fig5_panels_clean.R:
  ## the journal convention puts descriptive text in the caption, not on the plot.
  `5` = list(p = c("results/wgcna_real/clean/dendrogram_clean.png"   = "a",
                   "results/wgcna_real/clean/module_trait_clean.png" = "b",
                   "results/wgcna_real/clean/power_clean.png"        = "c"),
             ncol = 2, mm = 190, design = "AB\nCC", aspect = 0.72),
  ## Fig6 uses an explicit layout: (a) and (b) share the top row, (c) is a wide
  ## strip beneath. A 1-column stack gives aspect 1.96, far outside the
  ## reference range (0.34-1.00) and unreadable at 190 mm.
  `6` = list(p = c("results/figures/clean/fig6a_volcano_clean.png"     = "a",
                   "results/figures/clean/fig6b_heatmap_clean.png"           = "b",
                   "results/figures/deg_concordance_panel.png"         = "c"),
             ncol = 2, mm = 190, design = "AB\nCC", aspect = 0.62),
  ## (c)-(f): the four cohort KM curves are used INDIVIDUALLY. The old
  ## s16_km.png was a pre-composited montage carrying its own nested (A)-(D)
  ## labels, which collided with this figure's panel tags.
  `7` = list(p = c("results/figures/clean/fig7a_coefficients_clean.png"  = "a",
                   "results/validation_multi/forest_HR.png"              = "b",
                   "results/validation/KM_TCGA.png"                      = "c",
                   "results/validation/KM_ACRG.png"                      = "d",
                   "results/validation_multi/KM_GSE15459.png"            = "e",
                   "results/validation_multi/KM_GSE84437.png"            = "f"),
             ncol = 2, mm = 190, design = "AB\nCD\nEF", aspect = 1.0),
  ## Built from the six per-exposure scatters rather than the pre-composited
  ## s18_mr_scatter_all.png, which carried a baked-in overall title.
  `8` = list(p = c("results/figures/clean/mr/scatter_H__pylori_IgG_seropositivity.png" = "a",
                   "results/figures/clean/mr/scatter_Streptococcus__genus_.png"        = "b",
                   "results/figures/clean/mr/scatter_Fusobacterium.png"                = "c",
                   "results/figures/clean/mr/scatter_Prevotella.png"                   = "d",
                   "results/figures/clean/mr/scatter_Veillonella.png"                  = "e",
                   "results/figures/clean/mr/scatter_Lactobacillus.png"                = "f"),
             ncol = 3, mm = 190, design = "ABC\nDEF", aspect = 0.68)
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
  ## Fig8: the identical per-panel y-title was cropped off upstream (it was
  ## clipped by the source canvas); draw it once for the whole figure instead.
  if (nm == "8") {
    comb <- wrap_elements(comb) +
      labs(tag = YLAB_FIG8) +
      theme(plot.tag = element_text(size = 9, angle = 90, vjust = 1),
            plot.tag.position = "left")
  }
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

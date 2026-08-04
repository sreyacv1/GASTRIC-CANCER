#!/usr/bin/env Rscript
# ============================================================================
#  25_fig5_panels_clean.R
#
#  Re-renders the three Figure-5 panels WITHOUT in-panel titles.
#  The target journal's convention (see FIGURE_CAPTIONS.md) puts all
#  descriptive text in the caption, not on the plot. The original
#  14_wgcna_real.R panels carry main="..." titles, so they are re-plotted
#  here from the saved network object rather than cropped.
# ============================================================================
suppressPackageStartupMessages({ library(WGCNA) })
ROOT <- Sys.getenv("GC_ROOT", "/nfsshare/users/P126156127/workspace/bioinf/gastric_cancer")
setwd(ROOT)
OUT <- "results/wgcna_real/clean"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
e <- new.env(); load("results/wgcna_real/wgcna_real.RData", envir = e)

## ---- (a) dendrogram, no main title ----
png(file.path(OUT, "dendrogram_clean.png"), width = 2000, height = 1150, res = 200)
par(mar = c(1.2, 4.2, 0.6, 0.6))
plotDendroAndColors(e$net$dendrograms[[1]],
  e$moduleColors[e$net$blockGenes[[1]]],
  "Module", dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

## ---- (b) module-trait heatmap, no main title ----
png(file.path(OUT, "module_trait_clean.png"), width = 2000, height = 1500, res = 200)
par(mar = c(8.5, 8.5, 0.8, 1.2))
txt <- paste0(signif(e$mtCor, 2), "\n(", signif(e$mtP, 1), ")")
dim(txt) <- dim(e$mtCor)
labeledHeatmap(Matrix = e$mtCor, xLabels = colnames(e$mtCor),
  yLabels = rownames(e$mtCor), ySymbols = rownames(e$mtCor),
  colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt,
  setStdMargins = FALSE, cex.text = 0.42, cex.lab = 0.72, zlim = c(-1, 1),
  main = "")
dev.off()

## ---- (c) soft-power robustness, no per-panel titles ----
## The displayed panel is the power-ROBUSTNESS check (does the prognostic
## module survive a change of soft power?), not the raw scale-free fit.
## Source: results/wgcna_real/power_robustness_summary.csv
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
pr <- read.csv("results/wgcna_real/power_robustness_summary.csv")
th <- theme_classic(base_size = 9) +
      theme(plot.title = element_blank(), plot.subtitle = element_blank())
p1 <- ggplot(pr, aes(power, HR_perSD)) +
  geom_errorbar(aes(ymin = HR_lo, ymax = HR_hi), width = 0.35, colour = "#B2182B") +
  geom_point(size = 1.9, colour = "#B2182B") + geom_line(colour = "#B2182B") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  labs(x = "Soft-thresholding power", y = "Module HR per SD") + th
p2 <- ggplot(pr, aes(power, -log10(cox_p))) +
  geom_point(size = 1.9, colour = "#2166AC") + geom_line(colour = "#2166AC") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  labs(x = "Soft-thresholding power", y = expression(-log[10]~italic(p)~"(Cox)")) + th
p3 <- ggplot(pr, aes(power, hub8_comembership)) +
  geom_point(size = 1.9, colour = "#1B7837") + geom_line(colour = "#1B7837") +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(x = "Soft-thresholding power", y = "Hub co-membership fraction") + th
ggsave(file.path(OUT, "power_clean.png"), p1 | p2 | p3,
       width = 10, height = 3.0, dpi = 200, bg = "white")
message("wrote 3 clean Fig5 panels -> ", OUT)

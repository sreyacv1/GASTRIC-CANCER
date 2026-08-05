#!/usr/bin/env Rscript
## Committed build script for six supplementary figures that previously had NO
## build script (produced by ad-hoc inline code in earlier sessions):
##   S4  results/enrichment/path_ORA_GO_KEGG.png      (GO:BP + KEGG ORA montage)
##   S7  results/microbiome_biomarker/da_clr_barplot.png
##   S8  results/composite_figures/s15_immune.png
##   S9  results/composite_figures/s17_nomogram.png
##   S10 results/composite_figures/s19_mr_loo_all.png
##   S11 results/microbiome_biomarker/rf_batch_classifier.png
##
## Values are read from the committed result tables and never recomputed;
## assertions below fail loudly if a table changes under the caption.
suppressMessages({library(ggplot2); library(patchwork); library(magick); library(dplyr)})
## Ghostscript is absent on this host, so magick cannot rasterise PDF panels;
## image_read_pdf (pdftools backend) is used for the vector immune panels.
ROOT <- "."
rd <- function(p) read.csv(file.path(ROOT, p), check.names = FALSE)
ci <- function(p, dens = 300) {
  f <- file.path(ROOT, p); stopifnot(file.exists(f))
  if (grepl("\\.pdf$", f)) image_read_pdf(f, density = dens) else image_read(f)
}
pan <- function(img) ggplot() + annotation_raster(as.raster(img), -Inf, Inf, -Inf, Inf) +
  theme_void() + coord_fixed(ratio = image_info(img)$height / image_info(img)$width)

## ---- S4: ORA montage (GO:BP UP | KEGG UP) --------------------------------
s4 <- (pan(ci("results/enrichment/dotplot_GO_BP_UP.png")) |
       pan(ci("results/enrichment/dotplot_KEGG_UP.png"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave("results/enrichment/path_ORA_GO_KEGG.png", s4, width = 9, height = 4,
       dpi = 300, bg = "white")

## ---- S7: differential abundance, control vs cancer-adjacent (CLR) --------
da <- rd("results/microbiome_biomarker/04a_DA_control_vs_GCN.csv")
stopifnot(nrow(da) > 0, all(c("genus","effect_clr","ci_lo","ci_hi","q") %in% names(da)))
sig <- da %>% filter(q < 0.05) %>% arrange(desc(abs(effect_clr))) %>% head(30)
cat(sprintf("S7: %d genera q<0.05; plotting top %d by |CLR effect|\n",
            sum(da$q < 0.05), nrow(sig)))
p7 <- ggplot(sig, aes(reorder(genus, effect_clr), effect_clr,
                      fill = effect_clr > 0)) +
  geom_col() + geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25,
                             linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "#4C72B0", "TRUE" = "#C44E52"),
                    labels = c("depleted", "enriched"), name = NULL) +
  labs(x = NULL, y = "CLR effect size (cancer-adjacent vs control)") +
  theme_bw(base_size = 10)
ggsave("results/microbiome_biomarker/da_clr_barplot.png", p7, width = 8, height = 7,
       dpi = 300, bg = "white")

## ---- S8: immune montage --------------------------------------------------
s8 <- (pan(ci("results/plots/immune/Immune_Cell_Boxplots_Status.pdf")) |
       pan(ci("results/plots/immune/Immune_Pathway_Correlation_Heatmap.pdf"))) /
      pan(ci("results/plots/immune/EMT_vs_CD8_Correlation.pdf")) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave("results/composite_figures/s15_immune.png", s8, width = 9, height = 7.5,
       dpi = 300, bg = "white")

## ---- S9: nomogram + calibration -----------------------------------------
s9 <- (pan(ci("results/nomogram/Nomogram_GC_OS_real.png")) |
       pan(ci("results/nomogram/Calibration_curves_real.png"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave("results/composite_figures/s17_nomogram.png", s9, width = 10, height = 3.2,
       dpi = 300, bg = "white")

## ---- S10: all six MR leave-one-out panels -------------------------------
loo <- sort(Sys.glob(file.path(ROOT, "results/mr_real/loo_*.png")))
stopifnot(length(loo) == 6)
s10 <- wrap_plots(lapply(loo, function(f) pan(ci(sub("^\\./", "", f)))), ncol = 3) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave("results/composite_figures/s19_mr_loo_all.png", s10, width = 10, height = 6,
       dpi = 300, bg = "white")

## ---- S11: microbiome RF importance + batch sanity -----------------------
imp <- rd("results/microbiome_biomarker/05_rf_importance.csv") %>%
  arrange(desc(MeanDecreaseGini)) %>% head(10)
met <- rd("results/microbiome_biomarker/05_rf_metrics_and_batch_sanity.csv")
mv  <- setNames(met$value, met$metric)
stopifnot(abs(mv["cancer_vs_control_AUC"]      - 0.916119204847165) < 1e-9,
          abs(mv["flowcell_pred_accuracy"]     - 0.775891341256367) < 1e-9,
          abs(mv["flowcell_majority_baseline"] - 0.544991511035654) < 1e-9)
cat(sprintf("S11: cancer AUC %.3f | flowcell acc %.3f vs baseline %.3f\n",
            mv["cancer_vs_control_AUC"], mv["flowcell_pred_accuracy"],
            mv["flowcell_majority_baseline"]))
p11a <- ggplot(imp, aes(reorder(genus, MeanDecreaseGini), MeanDecreaseGini)) +
  geom_col(fill = "#4C72B0") + coord_flip() +
  labs(x = NULL, y = "Mean decrease in Gini") + theme_bw(base_size = 10)
bar <- data.frame(
  what = factor(c("Cancer vs control\n(AUC)", "Flowcell\n(accuracy)",
                  "Flowcell majority\nbaseline"),
                levels = c("Cancer vs control\n(AUC)", "Flowcell\n(accuracy)",
                           "Flowcell majority\nbaseline")),
  v = as.numeric(mv[c("cancer_vs_control_AUC", "flowcell_pred_accuracy",
                      "flowcell_majority_baseline")]))
p11b <- ggplot(bar, aes(what, v)) +
  geom_col(fill = c("#C44E52", "#DD8452", "grey65")) +
  geom_text(aes(label = sprintf("%.2f", v)), vjust = -0.4, size = 3) +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  labs(x = NULL, y = "Performance") + theme_bw(base_size = 10)
s11 <- (p11a | p11b) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave("results/microbiome_biomarker/rf_batch_classifier.png", s11, width = 9,
       height = 4, dpi = 300, bg = "white")

cat("wrote S4, S7, S8, S9, S10, S11 from committed sources\n")

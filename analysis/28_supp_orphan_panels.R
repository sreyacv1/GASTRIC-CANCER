#!/usr/bin/env Rscript
## Committed build script for six supplementary figures that previously had NO
## build script (they were produced by ad-hoc inline code in earlier sessions):
##   S4  results/enrichment/path_ORA_GO_KEGG.png       GO:BP + KEGG ORA montage
##   S7  results/microbiome_biomarker/da_clr_barplot.png
##   S8  results/composite_figures/s15_immune.png
##   S9  results/composite_figures/s17_nomogram.png
##   S10 results/composite_figures/s19_mr_loo_all.png
##   S11 results/microbiome_biomarker/rf_batch_classifier.png
##
## The panel layout and content here reproduce the FIGURES AS PUBLISHED and as
## described in their PAPER.md captions (verified panel-by-panel against the
## committed PNGs). Values come from the committed result tables and are never
## recomputed; the assertions fail loudly if a table changes under a caption.
##
## Ghostscript is absent on this host so magick cannot rasterise PDFs;
## image_read_pdf (pdftools backend) handles the vector panels.
suppressMessages({library(ggplot2); library(patchwork); library(magick)
                  library(dplyr); library(pdftools)})
ROOT <- "."
rd <- function(p) read.csv(file.path(ROOT, p), check.names = FALSE)
ci <- function(p, dens = 300) {
  f <- file.path(ROOT, p); stopifnot(file.exists(f))
  if (grepl("\\.pdf$", f)) image_read_pdf(f, density = dens) else image_read(f)
}
pan <- function(img) ggplot() +
  annotation_raster(as.raster(img), -Inf, Inf, -Inf, Inf) + theme_void() +
  coord_fixed(ratio = image_info(img)$height / image_info(img)$width)
tag <- function(p) p + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

## Some component panels were written by scripts whose deconvolution step needs
## MCPcounter (GitHub-only, not installable here), so their in-panel titles are
## removed at source in analysis/08_immune_deconvolution.R and
## analysis/17_external_utility_ACRG.R for future runs, and cropped here so the
## current montages already follow the no-in-panel-title convention. The crop
## fraction is measured per file by locating the contiguous ink bands at the top
## of each canvas: most have a single title band ending ~4-5% down, but the
## validation scatter has title AND subtitle bands (ending 3.9% and 7.1%), so it
## needs the deeper crop. Measure, do not guess, when swapping a component.
crop_title <- function(img, frac) {
  inf <- image_info(img)
  off <- round(inf$height * frac)
  image_crop(img, sprintf("%dx%d+0+%d", inf$width, inf$height - off, off))
}
cit <- function(p, frac, dens = 300) crop_title(ci(p, dens), frac)

## ---- S4: ORA montage, GO:BP UP | KEGG UP ---------------------------------
ggsave("results/enrichment/path_ORA_GO_KEGG.png",
       tag(pan(ci("results/enrichment/dotplot_GO_BP_UP.png")) |
           pan(ci("results/enrichment/dotplot_KEGG_UP.png"))),
       width = 9, height = 4, dpi = 300, bg = "white")

## ---- S7: two-panel CLR differential abundance ----------------------------
## Caption: left = control vs cancer-adjacent (44/61 at q<0.05); right = paired
## cancer-adjacent vs tumour (18/61); 12 most enriched + 12 most depleted shown.
da_a <- rd("results/microbiome_biomarker/04a_DA_control_vs_GCN.csv")
da_b <- rd("results/microbiome_biomarker/04b_DA_GCN_vs_GCT_paired.csv")
stopifnot(sum(da_a$q < 0.05) == 44L, sum(da_b$q < 0.05) == 18L,
          nrow(da_a) == 61L, nrow(da_b) == 61L)
cat(sprintf("S7: control-vs-GCN %d/61 sig | paired GCN-vs-GCT %d/61 sig\n",
            sum(da_a$q < 0.05), sum(da_b$q < 0.05)))
da_panel <- function(d, eff, ttl) {
  d$effect <- d[[eff]]; d$sig <- d$q < 0.05
  top <- bind_rows(d %>% arrange(desc(effect)) %>% head(12),
                   d %>% arrange(effect)       %>% head(12)) %>% distinct(genus, .keep_all = TRUE)
  ggplot(top, aes(reorder(genus, effect), effect,
                  fill = ifelse(!sig, "ns", ifelse(effect > 0, "up", "dn")))) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c(up = "#C44E52", dn = "#4C72B0", ns = "grey75"),
                      breaks = c("up", "dn", "ns"),
                      labels = c("enriched (q<0.05)", "depleted (q<0.05)", "n.s."),
                      name = NULL) +
    labs(x = NULL, y = ttl) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
}
ggsave("results/microbiome_biomarker/da_clr_barplot.png",
       (da_panel(da_a, "effect_clr", "CLR effect: control vs cancer-adjacent") |
        da_panel(da_b, "effect_clr_paired", "CLR effect: cancer-adjacent vs tumour (paired)")) +
         plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 12, height = 6.8, dpi = 300, bg = "white")

## LEGIBILITY, THE ONE RULE THAT MATTERS FOR MONTAGES
## A montage shrinks every component to its slot, so text that reads fine in the
## standalone PNG can become unreadable in the figure. Two failures were shipped
## before this note existed: a 7-facet panel whose strip and tick text vanished at
## montage scale, and an x-tick clipped mid-word to "Norma" because its facet was
## narrower than the label. Neither is detectable from the montage raster (see the
## note in analysis/check_paper_integrity.py), so it is enforced here:
##   1. Author each component at close to its FINAL printed size - a panel destined
##      for a half-width slot should be ~4-4.5 in wide, not 6-8 in shrunk down.
##   2. Keep facets per panel low (<=4). Seven facets in one quadrant cannot work.
##   3. After any change, VIEW the montage at full size and read every tick label.
##      Edge-ink and aspect-ratio checks pass on an illegible figure.

## ---- S8 panel (b): rebuilt from the committed score table ---------------
## The original results/plots/Immune_tumor_vs_normal.png packed 7 facets into a
## 1800x900 canvas; at montage scale its strip/tick text became illegible and the
## Macrophages facet's x-tick was clipped to "Norma". Rebuilt here from
## results/immune/deconvolution_scores.csv (long form) with a 2-row layout,
## readable text and no clipped labels. Wilcoxon adjusted p values are read from
## results/immune/tumor_vs_normal_stats.csv, never recomputed.
sc <- rd("results/immune/deconvolution_scores.csv")
st <- rd("results/immune/tumor_vs_normal_stats.csv")
rownames(sc) <- sc$feature; sc$feature <- NULL
## TCGA barcode position 14-15: "01"-"09" tumour, "10"-"19" normal.
tcode <- substr(sub("^([^-]+-[^-]+-[^-]+)-.*$", "\\4", colnames(sc)), 1, 2)
tcode <- substr(sapply(strsplit(colnames(sc), "-"), `[`, 4), 1, 2)
grp <- ifelse(as.integer(tcode) <= 9, "Tumour", "Normal")
stopifnot(sum(grp == "Tumour") > 300, sum(grp == "Normal") > 20)
cat(sprintf("S8b: %d tumour / %d normal samples\n",
            sum(grp == "Tumour"), sum(grp == "Normal")))
## Map the seven populations in the stats table to their score rows.
pop_map <- c("CD8 T cells (MCP)"          = "MCP_CD8 T cells",
             "T cells (MCP)"              = "MCP_T cells",
             "Cytotoxic lymphocytes (MCP)" = "MCP_Cytotoxic lymphocytes",
             "Monocytic lineage (MCP)"    = "MCP_Monocytic lineage",
             "CD8+ T-cells (xCell)"       = "xCell_CD8+ T-cells",
             "Macrophages (xCell)"        = "xCell_Macrophages",
             "Monocytes (xCell)"          = "xCell_Monocytes")
stopifnot(all(pop_map %in% rownames(sc)), all(names(pop_map) %in% st$population))
long <- do.call(rbind, lapply(names(pop_map), function(pp) {
  data.frame(population = pp, group = grp,
             score = as.numeric(sc[pop_map[[pp]], ]), row.names = NULL)
}))
lab <- setNames(sprintf("%s\np_adj = %.2g", st$population, st$p_adj_BH), st$population)
long$facet <- factor(lab[long$population], levels = lab[names(pop_map)])
p8b <- ggplot(long, aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.35, linewidth = 0.3) +
  facet_wrap(~ facet, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c(Normal = "#4C72B0", Tumour = "#C44E52")) +
  labs(x = NULL, y = "Deconvolution score") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(size = 7.2, lineheight = 1.05),
        axis.text.x = element_text(size = 7.5))
ggsave("results/plots/Immune_tumor_vs_normal_clean.png", p8b,
       width = 8.2, height = 4.2, dpi = 300, bg = "white")

## ---- S8: four-panel immune montage --------------------------------------
## (a) validation scatter (b) tumour vs normal (c) by subtype (d) CD8 KM
ggsave("results/composite_figures/s15_immune.png",
       tag((pan(cit("results/plots/Immune_validation_scatter.png", 0.080)) |
            pan(ci("results/plots/Immune_tumor_vs_normal_clean.png"))) /
           (pan(cit("results/plots/Immune_by_subtype.png",          0.058)) |
            pan(cit("results/plots/Immune_CD8_survival_KM.png",     0.050)))),
       width = 9, height = 7.1, dpi = 300, bg = "white")

## ---- S9: three-panel nomogram / calibration / external DCA --------------
## Layout note: a 3-across row gives each panel ~2.4 in at 180 mm, which shrank
## the calibration and DCA tick text below the legibility screen's reference
## value. The nomogram is wide and short, so it takes the full width and the two
## square panels sit below at ~3.5 in each.
## Row heights follow the content aspects, not patchwork's equal default: the
## nomogram is 2400x1500 (h/w 0.62) so at 9 in wide it needs ~5.6 in, while the
## lower row holds two ~0.8-aspect panels at 4.5 in wide each, needing ~3.6 in.
## Equal rows stretched the nomogram's canvas and left a large empty band under it.
ggsave("results/composite_figures/s17_nomogram.png",
       tag(pan(ci("results/nomogram_combined/combined_nomogram.png")) /
           (pan(cit("results/nomogram_combined/calibration_combined.png", 0.078)) |
            pan(cit("results/external_utility_ACRG/DCA_external.png", 0.050))) +
           plot_layout(heights = c(0.62 * 9, 3.6))),
       width = 9, height = 9.4, dpi = 300, bg = "white")

## ---- S10: all six MR leave-one-out panels -------------------------------
loo <- sort(Sys.glob(file.path(ROOT, "results/mr_real/loo_*.png")))
stopifnot(length(loo) == 6)
## 2 columns, not 3: the source loo_*.png are 900x750 with text sized for a
## standalone panel, so a 3-across row rendered their axis labels below the
## legibility screen's reference value. They cannot be re-rendered at larger text
## (the harmonised MR data was never persisted and needs a live OpenGWAS token),
## so the layout carries the fix.
ggsave("results/composite_figures/s19_mr_loo_all.png",
       tag(wrap_plots(lapply(loo, function(f) pan(ci(sub("^\\./", "", f)))), ncol = 2)),
       width = 9, height = 10.2, dpi = 300, bg = "white")

## ---- S11: RF importance (15 genera, contaminants red) + batch sanity ----
imp <- rd("results/microbiome_biomarker/05_rf_importance.csv") %>%
  arrange(desc(MeanDecreaseGini)) %>% head(15)
met <- rd("results/microbiome_biomarker/05_rf_metrics_and_batch_sanity.csv")
mv  <- setNames(met$value, met$metric)
stopifnot(abs(mv["cancer_vs_control_AUC"]      - 0.916119204847165) < 1e-9,
          abs(mv["flowcell_pred_accuracy"]     - 0.775891341256367) < 1e-9,
          abs(mv["flowcell_majority_baseline"] - 0.544991511035654) < 1e-9)
## Contaminant genera named in the S11 caption; all six must be in the top 15.
contam <- c("Dietzia", "Serinicoccus", "Methylobacterium-Methylorubrum",
            "Microbacterium", "Sphingomonas", "Serratia")
stopifnot(all(contam %in% imp$genus))
cat(sprintf("S11: cancer AUC %.3f | flowcell %.3f vs baseline %.3f | %d/6 contaminants in top 15\n",
            mv["cancer_vs_control_AUC"], mv["flowcell_pred_accuracy"],
            mv["flowcell_majority_baseline"], sum(contam %in% imp$genus)))
imp$is_contam <- imp$genus %in% contam
p11a <- ggplot(imp, aes(reorder(genus, MeanDecreaseGini), MeanDecreaseGini,
                        fill = is_contam)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey55", "TRUE" = "#C44E52"),
                    labels = c("gastric/oral taxon", "known contaminant"), name = NULL) +
  labs(x = NULL, y = "Mean decrease in Gini") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(colour = ifelse(
          imp$genus[order(imp$MeanDecreaseGini)] %in% contam, "#C44E52", "grey20"),
          face = ifelse(imp$genus[order(imp$MeanDecreaseGini)] %in% contam,
                        "italic", "plain")))
bar <- data.frame(
  what = factor(c("Cancer vs control\n(AUC)", "Sequencing flowcell\n(accuracy)"),
                levels = c("Cancer vs control\n(AUC)", "Sequencing flowcell\n(accuracy)")),
  v = as.numeric(mv[c("cancer_vs_control_AUC", "flowcell_pred_accuracy")]))
p11b <- ggplot(bar, aes(what, v)) +
  geom_col(fill = c("#4C72B0", "#DD8452"), width = 0.6) +
  geom_hline(yintercept = mv["flowcell_majority_baseline"], linetype = 3) +
  annotate("text", x = 2.42, y = mv["flowcell_majority_baseline"], vjust = -0.5,
           hjust = 1, size = 2.6, label = "majority-class baseline 0.55") +
  geom_text(aes(label = sprintf("%.2f", v)), vjust = -0.5, size = 3) +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(x = NULL, y = "Performance") + theme_bw(base_size = 9)
ggsave("results/microbiome_biomarker/rf_batch_classifier.png",
       tag(p11a | p11b), width = 9.5, height = 4.3, dpi = 300, bg = "white")

cat("wrote S4, S7, S8, S9, S10, S11 from committed sources\n")

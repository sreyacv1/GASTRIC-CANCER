#!/usr/bin/env Rscript
## Re-render Supplementary Figure S2 (WGCNA stromal-module hub genes by cell
## type) from the cached Seurat object.
##
## Two defects in the committed PNG:
##   1. 1050x1200 px at 150 dpi -- the lowest-resolution figure in the paper.
##   2. It carries an in-panel title, which every other figure had removed to
##      match the journal convention (descriptive text belongs in the caption).
## Re-rendering from the cached object (not re-clustering) leaves the underlying
## cell assignments untouched, so the gene/cell-type content is unchanged.
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
ROOT <- "/nfsshare/users/P126156127/workspace/bioinf/gastric_cancer"
setwd(ROOT)
outdir <- "results/scrna"
rds <- file.path(outdir, "GSE134520_seurat.rds")
stopifnot(file.exists(rds))
merged <- readRDS(rds)

hub <- read.csv("results/wgcna_real/hub_genes_prognostic_module.csv",
                stringsAsFactors = FALSE)
gcol <- intersect(c("gene", "Gene", "hub_gene"), names(hub))[1]
hub_genes <- hub[[gcol]]
## Same panel definition as analysis/15_scrna_validation.R L126-129: the 13
## canonical stromal markers plus the top 25 hub genes by hubScore. Using the
## full 263-gene table would plot a DIFFERENT figure from the published one.
key_stromal <- c("POSTN","FAP","COL1A2","CDH11","SPARC","LUM","BGN","COL1A1",
                 "COL3A1","DCN","THBS2","FN1","VCAN")
present_hub <- intersect(unique(c(key_stromal, head(hub_genes, 25))), rownames(merged))
cat("hub genes present:", length(present_hub), "of", length(hub_genes), "\n")

lvl <- names(sort(table(merged$cell_type), decreasing = TRUE))
merged$cell_type <- factor(merged$cell_type, levels = lvl)

p <- DotPlot(merged, features = present_hub, group.by = "cell_type") +
  coord_flip() + RotatedAxis() +
  ggtitle(NULL) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_text(size = 7, face = "italic"),
        axis.text.x = element_text(size = 7.5),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.34, "cm"))
ggsave(file.path(outdir, "DotPlot_stromal_module_hub.png"), p,
       width = 5.6, height = 6.4, dpi = 600, bg = "white")
cat("wrote DotPlot_stromal_module_hub.png at 5.6 x 6.4 in, 600 dpi\n")
cat("cells:", ncol(merged), "| cell types:", nlevels(merged$cell_type), "\n")

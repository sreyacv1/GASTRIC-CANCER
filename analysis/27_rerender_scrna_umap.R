#!/usr/bin/env Rscript
## Re-render Supplementary Figure S1 (cell-type UMAP) from the cached Seurat
## object written by analysis/15_scrna_validation.R (results/scrna/GSE134520_seurat.rds).
##
## Why a separate script: 15_scrna_validation.R re-runs PCA/clustering/UMAP from
## the raw counts, and Seurat's clustering and UMAP embedding are stochastic, so
## re-running it to change only a plot label risks altering cluster assignments
## and cell-type proportions that the caption quotes. Reading the cached object
## re-renders the identical embedding. Plot settings here are kept byte-identical
## to the ggsave call in 15_scrna_validation.R.
suppressMessages({library(Seurat); library(ggplot2); library(scales)})
rds <- "results/scrna/GSE134520_seurat.rds"
stopifnot(file.exists(rds))
m <- readRDS(rds)
stopifnot("cell_type" %in% colnames(m@meta.data))

## Assertions on the quantities the S1 caption states.
stopifnot(ncol(m) == 43992L)
prop <- round(100 * table(m$cell_type) / ncol(m), 1)
expected <- c(Epithelial = 68.6, Endothelial = 7.6, Myeloid = 6.4, Plasma = 5.8,
              Tcell = 5.8, Fibroblast = 4.2, Mast = 1.0, Bcell = 0.6)
stopifnot(identical(sort(names(prop)), sort(names(expected))))
stopifnot(all(abs(prop[names(expected)] - expected) < 0.05))
cat("cells:", ncol(m), "| cell types:", length(unique(m$cell_type)),
    "| proportions match caption\n")

## No in-panel title: the dataset and colouring are stated in the caption.
## set.seed is load-bearing: repel = TRUE places the cell-type labels with
## ggrepel, whose nudging is stochastic, so unseeded renders differ by ~0.7% of
## pixels (label positions only). Seeding makes the figure byte-reproducible.
set.seed(42)
ctcols <- setNames(hue_pal()(length(unique(m$cell_type))), sort(unique(m$cell_type)))
p_umap <- DimPlot(m, group.by = "cell_type", label = TRUE, repel = TRUE,
                  cols = ctcols) + ggtitle(NULL)
ggsave("results/scrna/UMAP_celltypes.png", p_umap, width = 8, height = 6, dpi = 300,
       bg = "white")
cat("wrote results/scrna/UMAP_celltypes.png (300 dpi, no title)\n")

#!/usr/bin/env Rscript
# 15_scrna_validation.R
# Single-cell validation of the GC prognostic signature + WGCNA stromal module.
# Dataset: GSE134520 (Zhang et al. 2019, Cell Reports) — human gastric mucosa
# single-cell atlas spanning the premalignant -> early gastric cancer cascade
# (NAG, CAG, IMW/IMS intestinal metaplasia, EGC). 13 samples, processed UMI
# count matrices (genes x cells) from GEO supplementary files. REAL DATA ONLY.
#
# Goal: localise the 25-gene prognostic signature and the WGCNA "red" stromal/
# EMT module hub genes to cell types, testing whether the stromal module
# (POSTN, FAP, COL1A2, CDH11, SPARC ...) is CAF/fibroblast-restricted.

suppressPackageStartupMessages({
  library(data.table); library(Seurat); library(Matrix)
  library(ggplot2); library(dplyr)
})
set.seed(42)

root   <- Sys.getenv("GC_ROOT", getwd())
datdir <- file.path(root, "data/scrna/GSE134520")
outdir <- file.path(root, "results/scrna")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- sample -> disease-stage metadata (from GEO sample titles) -------------
files <- list.files(datdir, pattern = "_processed_.*\\.txt\\.gz$", full.names = TRUE)
sample_id <- sub(".*_processed_(.*)\\.txt\\.gz$", "\\1", files)
# stage: strip trailing replicate digits
stage <- gsub("[0-9]+$", "", sample_id)
stage_lab <- c(NAG="Non-atrophic gastritis", CAG="Chronic atrophic gastritis",
               IMW="Intestinal metaplasia (mild)", IMS="Intestinal metaplasia (severe)",
               EGC="Early gastric cancer")[stage]
cat("Found", length(files), "samples:\n")
print(data.frame(sample_id, stage, stage_lab))

# ---- read each matrix, build per-sample Seurat object ----------------------
obj_list <- list()
for (i in seq_along(files)) {
  # header line = cell barcodes (no corner label); data rows = gene<TAB>counts.
  # read header separately to avoid the off-by-one column misalignment.
  bc <- strsplit(readLines(gzfile(files[i]), n = 1), "\t")[[1]]
  bc <- bc[nzchar(bc)]
  m  <- fread(files[i], skip = 1, header = FALSE, data.table = FALSE)
  g  <- m[[1]]; m[[1]] <- NULL
  m  <- as.matrix(m)
  stopifnot(ncol(m) == length(bc))            # genes x cells alignment check
  if (anyDuplicated(g)) m <- rowsum(m, group = g) else rownames(m) <- g
  colnames(m) <- paste0(sample_id[i], "_", make.unique(bc))
  so <- CreateSeuratObject(counts = as(m, "dgCMatrix"), project = sample_id[i],
                           min.cells = 3, min.features = 200)
  so$sample <- sample_id[i]; so$stage <- stage[i]
  so$stage_lab <- unname(stage_lab[i])
  obj_list[[sample_id[i]]] <- so
  cat(sprintf("  %-6s %5d genes x %5d cells (post min.filter)\n",
              sample_id[i], nrow(so), ncol(so)))
}

merged <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = NULL)
merged <- JoinLayers(merged)
cat("Merged:", nrow(merged), "genes x", ncol(merged), "cells\n")

# ---- QC --------------------------------------------------------------------
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
before <- ncol(merged)
merged <- subset(merged, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                                  percent.mt < 20)
cat("QC: kept", ncol(merged), "of", before, "cells\n")

# ---- normalise / HVG / PCA / cluster / UMAP --------------------------------
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)
cat("Clusters:", nlevels(merged$seurat_clusters), "\n")

# ---- cell-type annotation by canonical markers -----------------------------
markers <- list(
  Epithelial   = c("EPCAM","KRT8","KRT18","KRT19","MUC1"),
  Fibroblast   = c("DCN","LUM","COL1A1","COL1A2","PDGFRB","FAP","COL3A1"),
  Tcell        = c("CD3D","CD3E","CD2","CD8A","IL7R"),
  Bcell        = c("CD79A","CD79B","MS4A1"),
  Plasma       = c("MZB1","IGHG1","JCHAIN"),
  Myeloid      = c("CD68","LYZ","CD14","AIF1","C1QA"),
  Endothelial  = c("PECAM1","VWF","CLDN5"),
  Mast         = c("TPSAB1","TPSB2","CPA3"))
markers <- lapply(markers, function(g) intersect(g, rownames(merged)))

# score each cluster by mean scaled expression of each marker set, assign max
avg <- AverageExpression(merged, features = unlist(markers),
                         group.by = "seurat_clusters", layer = "data")$RNA
zscore <- t(scale(t(log1p(avg))))                 # gene-wise z across clusters
setscore <- sapply(markers, function(gs) colMeans(zscore[gs, , drop = FALSE]))
rownames(setscore) <- sub("^g", "", rownames(setscore))  # AverageExpression prepends 'g'
celltype_by_cluster <- colnames(setscore)[max.col(setscore, ties.method = "first")]
names(celltype_by_cluster) <- rownames(setscore)  # cluster ids
cat("\nCluster -> cell type:\n"); print(celltype_by_cluster)

merged$cell_type <- unname(celltype_by_cluster[as.character(merged$seurat_clusters)])

# ---- outputs: UMAP ---------------------------------------------------------
ctcols <- setNames(scales::hue_pal()(length(unique(merged$cell_type))),
                   sort(unique(merged$cell_type)))
p_umap <- DimPlot(merged, group.by = "cell_type", label = TRUE, repel = TRUE,
                  cols = ctcols) +
  ggtitle("GSE134520 gastric mucosa scRNA-seq — cell types")
ggsave(file.path(outdir, "UMAP_celltypes.png"), p_umap, width = 8, height = 6, dpi = 150)

p_umap_stage <- DimPlot(merged, group.by = "stage") +
  ggtitle("GSE134520 — disease stage")
ggsave(file.path(outdir, "UMAP_stage.png"), p_umap_stage, width = 8, height = 6, dpi = 150)

# ---- gene panels -----------------------------------------------------------
hub <- read.csv(file.path(root, "results/wgcna_real/hub_genes_prognostic_module.csv"))
sig <- read.csv(file.path(root, "results/validation/signature_coefficients.csv"))
hub_genes <- hub$gene
sig_genes <- sig$gene

key_stromal <- c("POSTN","FAP","COL1A2","CDH11","SPARC","LUM","BGN","COL1A1",
                 "COL3A1","DCN","THBS2","FN1","VCAN")

present_hub <- intersect(unique(c(key_stromal, head(hub_genes, 25))), rownames(merged))
present_sig <- intersect(sig_genes, rownames(merged))
cat("\nHub genes present:", length(present_hub), "/", "\n")
cat("Signature genes present:", length(present_sig), "of", length(sig_genes), "\n")

lvl <- names(sort(table(merged$cell_type), decreasing = TRUE))
merged$cell_type <- factor(merged$cell_type, levels = lvl)

p_hub <- DotPlot(merged, features = present_hub, group.by = "cell_type") +
  coord_flip() + RotatedAxis() +
  ggtitle("WGCNA stromal module hub genes by cell type") +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(outdir, "DotPlot_stromal_module_hub.png"), p_hub,
       width = 7, height = 8, dpi = 150)

p_sig <- DotPlot(merged, features = present_sig, group.by = "cell_type") +
  coord_flip() + RotatedAxis() +
  ggtitle("25-gene prognostic signature by cell type") +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(outdir, "DotPlot_signature_genes.png"), p_sig,
       width = 7, height = 8, dpi = 150)

# ---- tables ----------------------------------------------------------------
# cell-type composition
comp <- as.data.frame(table(merged$cell_type)); colnames(comp) <- c("cell_type","n_cells")
comp$pct <- round(100 * comp$n_cells / sum(comp$n_cells), 2)
write.csv(comp, file.path(outdir, "celltype_composition.csv"), row.names = FALSE)
cat("\nCell-type composition:\n"); print(comp)

# per-cluster annotation table with scores
ann <- data.frame(cluster = rownames(setscore),
                  cell_type = celltype_by_cluster,
                  n_cells = as.integer(table(merged$seurat_clusters)[rownames(setscore)]),
                  round(setscore, 3))
write.csv(ann, file.path(outdir, "cluster_annotation.csv"), row.names = FALSE)

# per-gene dominant cell type: mean expression across cell types -> argmax
dominant_for <- function(genes) {
  ae <- AverageExpression(merged, features = genes, group.by = "cell_type",
                          layer = "data")$RNA
  ae <- as.matrix(ae)
  dom <- colnames(ae)[max.col(ae, ties.method = "first")]
  frac <- ae[cbind(seq_len(nrow(ae)), max.col(ae, ties.method = "first"))] / rowSums(ae)
  data.frame(gene = rownames(ae), dominant_cell_type = dom,
             frac_in_dominant = round(frac, 3), row.names = NULL)
}
dom_hub <- dominant_for(present_hub); dom_hub$panel <- "stromal_hub"
dom_sig <- dominant_for(present_sig); dom_sig$panel <- "signature"
dom_all <- rbind(dom_hub, dom_sig)
write.csv(dom_all, file.path(outdir, "gene_dominant_celltype.csv"), row.names = FALSE)
cat("\nStromal hub genes dominant cell type:\n"); print(dom_hub)
cat("\nSignature genes dominant cell type:\n"); print(dom_sig)

# ---- CAF localisation verdict ----------------------------------------------
caf_check <- dom_hub[dom_hub$gene %in% key_stromal, ]
n_fib <- sum(caf_check$dominant_cell_type == "Fibroblast")
cat(sprintf("\n=== VERDICT: %d/%d key stromal genes (POSTN,FAP,COL1A2,CDH11,SPARC,...) dominant in Fibroblasts ===\n",
            n_fib, nrow(caf_check)))

saveRDS(merged, file.path(outdir, "GSE134520_seurat.rds"))
cat("\nDone. Outputs in", outdir, "\n")

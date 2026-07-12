#!/usr/bin/env Rscript
# 09_functional_enrichment.R
# Real ORA + GSEA for gastric cancer DEGs. No fabricated data.
# ponytail: single script, sequential; splitting into modules buys nothing here.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(msigdbr)
  library(fgsea)
  library(DOSE)
  library(ggplot2)
})

set.seed(42)
outdir <- "results/enrichment"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
options(stringsAsFactors = FALSE)

log <- function(...) cat(sprintf(...), "\n")

# ---- helpers ---------------------------------------------------------------
sym2entrez <- function(symbols) {
  mp <- suppressWarnings(
    AnnotationDbi::select(org.Hs.eg.db, keys = unique(symbols),
                          keytype = "SYMBOL", columns = "ENTREZID"))
  mp <- mp[!is.na(mp$ENTREZID), ]
  mp[!duplicated(mp$SYMBOL), ]
}

save_enrich <- function(obj, path) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    log("  (no significant terms) %s", path); return(invisible())
  }
  write.csv(as.data.frame(obj), path, row.names = FALSE)
  log("  wrote %s (%d terms)", path, nrow(as.data.frame(obj)))
}

dotplot_safe <- function(obj, title, path, n = 20) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(invisible())
  p <- dotplot(obj, showCategory = min(n, nrow(as.data.frame(obj)))) +
    ggtitle(title)
  ggsave(path, p, width = 9, height = 8, dpi = 150)
  log("  plot %s", path)
}

# ============================================================================
# PART 1: ORA  (Tumor vs Normal)
# ============================================================================
log("== Loading TCGA Tumor-vs-Normal DEGs ==")
tn <- read.csv("results/tables/TCGA_DEG_results_symbols.csv")
tn <- tn[!is.na(tn$padj) & !is.na(tn$log2FoldChange) & tn$gene_symbol != "", ]

up_sym   <- unique(tn$gene_symbol[tn$padj < 0.05 & tn$log2FoldChange >  1])
down_sym <- unique(tn$gene_symbol[tn$padj < 0.05 & tn$log2FoldChange < -1])
log("UP genes: %d | DOWN genes: %d", length(up_sym), length(down_sym))

# universe = all tested genes mapped to entrez
uni_map <- sym2entrez(tn$gene_symbol)
universe <- unique(uni_map$ENTREZID)
up_ent   <- uni_map$ENTREZID[match(up_sym,   uni_map$SYMBOL)]
up_ent   <- up_ent[!is.na(up_ent)]
down_ent <- uni_map$ENTREZID[match(down_sym, uni_map$SYMBOL)]
down_ent <- down_ent[!is.na(down_ent)]
log("Mapped -> UP entrez: %d | DOWN entrez: %d | universe: %d",
    length(up_ent), length(down_ent), length(universe))

run_ora <- function(ent, tag) {
  for (ont in c("BP", "MF", "CC")) {
    ego <- tryCatch(enrichGO(gene = ent, universe = universe,
                             OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                             ont = ont, pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                             readable = TRUE),
                    error = function(e) { log("  enrichGO %s err: %s", ont, e$message); NULL })
    save_enrich(ego, file.path(outdir, sprintf("ORA_GO_%s_%s.csv", ont, tag)))
    dotplot_safe(ego, sprintf("GO:%s enrichment (%s)", ont, tag),
                 file.path(outdir, sprintf("dotplot_GO_%s_%s.png", ont, tag)))
  }
  ekg <- tryCatch(enrichKEGG(gene = ent, universe = universe,
                             organism = "hsa", pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                  error = function(e) { log("  enrichKEGG err: %s", e$message); NULL })
  if (!is.null(ekg) && nrow(as.data.frame(ekg)) > 0)
    ekg <- setReadable(ekg, org.Hs.eg.db, keyType = "ENTREZID")
  save_enrich(ekg, file.path(outdir, sprintf("ORA_KEGG_%s.csv", tag)))
  dotplot_safe(ekg, sprintf("KEGG enrichment (%s)", tag),
               file.path(outdir, sprintf("dotplot_KEGG_%s.png", tag)))
}

log("== ORA: UP genes ==");   run_ora(up_ent,   "UP")
log("== ORA: DOWN genes =="); run_ora(down_ent, "DOWN")

# ============================================================================
# PART 2: GSEA (Tumor vs Normal) — rank by DESeq2 stat
# ============================================================================
log("== GSEA Tumor-vs-Normal ==")
tn$entrez <- uni_map$ENTREZID[match(tn$gene_symbol, uni_map$SYMBOL)]
tn_r <- tn[!is.na(tn$entrez) & !is.na(tn$stat), ]
tn_r <- tn_r[!duplicated(tn_r$entrez), ]
ranks <- sort(setNames(tn_r$stat, tn_r$entrez), decreasing = TRUE)
log("ranked genes: %d", length(ranks))

# Hallmark (entrez = ncbi_gene in msigdbr >=10)
H <- msigdbr(species = "Homo sapiens", collection = "H")
H_list <- split(H$ncbi_gene, H$gs_name)
# C2:CP (canonical pathways: KEGG legacy + Reactome)
C2  <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
C2b <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
C2all <- rbind(C2, C2b)
C2_list <- split(C2all$ncbi_gene, C2all$gs_name)

run_fgsea <- function(pw, ranks, path, minSize = 10, maxSize = 500) {
  res <- fgsea(pathways = pw, stats = ranks, minSize = minSize,
               maxSize = maxSize, eps = 0.0)
  res <- res[order(res$padj), ]
  df <- as.data.frame(res)
  df$leadingEdge <- vapply(df$leadingEdge, function(x) paste(head(x, 30), collapse = ";"), "")
  write.csv(df, path, row.names = FALSE)
  log("  wrote %s (%d sets, %d sig padj<0.05)", path, nrow(df), sum(df$padj < 0.05, na.rm = TRUE))
  res
}

hall_tn <- run_fgsea(H_list,  ranks, file.path(outdir, "GSEA_Hallmark_TumorVsNormal.csv"))
c2_tn   <- run_fgsea(C2_list, ranks, file.path(outdir, "GSEA_C2CP_TumorVsNormal.csv"))

# GSEA enrichment plots for top hallmarks by |NES| among significant
sigH <- hall_tn[hall_tn$padj < 0.05, ]
sigH <- sigH[order(-abs(sigH$NES)), ]
topH <- head(sigH$pathway, 4)
for (i in seq_along(topH)) {
  p <- plotEnrichment(H_list[[topH[i]]], ranks) +
    ggtitle(sprintf("%s\nNES=%.2f  padj=%.1e", topH[i],
                    sigH$NES[i], sigH$padj[i]))
  ggsave(file.path(outdir, sprintf("GSEA_enrich_Hallmark_TN_%d_%s.png", i,
                                   gsub("[^A-Za-z0-9]+", "_", topH[i]))),
         p, width = 7, height = 5, dpi = 150)
}

# Hallmark NES barplot (all significant)
nes_barplot <- function(res, title, path, top = 25) {
  d <- as.data.frame(res)
  d <- d[!is.na(d$padj) & d$padj < 0.05, ]
  if (nrow(d) == 0) { log("  no sig hallmarks for barplot"); return(invisible()) }
  d <- d[order(-abs(d$NES)), ]
  d <- head(d, top)
  d$pathway <- gsub("HALLMARK_", "", d$pathway)
  d$dir <- ifelse(d$NES > 0, "Up in Tumor/Diffuse", "Down")
  p <- ggplot(d, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#2980b9"),
                      guide = "none") +
    labs(x = NULL, y = "NES", title = title) + theme_bw(base_size = 11)
  ggsave(path, p, width = 8, height = 7, dpi = 150)
  log("  plot %s", path)
}
nes_barplot(hall_tn, "Hallmark GSEA — Tumor vs Normal",
            file.path(outdir, "GSEA_Hallmark_NES_barplot_TumorVsNormal.png"))

# ============================================================================
# PART 3: GSEA Hallmark — Lauren Diffuse vs Intestinal
# ============================================================================
log("== GSEA Lauren Diffuse-vs-Intestinal ==")
lr <- read.csv("results/tables/Lauren_DEG_symbols.csv")
lr <- lr[!is.na(lr$t) & lr$gene_symbol != "", ]
lr_map <- sym2entrez(lr$gene_symbol)
lr$entrez <- lr_map$ENTREZID[match(lr$gene_symbol, lr_map$SYMBOL)]
lr_r <- lr[!is.na(lr$entrez), ]
lr_r <- lr_r[!duplicated(lr_r$entrez), ]
# positive t = up in Diffuse (per sig labels Up_Diffuse)
lranks <- sort(setNames(lr_r$t, lr_r$entrez), decreasing = TRUE)
log("Lauren ranked genes: %d", length(lranks))

hall_lr <- run_fgsea(H_list, lranks,
                     file.path(outdir, "GSEA_Hallmark_DiffuseVsIntestinal.csv"))
nes_barplot(hall_lr, "Hallmark GSEA — Diffuse vs Intestinal (NES>0 = up in Diffuse)",
            file.path(outdir, "GSEA_Hallmark_NES_barplot_DiffuseVsIntestinal.png"))
sigL <- hall_lr[hall_lr$padj < 0.05, ]
sigL <- sigL[order(-abs(sigL$NES)), ]
topL <- head(sigL$pathway, 4)
for (i in seq_along(topL)) {
  p <- plotEnrichment(H_list[[topL[i]]], lranks) +
    ggtitle(sprintf("%s\nNES=%.2f  padj=%.1e", topL[i], sigL$NES[i], sigL$padj[i]))
  ggsave(file.path(outdir, sprintf("GSEA_enrich_Hallmark_Lauren_%d_%s.png", i,
                                   gsub("[^A-Za-z0-9]+", "_", topL[i]))),
         p, width = 7, height = 5, dpi = 150)
}

# ============================================================================
# SUMMARY to console
# ============================================================================
topn <- function(res, n = 10, pos = TRUE) {
  d <- as.data.frame(res)
  d <- d[!is.na(d$padj) & d$padj < 0.05, ]
  d <- d[order(if (pos) -d$NES else d$NES), ]
  head(d[, c("pathway", "NES", "padj", "size")], n)
}
log("\n===== SUMMARY =====")
log("\n-- Top Hallmarks UP in Tumor --"); print(topn(hall_tn, 10, TRUE))
log("\n-- Top Hallmarks DOWN in Tumor --"); print(topn(hall_tn, 10, FALSE))
log("\n-- Top C2:CP UP in Tumor --"); print(topn(c2_tn, 10, TRUE))
log("\n-- Top Hallmarks UP in Diffuse --"); print(topn(hall_lr, 10, TRUE))
log("\nDONE")

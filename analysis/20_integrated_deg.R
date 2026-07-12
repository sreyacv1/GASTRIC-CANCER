#!/usr/bin/env Rscript
# =============================================================================
#  INTEGRATED tumour-vs-normal DEG on the harmonised TCGA + GTEx + GEO matrix.
#
#  Rationale: leverage GTEx's 407 normal-stomach samples as a robust normal
#  reference alongside TCGA/GEO. Discovery = TCGA + GTEx + GSE27342 + GSE63089;
#  GSE62254 is held out as a validation cohort (excluded here).
#
#  NOTE ON SCALE: combined_expr_bc is ComBat-corrected AND z-scored/scale-
#  compressed, so |logFC| is tiny and NOT meaningful. Significance is defined
#  by adj.P.Val < 0.05 and genes are ranked by the moderated t-statistic.
#
#  Model: ~ sample_type + dataset  (tumour effect adjusted for cohort).
#  Honest caveat (printed below): tumour samples live mainly in TCGA/GEO and
#  GTEx is all-normal, so batch and biology are partially confounded. The
#  dataset covariate + the fact that TCGA, GSE27342 and GSE63089 each contain
#  BOTH tumour and normal within-cohort mitigate but do not eliminate this.
# =============================================================================
suppressPackageStartupMessages({
  library(limma)
})

OUT     <- "results/tables"
ENR     <- "results/enrichment_integrated"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(ENR, showWarnings = FALSE, recursive = TRUE)

load("data/processed/combined_transcriptome.RData")  # combined_expr_bc, combined_meta
stopifnot(all(colnames(combined_expr_bc) == rownames(combined_meta)))
meta <- combined_meta
expr <- combined_expr_bc
cat(sprintf("Combined matrix: %d genes x %d samples\n", nrow(expr), ncol(expr)))

## ---------------------------------------------------------------------------
## 0. FIX GSE27342 label bug: meta tags all 160 as Tumor, but it is 80T/80N.
##    Real labels: phenotype title containing "control" -> Normal, else Tumor.
## ---------------------------------------------------------------------------
ph <- read.csv("data/host/GSE27342_phenotype.csv", stringsAsFactors = FALSE)
ph$gsm      <- ph$X
ph$is_norm  <- grepl("control", ph$title, ignore.case = TRUE)
lab_map     <- setNames(ifelse(ph$is_norm, "Normal", "Tumor"), ph$gsm)

is27 <- meta$dataset == "GSE27342"
mapped <- lab_map[rownames(meta)[is27]]
stopifnot(!any(is.na(mapped)))                 # every GSE27342 sample matched
meta$sample_type_simple[is27] <- unname(mapped)
cat(sprintf("GSE27342 labels fixed: %d Normal / %d Tumor (was 0/160)\n",
            sum(meta$sample_type_simple[is27]=="Normal"),
            sum(meta$sample_type_simple[is27]=="Tumor")))

cat("\nDataset x sample_type (after fix):\n")
print(table(meta$dataset, meta$sample_type_simple))

## ---------------------------------------------------------------------------
## 1. INTEGRATED DEG: discovery = TCGA + GTEx ONLY.
##    Tumour = TCGA 412; Normal = TCGA 36 + GTEx 407 (=443). GTEx augments the
##    normal baseline. Model ~ sample_type + dataset (dataset = TCGA vs GTEx);
##    tumour effect is identified from TCGA's within-cohort T-vs-N, adjusted for
##    the TCGA/GTEx mean shift. GEO cohorts are VALIDATION, not discovery.
## ---------------------------------------------------------------------------
disc <- meta$dataset %in% c("TCGA", "GTEx")
e <- expr[, disc]
m <- meta[disc, ]
sample_type <- factor(m$sample_type_simple, levels = c("Normal", "Tumor"))
dataset     <- factor(m$dataset)

cat(sprintf("\nIntegrated discovery (TCGA+GTEx): %d tumour / %d normal | %d genes\n",
            sum(sample_type=="Tumor"), sum(sample_type=="Normal"), nrow(e)))
cat("Discovery datasets:\n"); print(table(m$dataset, m$sample_type_simple))

design <- model.matrix(~ sample_type + dataset)
fit    <- eBayes(lmFit(e, design))
res    <- topTable(fit, coef = "sample_typeTumor", number = Inf, sort.by = "t")
res$gene <- rownames(res)
# Significance by adj.P only (scale is compressed -> logFC not usable).
# Direction from the moderated t-statistic.
res$sig <- ifelse(res$adj.P.Val < 0.05 & res$t > 0, "Up",
           ifelse(res$adj.P.Val < 0.05 & res$t < 0, "Down", "NS"))
res <- res[order(-res$t), ]
write.csv(res[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","sig")],
          file.path(OUT, "DEG_integrated_TCGA_GTEx.csv"), row.names = FALSE)

n_up <- sum(res$sig=="Up"); n_dn <- sum(res$sig=="Down")
cat(sprintf("\nIntegrated DEG (adj.P<0.05): %d up / %d down of %d genes\n",
            n_up, n_dn, nrow(res)))
cat("\nTop 20 UP (by t):\n")
print(head(res[res$t>0, c("gene","t","adj.P.Val")], 20), row.names=FALSE)
cat("\nTop 20 DOWN (by t):\n")
print(head(res[order(res$t), c("gene","t","adj.P.Val")], 20), row.names=FALSE)

## ---------------------------------------------------------------------------
## 2. WITHIN-COHORT GEO DEGs: GSE27342 (80/80) and GSE63089 (45/45).
## ---------------------------------------------------------------------------
geo_deg <- function(ds) {
  sel <- meta$dataset == ds
  if (sum(sel)==0 || length(unique(meta$sample_type_simple[sel]))<2) return(NULL)
  eg <- expr[, sel]
  gg <- factor(meta$sample_type_simple[sel], levels=c("Normal","Tumor"))
  fg <- eBayes(lmFit(eg, model.matrix(~ gg)))
  rg <- topTable(fg, coef="ggTumor", number=Inf, sort.by="t"); rg$gene <- rownames(rg)
  rg$sig <- ifelse(rg$adj.P.Val<0.05 & rg$t>0, "Up",
             ifelse(rg$adj.P.Val<0.05 & rg$t<0, "Down","NS"))
  write.csv(rg[,c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","sig")],
            file.path(OUT, sprintf("DEG_GEO_%s.csv", ds)), row.names=FALSE)
  cat(sprintf("\nGSE %s DEG (%d T / %d N): %d up / %d down (adj.P<0.05)\n",
              ds, sum(gg=="Tumor"), sum(gg=="Normal"),
              sum(rg$sig=="Up"), sum(rg$sig=="Down")))
  cat("Top 10 UP:\n"); print(head(rg[rg$t>0,c("gene","t","adj.P.Val")],10), row.names=FALSE)
  rg
}
rg27 <- geo_deg("GSE27342")
rg63 <- geo_deg("GSE63089")

## ---------------------------------------------------------------------------
## 3. CROSS-COHORT CONCORDANCE.
##    Integrated moderated-t  vs  TCGA-only log2FC, and vs each GEO cohort t.
## ---------------------------------------------------------------------------
conc <- function(name, gA, sA, gB, sB) {
  cg <- intersect(gA, gB)
  a  <- sA[match(cg, gA)]; b <- sB[match(cg, gB)]
  ok <- is.finite(a) & is.finite(b)
  r  <- cor(a[ok], b[ok], method="pearson")
  cat(sprintf("Concordance %-32s Pearson r = %+.3f (n=%d)\n", name, r, sum(ok)))
  data.frame(comparison=name, pearson_r=r, n=sum(ok))
}
cat("\n--- Cross-cohort concordance ---\n")
ctab <- list()
tcga <- read.csv("results/tables/TCGA_DEG_results_symbols.csv", stringsAsFactors=FALSE)
ctab[[1]] <- conc("integrated_t vs TCGA_log2FC", res$gene, res$t,
                  tcga$gene_symbol, tcga$log2FoldChange)
ctab[[2]] <- conc("integrated_t vs GSE27342_t", res$gene, res$t, rg27$gene, rg27$t)
ctab[[3]] <- conc("integrated_t vs GSE63089_t", res$gene, res$t, rg63$gene, rg63$t)
ctab[[4]] <- conc("GSE27342_t vs GSE63089_t",   rg27$gene, rg27$t, rg63$gene, rg63$t)
concord <- do.call(rbind, ctab)
write.csv(concord, file.path(OUT,"DEG_integrated_concordance.csv"), row.names=FALSE)

## ---------------------------------------------------------------------------
## 4. FUNCTIONAL ENRICHMENT on the INTEGRATED ranking.
##    (a) fgsea vs MSigDB Hallmark, ranked by integrated t.
##    (b) clusterProfiler GO-BP + KEGG ORA on adj.P<0.05 up / down sets.
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(fgsea); library(msigdbr)
  library(clusterProfiler); library(org.Hs.eg.db)
})
set.seed(1)

# ---- (a) fgsea Hallmark ----
# msigdbr API changed across versions; try new then old signature.
hm <- tryCatch(msigdbr(species="Homo sapiens", collection="H"),
               error=function(e) msigdbr(species="Homo sapiens", category="H"))
gs_col <- if ("gs_name" %in% names(hm)) "gs_name" else "gs_collection"
pathways <- split(hm$gene_symbol, hm[[gs_col]])
ranks <- setNames(res$t, res$gene)
ranks <- ranks[is.finite(ranks)]
ranks <- sort(ranks, decreasing=TRUE)
fg <- fgseaMultilevel(pathways, ranks, minSize=10, maxSize=500, eps=0)
fg <- fg[order(fg$padj, -abs(fg$NES)), ]
fg_out <- as.data.frame(fg)
fg_out$leadingEdge <- vapply(fg_out$leadingEdge, function(x) paste(head(x,20),collapse=";"), "")
write.csv(fg_out, file.path(ENR,"fgsea_Hallmark_integrated.csv"), row.names=FALSE)
cat(sprintf("\nfgsea Hallmark: %d/%d pathways at padj<0.05\n",
            sum(fg$padj<0.05, na.rm=TRUE), nrow(fg)))
cat("Top 15 Hallmark by padj (NES sign = up in tumour if +):\n")
print(head(fg_out[,c("pathway","NES","pval","padj","size")],15), row.names=FALSE)

# ---- (b) ORA: GO-BP + KEGG on up / down sets ----
up_g <- res$gene[res$sig=="Up"]; dn_g <- res$gene[res$sig=="Down"]
univ <- res$gene
sym2ent <- function(s) {
  m <- suppressWarnings(bitr(s, "SYMBOL","ENTREZID", org.Hs.eg.db))
  unique(m$ENTREZID)
}
uE <- sym2ent(up_g); dE <- sym2ent(dn_g); bgE <- sym2ent(univ)

run_ora <- function(entrez, bg, tag) {
  go <- tryCatch(enrichGO(entrez, org.Hs.eg.db, ont="BP", universe=bg,
                          pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1,
                          readable=TRUE), error=function(e) NULL)
  if (!is.null(go) && nrow(as.data.frame(go))>0) {
    write.csv(as.data.frame(go), file.path(ENR,sprintf("ORA_GO_BP_%s.csv",tag)), row.names=FALSE)
    cat(sprintf("\nGO-BP %s: %d terms padj<0.05. Top 8:\n", tag, sum(as.data.frame(go)$p.adjust<0.05)))
    print(head(as.data.frame(go)[,c("ID","Description","p.adjust","Count")],8), row.names=FALSE)
  } else cat(sprintf("\nGO-BP %s: no enriched terms\n", tag))
  kg <- tryCatch(enrichKEGG(entrez, organism="hsa", universe=bg,
                            pAdjustMethod="BH", pvalueCutoff=0.05),
                 error=function(e) NULL)
  if (!is.null(kg) && nrow(as.data.frame(kg))>0) {
    write.csv(as.data.frame(kg), file.path(ENR,sprintf("ORA_KEGG_%s.csv",tag)), row.names=FALSE)
    cat(sprintf("KEGG %s: %d pathways padj<0.05. Top 8:\n", tag, sum(as.data.frame(kg)$p.adjust<0.05)))
    print(head(as.data.frame(kg)[,c("ID","Description","p.adjust","Count")],8), row.names=FALSE)
  } else cat(sprintf("KEGG %s: no enriched pathways\n", tag))
}
run_ora(uE, bgE, "UP")
run_ora(dE, bgE, "DOWN")

cat("\n================ SAVED ================\n")
cat("results/tables/DEG_integrated_TCGA_GTEx.csv   [DISCOVERY: TCGA+GTEx]\n")
cat("results/tables/DEG_GEO_GSE27342.csv           [VALIDATION only]\n")
cat("results/tables/DEG_GEO_GSE63089.csv           [VALIDATION only]\n")
cat("results/tables/DEG_integrated_concordance.csv\n")
cat("results/enrichment_integrated/fgsea_Hallmark_integrated.csv\n")
cat("results/enrichment_integrated/ORA_{GO_BP,KEGG}_{UP,DOWN}.csv\n")
cat("\nCAVEAT: discovery integrates TCGA + GTEx only. All tumours are TCGA and\n")
cat("all GTEx samples are normal, so batch and biology are partially confounded.\n")
cat("~ sample_type + dataset identifies the tumour effect from TCGA's own\n")
cat("within-cohort tumour-vs-normal (36 TCGA normals present), with GTEx\n")
cat("augmenting the normal baseline; this mitigates but does not eliminate the\n")
cat("confound. GEO cohorts (GSE27342, GSE63089) are external VALIDATION only.\n")

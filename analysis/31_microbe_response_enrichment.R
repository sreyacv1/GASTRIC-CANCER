#!/usr/bin/env Rscript
# =============================================================================
#  31_microbe_response_enrichment.R  --  Route #2 of the integration bridge.
#  Test whether MICROBE-RESPONSE host pathways (the specific programmes gastric
#  dysbiosis / H. pylori are known to trigger) are enriched in the TCGA
#  tumour-vs-normal transcriptome. Turns "inflammation is up" into
#  "the pathways H. pylori specifically triggers are up".
# =============================================================================
suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
set.seed(1105)
OUT <- "results/integration"; dir.create(OUT, showWarnings=FALSE)
load("results/rdata/tcga_processed.RData")   # res: Ensembl rownames, stat

# --- ranked gene list (Wald stat), Ensembl -> Entrez -------------------------
ens <- sub("\\..*$", "", rownames(res))
map <- suppressMessages(mapIds(org.Hs.eg.db, keys=ens, column="ENTREZID",
                               keytype="ENSEMBL", multiVals="first"))
gl <- res$stat; names(gl) <- map[ens]
gl <- gl[!is.na(names(gl)) & is.finite(gl)]
gl <- gl[!duplicated(names(gl))]
gl <- sort(gl, decreasing=TRUE)
cat(sprintf("[RANK] %d genes ranked by Wald stat\n", length(gl)))

# --- KEGG GSEA: pull the H. pylori pathway + bacterial/immune pathways --------
kk <- gseKEGG(gl, organism="hsa", pvalueCutoff=1, verbose=FALSE, seed=TRUE)
kres <- as.data.frame(kk)
micro_kegg <- kres[grepl("Helicobacter|bacter|NF-kappa|Toll|NOD|IL-17|TNF|Cytokine|Chemokine|JAK|infection",
                         kres$Description, ignore.case=TRUE),
                   c("ID","Description","setSize","NES","pvalue","p.adjust")]
write.csv(micro_kegg, file.path(OUT,"microbe_response_KEGG.csv"), row.names=FALSE)
cat("\n[KEGG] microbe/inflammation pathways (tumour-vs-normal GSEA):\n")
print(micro_kegg[order(-micro_kegg$NES),], row.names=FALSE)
hp <- kres[kres$ID=="hsa05120",]
if(nrow(hp)) cat(sprintf("\n>>> KEGG 'Epithelial cell signaling in H. pylori infection' (hsa05120): NES=%.2f, p.adj=%.2g\n", hp$NES, hp$p.adjust))

# --- GO BP GSEA: response-to-bacterium / LPS / interferon --------------------
gg <- gseGO(gl, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=1, verbose=FALSE, seed=TRUE,
            minGSSize=10, maxGSSize=800)
gres <- as.data.frame(gg)
micro_go <- gres[grepl("response to bacter|response to lipopolysacc|response to molecule of bacterial|defense response to bacter|response to interferon|NF-kappaB|response to tumor necrosis|inflammatory response|response to interleukin-6",
                       gres$Description, ignore.case=TRUE),
                 c("ID","Description","setSize","NES","pvalue","p.adjust")]
micro_go <- micro_go[order(-micro_go$NES),]
write.csv(micro_go, file.path(OUT,"microbe_response_GO.csv"), row.names=FALSE)
cat("\n[GO BP] microbe/inflammation response terms (top 15 by NES):\n")
print(head(micro_go,15), row.names=FALSE)
cat(sprintf("\n[SUMMARY] microbe-response GO terms enriched UP (NES>0, p.adj<0.05): %d\n",
            sum(micro_go$NES>0 & micro_go$p.adjust<0.05, na.rm=TRUE)))
cat("[DONE] Route #2 microbe-response enrichment ->", OUT, "\n")

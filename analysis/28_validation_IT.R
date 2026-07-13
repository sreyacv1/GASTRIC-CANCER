#!/usr/bin/env Rscript
# =============================================================================
#  28_validation_IT.R -- INDEPENDENT cross-cohort validation of the gastric
#  dysbiosis biomarker in PRJNA641258 (Italy; Bologna/Rimini, IJMS 2020).
#
#  Why this cohort: gastric biopsy, V3-V4 (341F, same region as our PRJDB20660),
#  20 tumour (10 ADC + 10 SRCC) vs 20 matched controls, PAIRED per patient ->
#  tumour/normal sequenced together -> structurally CANNOT be flowcell-confounded
#  (the failure mode that invalidated PRJDB20660's tumour-vs-normal contrast).
#
#  Test: does the dysbiosis signal (diversity shift, oral-taxa "oralization",
#  Helicobacter) replicate HERE, in a batch-clean cohort? Cross-cohort transfer
#  cannot be batch (batch does not transfer across independent labs) -> a
#  positive result validates the biomarker as biology.
# =============================================================================
suppressMessages({ library(dada2); library(phyloseq); library(vegan) })
set.seed(1105)
RAW <- "data/microbiome/validation_IT/raw"
OUT <- "results/microbiome_biomarker/validation_IT"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
FILT<- "data/microbiome/validation_IT/filtered"; dir.create(FILT, showWarnings=FALSE)
man <- read.csv("data/microbiome/validation_IT/manifest.csv", stringsAsFactors=FALSE)
man$group <- ifelse(grepl("CTRL", man$disease), "Control", "Tumour")
THREADS <- 16

fnFs <- file.path(RAW, paste0(man$run, "_1.fastq.gz"))
fnRs <- file.path(RAW, paste0(man$run, "_2.fastq.gz"))
ok <- file.exists(fnFs) & file.exists(fnRs)
man <- man[ok,]; fnFs <- fnFs[ok]; fnRs <- fnRs[ok]
cat(sprintf("[IT] %d samples: %s\n", nrow(man), paste(names(table(man$group)),table(man$group),collapse=" ")))

filtFs <- file.path(FILT, paste0(man$run,"_F.fq.gz")); filtRs <- file.path(FILT, paste0(man$run,"_R.fq.gz"))
# V3-V4 341F/805R -> trimLeft primers; truncLen for 2x300 MiSeq
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(17,21), truncLen=c(270,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=THREADS)
keep <- file.exists(filtFs); filtFs<-filtFs[keep]; filtRs<-filtRs[keep]; man<-man[keep,]
errF <- learnErrors(filtFs, multithread=THREADS, verbose=0)
errR <- learnErrors(filtRs, multithread=THREADS, verbose=0)
ddF <- dada(filtFs, err=errF, multithread=THREADS, verbose=0)
ddR <- dada(filtRs, err=errR, multithread=THREADS, verbose=0)
CKPT <- file.path(OUT, "dada2_IT.rds")
if (file.exists(CKPT)) {
  ck <- readRDS(CKPT); seqtab <- ck$seqtab; tax <- ck$tax; man <- ck$man
  cat("[IT] loaded DADA2 checkpoint\n")
} else {
  mg  <- mergePairs(ddF, filtFs, ddR, filtRs, minOverlap=12, verbose=TRUE)
  seqtab <- removeBimeraDenovo(makeSequenceTable(mg), method="consensus", multithread=THREADS, verbose=TRUE)
  rownames(seqtab) <- man$run
  cat(sprintf("[IT] %d samples x %d ASVs\n", nrow(seqtab), ncol(seqtab)))
  tax <- assignTaxonomy(seqtab, "data/microbiome/ref/silva_nr99_v138.1_train_set.fa.gz",
                        multithread=THREADS, tryRC=TRUE)
  saveRDS(list(seqtab=seqtab, tax=tax, man=man), CKPT)
}

# phyloseq, clean host/off-target, genus
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax),
               sample_data(data.frame(man, row.names=man$run)))
tt <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors=FALSE)
bad <- is.na(tt$Kingdom) | tt$Kingdom != "Bacteria" |
       (!is.na(tt$Family) & tt$Family == "Mitochondria") |
       (!is.na(tt$Order)  & tt$Order  == "Chloroplast")
ps <- prune_taxa(rownames(tt)[!bad], ps); ps <- prune_samples(sample_sums(ps)>=1000, ps)
psG <- tax_glom(ps, "Genus", NArm=TRUE)
gn <- make.unique(as.character(tax_table(psG)[,"Genus"])); taxa_names(psG) <- gn
cat(sprintf("[IT] final %d samples x %d genera\n", nsamples(psG), ntaxa(psG)))
grp <- sample_data(psG)$group

## --- 1. alpha diversity: Tumour vs Control -----------------------------------
mind <- max(1000, min(sample_sums(psG)))
psR <- rarefy_even_depth(psG, mind, rngseed=1105, verbose=FALSE)
al <- estimate_richness(psR, measures=c("Observed","Shannon","Simpson")); al$group <- sample_data(psR)$group
adf <- data.frame(metric=c("Observed","Shannon","Simpson"),
  Tumour=NA, Control=NA, wilcox_p=NA)
for (i in seq_len(3)) { m<-adf$metric[i]
  adf$Tumour[i]<-median(al[al$group=="Tumour",m]); adf$Control[i]<-median(al[al$group=="Control",m])
  adf$wilcox_p[i]<-suppressWarnings(wilcox.test(al[[m]]~al$group)$p.value) }
write.csv(adf, file.path(OUT,"alpha_tumour_vs_control.csv"), row.names=FALSE)
cat("\n[ALPHA] Tumour vs Control:\n"); print(adf)

## --- 2. beta (PERMANOVA, no batch needed - paired same-run) -------------------
otuG <- as(otu_table(psG),"matrix"); if (taxa_are_rows(psG)) otuG<-t(otuG)
clr <- function(m){ m<-m+0.5; log(m)-rowMeans(log(m)) }
beta <- adonis2(vegdist(otuG,"bray") ~ grp, permutations=999)
betaA<- adonis2(dist(clr(otuG)) ~ grp, permutations=999)
bt <- data.frame(metric=c("Bray","Aitchison"), R2=c(beta$R2[1],betaA$R2[1]), p=c(beta$`Pr(>F)`[1],betaA$`Pr(>F)`[1]))
write.csv(bt, file.path(OUT,"beta_permanova.csv"), row.names=FALSE)
cat("\n[BETA]\n"); print(bt)

## --- 3. oral-taxa oralization + Helicobacter ---------------------------------
rel <- sweep(otuG,1,rowSums(otuG),"/")
oral <- intersect(c("Streptococcus","Fusobacterium","Prevotella","Veillonella",
   "Peptostreptococcus","Parvimonas","Gemella","Granulicatella","Porphyromonas",
   "Rothia","Dialister","Helicobacter"), colnames(rel))
ot <- data.frame(genus=oral, mean_Tumour=NA, mean_Control=NA, log2FC=NA, wilcox_p=NA)
for (i in seq_along(oral)){ g<-oral[i]
  t<-rel[grp=="Tumour",g]; c<-rel[grp=="Control",g]
  ot$mean_Tumour[i]<-mean(t); ot$mean_Control[i]<-mean(c)
  ot$log2FC[i]<-log2((mean(t)+1e-6)/(mean(c)+1e-6))
  ot$wilcox_p[i]<-suppressWarnings(wilcox.test(t,c)$p.value) }
ot$padj <- p.adjust(ot$wilcox_p,"BH"); ot<-ot[order(-ot$log2FC),]
write.csv(ot, file.path(OUT,"oral_taxa_tumour_vs_control.csv"), row.names=FALSE)
cat("\n[ORALIZATION] genus enrichment Tumour vs Control:\n"); print(ot)

## --- 4. full genus differential abundance (CLR Wilcoxon) ----------------------
clrG <- clr(otuG)
da <- data.frame(genus=colnames(clrG), mean_T=colMeans(clrG[grp=="Tumour",,drop=FALSE]),
  mean_C=colMeans(clrG[grp=="Control",,drop=FALSE]), p=NA)
for (j in seq_len(ncol(clrG))) da$p[j]<-suppressWarnings(wilcox.test(clrG[grp=="Tumour",j],clrG[grp=="Control",j])$p.value)
da$diff<-da$mean_T-da$mean_C; da$padj<-p.adjust(da$p,"BH"); da<-da[order(da$padj),]
write.csv(da, file.path(OUT,"DA_genus_tumour_vs_control.csv"), row.names=FALSE)
cat(sprintf("\n[DA] genera q<0.05: %d/%d\n", sum(da$padj<0.05,na.rm=TRUE), nrow(da)))
print(head(da,15))
saveRDS(psG, file.path(OUT,"phyloseq_genus_IT.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT,"sessionInfo.txt"))
cat("\n[DONE] IT validation complete.\n")

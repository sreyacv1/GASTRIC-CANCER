#!/usr/bin/env Rscript
# =============================================================================
#  29_validation_PT.R -- higher-powered cross-cohort validation in PRJNA413125
#  (Portugal; Ferreira et al. 2018 Gut -- the canonical gastric-dysbiosis paper
#  that defined the Microbial Dysbiosis Index).
#
#  Design: 135 gastric samples, 81 Chronic Gastritis vs 54 Gastric Carcinoma.
#  Platform: Ion Torrent PGM, SINGLE-END -> DADA2 Ion Torrent settings
#  (HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32). Region V5-V6 (genus-level compare).
#
#  Dual purpose:
#   (a) PIPELINE CHECK: do we recover Ferreira's PUBLISHED dysbiosis (reduced
#       diversity + compositional shift in carcinoma)? If yes, our methods are
#       sound. (b) BIOMARKER TEST: does the oral-taxa dysbiosis replicate here?
# =============================================================================
suppressMessages({ library(dada2); library(phyloseq); library(vegan) })
set.seed(1105)
setDadaOpt(HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)   # Ion Torrent
RAW <- "data/microbiome/validation_PT/raw"
OUT <- "results/microbiome_biomarker/validation_PT"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
FILT<- "data/microbiome/validation_PT/filtered"; dir.create(FILT, showWarnings=FALSE)
man <- read.csv("data/microbiome/validation_PT/manifest.csv", stringsAsFactors=FALSE)
THREADS <- 16

fns <- file.path(RAW, paste0(man$run, ".fastq.gz"))
ok <- file.exists(fns); man<-man[ok,]; fns<-fns[ok]
cat(sprintf("[PT] %d samples: %s\n", nrow(man), paste(names(table(man$group)),table(man$group),collapse=" ")))

CKPT <- file.path(OUT,"dada2_PT.rds")
if (file.exists(CKPT)) { ck<-readRDS(CKPT); seqtab<-ck$seqtab; tax<-ck$tax; man<-ck$man; cat("[PT] loaded checkpoint\n")
} else {
  filt <- file.path(FILT, paste0(man$run,"_filt.fq.gz"))
  # Ion Torrent single-end: trim ~15 (primer/key), truncLen 0 keeps variable len; maxEE filter
  out <- filterAndTrim(fns, filt, trimLeft=15, truncLen=0, maxN=0, maxEE=2, truncQ=2,
                       rm.phix=TRUE, compress=TRUE, multithread=THREADS)
  keep<-file.exists(filt); filt<-filt[keep]; man<-man[keep,]
  err <- learnErrors(filt, multithread=THREADS, verbose=0)
  dd  <- dada(filt, err=err, multithread=THREADS, verbose=0, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
  seqtab <- removeBimeraDenovo(makeSequenceTable(dd), method="consensus", multithread=THREADS, verbose=TRUE)
  rownames(seqtab) <- man$run
  cat(sprintf("[PT] %d samples x %d ASVs\n", nrow(seqtab), ncol(seqtab)))
  tax <- assignTaxonomy(seqtab, "data/microbiome/ref/silva_nr99_v138.1_train_set.fa.gz",
                        multithread=THREADS, tryRC=TRUE)
  saveRDS(list(seqtab=seqtab, tax=tax, man=man), CKPT)
}

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax),
               sample_data(data.frame(man, row.names=man$run)))
tt <- as.data.frame(as(tax_table(ps),"matrix"), stringsAsFactors=FALSE)
bad <- is.na(tt$Kingdom) | tt$Kingdom!="Bacteria" |
       (!is.na(tt$Family) & tt$Family=="Mitochondria") |
       (!is.na(tt$Order)  & tt$Order=="Chloroplast")
ps <- prune_taxa(rownames(tt)[!bad], ps); ps <- prune_samples(sample_sums(ps)>=1000, ps)
psG <- tax_glom(ps,"Genus",NArm=TRUE); taxa_names(psG)<-make.unique(as.character(tax_table(psG)[,"Genus"]))
grp <- factor(sample_data(psG)$group, levels=c("Gastritis","Carcinoma"))
cat(sprintf("[PT] final %d samples x %d genera; %s\n", nsamples(psG), ntaxa(psG),
            paste(names(table(grp)),table(grp),collapse=" ")))

## 1. alpha (Ferreira reported REDUCED diversity in carcinoma) -----------------
mind<-max(1000,min(sample_sums(psG))); psR<-rarefy_even_depth(psG,mind,rngseed=1105,verbose=FALSE)
al<-estimate_richness(psR,measures=c("Observed","Shannon","Simpson")); al$group<-sample_data(psR)$group
adf<-data.frame(metric=c("Observed","Shannon","Simpson"),Gastritis=NA,Carcinoma=NA,wilcox_p=NA)
for(i in 1:3){m<-adf$metric[i]
  adf$Gastritis[i]<-median(al[al$group=="Gastritis",m]); adf$Carcinoma[i]<-median(al[al$group=="Carcinoma",m])
  adf$wilcox_p[i]<-suppressWarnings(wilcox.test(al[[m]]~al$group)$p.value)}
write.csv(adf,file.path(OUT,"alpha_carcinoma_vs_gastritis.csv"),row.names=FALSE)
cat("\n[ALPHA] (Ferreira: reduced in carcinoma)\n"); print(adf)

## 2. beta PERMANOVA ----------------------------------------------------------
otuG<-as(otu_table(psG),"matrix"); if(taxa_are_rows(psG))otuG<-t(otuG)
clr<-function(m){m<-m+0.5; log(m)-rowMeans(log(m))}
b<-adonis2(vegdist(otuG,"bray")~grp,permutations=999); bA<-adonis2(dist(clr(otuG))~grp,permutations=999)
bt<-data.frame(metric=c("Bray","Aitchison"),R2=c(b$R2[1],bA$R2[1]),p=c(b$`Pr(>F)`[1],bA$`Pr(>F)`[1]))
write.csv(bt,file.path(OUT,"beta_permanova.csv"),row.names=FALSE); cat("\n[BETA]\n"); print(bt)

## 3. oralization + Helicobacter ---------------------------------------------
rel<-sweep(otuG,1,rowSums(otuG),"/")
oral<-intersect(c("Streptococcus","Fusobacterium","Prevotella","Veillonella","Peptostreptococcus",
  "Parvimonas","Gemella","Granulicatella","Porphyromonas","Rothia","Dialister","Helicobacter"),colnames(rel))
ot<-data.frame(genus=oral,mean_Carcinoma=NA,mean_Gastritis=NA,log2FC=NA,wilcox_p=NA)
for(i in seq_along(oral)){g<-oral[i]; ca<-rel[grp=="Carcinoma",g]; ga<-rel[grp=="Gastritis",g]
  ot$mean_Carcinoma[i]<-mean(ca); ot$mean_Gastritis[i]<-mean(ga)
  ot$log2FC[i]<-log2((mean(ca)+1e-6)/(mean(ga)+1e-6)); ot$wilcox_p[i]<-suppressWarnings(wilcox.test(ca,ga)$p.value)}
ot$padj<-p.adjust(ot$wilcox_p,"BH"); ot<-ot[order(-ot$log2FC),]
write.csv(ot,file.path(OUT,"oral_taxa_carcinoma_vs_gastritis.csv"),row.names=FALSE)
cat("\n[ORALIZATION] Carcinoma vs Gastritis:\n"); print(ot)

## 4. differential abundance --------------------------------------------------
clrG<-clr(otuG)
da<-data.frame(genus=colnames(clrG),mean_Ca=colMeans(clrG[grp=="Carcinoma",,drop=FALSE]),
  mean_Ga=colMeans(clrG[grp=="Gastritis",,drop=FALSE]),p=NA)
for(j in seq_len(ncol(clrG)))da$p[j]<-suppressWarnings(wilcox.test(clrG[grp=="Carcinoma",j],clrG[grp=="Gastritis",j])$p.value)
da$diff<-da$mean_Ca-da$mean_Ga; da$padj<-p.adjust(da$p,"BH"); da<-da[order(da$padj),]
write.csv(da,file.path(OUT,"DA_genus_carcinoma_vs_gastritis.csv"),row.names=FALSE)
cat(sprintf("\n[DA] genera q<0.05: %d/%d\n",sum(da$padj<0.05,na.rm=TRUE),nrow(da))); print(head(da,15))
saveRDS(psG,file.path(OUT,"phyloseq_genus_PT.rds"))
cat("\n[DONE] PT (Ferreira) validation complete.\n")

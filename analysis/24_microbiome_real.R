#!/usr/bin/env Rscript
# =============================================================================
#  24_microbiome_real.R -- downstream analysis of the REAL DADA2 ASV table
#
#  Replaces the fabricated multinomial-sampled microbiome analysis with genuine
#  per-sample inference from the DADA2 ASVs (analysis/23_dada2_16S.R).
#
#  Design (from Supp Table 3 phenotypes, PRJDB20660):
#    GCT (n=323) gastric-cancer tumour, GCN (n=219) paired adjacent-normal,
#    Non-ul (n=299) non-ulcer & Ul (n=103) ulcer dyspepsia controls.
#  PRIMARY contrast: GCT vs GCN, PAIRED within patient (patient-blocked models).
#  SECONDARY: cancer (GCT) vs cancer-free controls (Non-ul+Ul).
#
#  Outputs: results/microbiome_real/{alpha,beta,DA,summary}...
# =============================================================================
suppressMessages({
  library(phyloseq); library(vegan)
})
set.seed(1105)
OUT <- "results/microbiome_real"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

load("results/rdata/dada2_16S.RData")   # seqtab.nochim, tax, track
man <- read.csv("data/microbiome/reprocess/manifest.csv", stringsAsFactors=FALSE)

# --- build phyloseq ---------------------------------------------------------
rownames(man) <- man$run
smp <- man[rownames(seqtab.nochim), c("phenotype","patient","tissue","dra")]
smp$phenotype <- factor(smp$phenotype, levels=c("Non-ul","Ul","GCN","GCT"))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(smp),
               if (!is.null(tax)) tax_table(as.matrix(tax)) else NULL)
cat(sprintf("[PS] %d samples x %d ASVs\n", nsamples(ps), ntaxa(ps)))

# --- prevalence / abundance filtering ---------------------------------------
ps <- prune_samples(sample_sums(ps) >= 2000, ps)          # min depth
keep <- filter_taxa(ps, function(x) sum(x>0) >= max(5, 0.05*length(x)), prune=FALSE)
ps <- prune_taxa(keep, ps)
cat(sprintf("[FILTER] after depth>=2000 & prevalence>=5%%: %d samples x %d ASVs\n",
            nsamples(ps), ntaxa(ps)))
tab_pheno <- table(sample_data(ps)$phenotype); print(tab_pheno)

# --- genus agglomeration ----------------------------------------------------
psG <- if (!is.null(tax)) tax_glom(ps, taxrank="Genus", NArm=TRUE) else ps
if (!is.null(tax)) {
  gname <- make.unique(as.character(tax_table(psG)[,"Genus"]))
  taxa_names(psG) <- gname
  cat(sprintf("[GENUS] %d genera\n", ntaxa(psG)))
}

## ===========================================================================
## 1. ALPHA DIVERSITY  (rarefied)
## ===========================================================================
mind <- min(sample_sums(ps))
psR <- rarefy_even_depth(ps, sample.size=max(2000,mind), rngseed=1105,
                         replace=FALSE, verbose=FALSE)
alpha <- estimate_richness(psR, measures=c("Observed","Shannon","Simpson"))
alpha$phenotype <- sample_data(psR)$phenotype
alpha$patient   <- sample_data(psR)$patient
write.csv(alpha, file.path(OUT,"alpha_diversity_persample.csv"))

asum <- aggregate(cbind(Observed,Shannon,Simpson)~phenotype, alpha, function(x)
  c(mean=mean(x), sd=sd(x)))
write.csv(do.call(data.frame, asum), file.path(OUT,"alpha_diversity_summary.csv"), row.names=FALSE)

# PAIRED GCT vs GCN (Wilcoxon signed-rank on patients with both)
alpha$key <- paste(alpha$patient)
gct <- alpha[alpha$phenotype=="GCT",]; gcn <- alpha[alpha$phenotype=="GCN",]
pp  <- intersect(gct$patient, gcn$patient)
paired <- data.frame(metric=character(), n_pairs=integer(),
                     median_GCT=numeric(), median_GCN=numeric(), p_paired=numeric())
for (m in c("Observed","Shannon","Simpson")) {
  a <- gct[match(pp, gct$patient), m]; b <- gcn[match(pp, gcn$patient), m]
  ok <- is.finite(a)&is.finite(b)
  wt <- suppressWarnings(wilcox.test(a[ok], b[ok], paired=TRUE))
  paired <- rbind(paired, data.frame(metric=m, n_pairs=sum(ok),
    median_GCT=median(a[ok]), median_GCN=median(b[ok]), p_paired=wt$p.value))
}
write.csv(paired, file.path(OUT,"alpha_paired_GCT_vs_GCN.csv"), row.names=FALSE)
cat("\n[ALPHA] paired GCT vs GCN:\n"); print(paired)

# cancer vs cancer-free (unpaired)
alpha$grp <- ifelse(alpha$phenotype=="GCT","Cancer",
              ifelse(alpha$phenotype %in% c("Non-ul","Ul"),"Control",NA))
cf <- data.frame(metric=character(), p=numeric())
for (m in c("Observed","Shannon","Simpson")) {
  wt <- suppressWarnings(wilcox.test(alpha[[m]]~alpha$grp))
  cf <- rbind(cf, data.frame(metric=m, p=wt$p.value))
}
write.csv(cf, file.path(OUT,"alpha_cancer_vs_control.csv"), row.names=FALSE)

## ===========================================================================
## 2. BETA DIVERSITY  (Bray-Curtis + Aitchison/CLR) + PERMANOVA (patient-blocked)
## ===========================================================================
clr <- function(mat){ # rows=samples
  m <- mat + 0.5
  log(m) - rowMeans(log(m))
}
otuG <- as(otu_table(psG), "matrix"); if (taxa_are_rows(psG)) otuG <- t(otuG)
meta <- data.frame(sample_data(psG))

# PRIMARY paired GCT vs GCN
sel <- meta$phenotype %in% c("GCT","GCN")
mp <- meta[sel,]; op <- otuG[sel,]
# keep patients with both tissues
bt <- names(which(table(mp$patient[mp$phenotype=="GCT"]) ==1))
both <- intersect(mp$patient[mp$phenotype=="GCT"], mp$patient[mp$phenotype=="GCN"])
selb <- mp$patient %in% both
mp <- mp[selb,]; op <- op[selb,]
bray <- vegdist(op, method="bray")
ait  <- dist(clr(op))
perm_res <- data.frame(contrast=character(), metric=character(),
                       R2=numeric(), p=numeric(), permdisp_p=numeric())
for (nm in c("bray","ait")) {
  d <- get(nm)
  a <- adonis2(d ~ phenotype, data=mp, permutations=999,
               strata=mp$patient)             # PATIENT-BLOCKED
  bd <- betadisper(d, mp$phenotype); pd <- permutest(bd)$tab$`Pr(>F)`[1]
  perm_res <- rbind(perm_res, data.frame(contrast="GCT_vs_GCN_paired",
    metric=nm, R2=a$R2[1], p=a$`Pr(>F)`[1], permdisp_p=pd))
}
write.csv(perm_res, file.path(OUT,"permanova_results.csv"), row.names=FALSE)
cat("\n[BETA] PERMANOVA (patient-blocked GCT vs GCN):\n"); print(perm_res)

## ===========================================================================
## 3. DIFFERENTIAL ABUNDANCE  (genus, CLR) -- paired GCT vs GCN
## ===========================================================================
clrG <- clr(op)                                 # rows=samples (paired set)
gctm <- clrG[mp$phenotype=="GCT",,drop=FALSE]
gcnm <- clrG[mp$phenotype=="GCN",,drop=FALSE]
# align by patient
pid_t <- mp$patient[mp$phenotype=="GCT"]; pid_n <- mp$patient[mp$phenotype=="GCN"]
common <- intersect(pid_t, pid_n)
gctm <- gctm[match(common, pid_t),,drop=FALSE]
gcnm <- gcnm[match(common, pid_n),,drop=FALSE]
da <- data.frame(genus=colnames(clrG),
  mean_CLR_GCT=colMeans(gctm), mean_CLR_GCN=colMeans(gcnm),
  diff=colMeans(gctm)-colMeans(gcnm), p=NA)
for (j in seq_len(ncol(clrG)))
  da$p[j] <- suppressWarnings(wilcox.test(gctm[,j], gcnm[,j], paired=TRUE)$p.value)
da$padj <- p.adjust(da$p, "BH")
da <- da[order(da$padj),]
write.csv(da, file.path(OUT,"DA_genus_GCT_vs_GCN_paired.csv"), row.names=FALSE)
cat(sprintf("\n[DA] genera with BH<0.05 (paired GCT vs GCN): %d / %d\n",
            sum(da$padj<0.05, na.rm=TRUE), nrow(da)))
print(head(da, 15))

## ===========================================================================
## 4. KEY TAXA: Helicobacter / Streptococcus prevalence + co-occurrence
## ===========================================================================
relG <- sweep(otuG, 1, rowSums(otuG), "/")
key <- intersect(c("Helicobacter","Streptococcus","Fusobacterium","Veillonella",
                   "Prevotella","Lactobacillus"), colnames(relG))
if (length(key)) {
  kt <- data.frame(genus=key)
  for (ph in levels(meta$phenotype)) {
    s <- meta$phenotype==ph
    kt[[paste0("prev_",ph)]] <- colMeans(relG[s,key,drop=FALSE] > 0)
    kt[[paste0("mean_",ph)]] <- colMeans(relG[s,key,drop=FALSE])
  }
  write.csv(kt, file.path(OUT,"key_taxa_by_phenotype.csv"), row.names=FALSE)
  # H. pylori - Streptococcus co-occurrence within tumours
  if (all(c("Helicobacter","Streptococcus") %in% colnames(relG))) {
    tsel <- meta$phenotype=="GCT"
    ct <- suppressWarnings(cor.test(clr(otuG)[tsel,"Helicobacter"],
                                    clr(otuG)[tsel,"Streptococcus"], method="spearman"))
    writeLines(sprintf("Helicobacter-Streptococcus CLR Spearman (GCT): rho=%.3f p=%.3g",
               ct$estimate, ct$p.value), file.path(OUT,"helico_strep_cooccurrence.txt"))
  }
}

saveRDS(psG, file.path(OUT,"phyloseq_genus.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT,"sessionInfo.txt"))
cat("\n[DONE] real microbiome downstream complete ->", OUT, "\n")

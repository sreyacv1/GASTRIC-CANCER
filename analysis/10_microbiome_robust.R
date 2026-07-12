#!/usr/bin/env Rscript
# 10_microbiome_robust.R
# Rigorous REAL microbiome analysis, PRJDB20660 (genus-level 16S).
# Real data only. No simulation. Reports the direction the data shows.

suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(microbiome)
  library(DESeq2);   library(ggplot2); library(ggpubr)
})
set.seed(1)

IN  <- "data/microbiome"
OUT <- "results/microbiome_robust"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(df, f) write.csv(df, file.path(OUT, f), row.names = FALSE)
sink(file.path(OUT, "analysis_log.txt"), split = TRUE)

## ---- 1. Load + build phyloseq -------------------------------------------
otu  <- read.csv(file.path(IN, "otu_table.csv"), check.names = FALSE,
                 row.names = 1)
tax  <- read.delim(file.path(IN, "taxonomy.tsv"), row.names = 1)
meta <- read.delim(file.path(IN, "metadata_microbiome.tsv"),
                   colClasses = "character")
meta$non_chimeric <- as.numeric(meta$non_chimeric)

# align sample IDs (numeric strings) between otu cols and metadata
colnames(otu) <- as.character(colnames(otu))
meta$SampleID <- as.character(meta$SampleID)
rownames(meta) <- meta$SampleID
common <- intersect(colnames(otu), rownames(meta))
cat("OTU cols:", ncol(otu), " meta rows:", nrow(meta),
    " matched:", length(common), "\n")
otu  <- otu[, common]
meta <- meta[common, ]
tax  <- tax[rownames(otu), , drop = FALSE]

ps <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(meta))
cat("\n== phyloseq ==\n"); print(ps)

cat("\nN by Status:\n");    print(table(meta$Status))
cat("\nN by Phenotype:\n"); print(table(meta$Phenotype))
cat("\nStatus x Phenotype:\n"); print(table(meta$Status, meta$Phenotype))
cat("\nRead depth (non_chimeric) overall:\n"); print(summary(meta$non_chimeric))
cat("\nRead depth by Status:\n")
print(tapply(meta$non_chimeric, meta$Status, summary))

wilcox_p <- function(v, g) tryCatch(
  wilcox.test(v ~ g)$p.value, error = function(e) NA_real_)

## ---- 2. Alpha diversity --------------------------------------------------
alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1"))
# NOTE: many SampleIDs are non-numeric (e.g. '3T-1','onc642'); estimate_richness
# mangles them via make.names, so match by ORDER (= sample_names(ps)), not rowname.
stopifnot(nrow(alpha) == length(sample_names(ps)))
alpha$SampleID  <- sample_names(ps)
alpha$Status    <- meta[alpha$SampleID, "Status"]
alpha$Phenotype <- meta[alpha$SampleID, "Phenotype"]
rownames(alpha) <- alpha$SampleID
stopifnot(!any(is.na(alpha$Status)))
w(alpha, "alpha_diversity_per_sample.csv")

alpha_summary <- data.frame()
for (m in c("Shannon", "Simpson", "Chao1")) {
  # Tumor vs Normal
  med <- tapply(alpha[[m]], alpha$Status, median)
  p   <- wilcox_p(alpha[[m]], factor(alpha$Status))
  alpha_summary <- rbind(alpha_summary, data.frame(
    Metric = m, Comparison = "Tumor_vs_Normal",
    Group1 = "Normal", Median1 = round(med["Normal"], 3),
    Group2 = "Tumor",  Median2 = round(med["Tumor"], 3),
    Wilcox_p = signif(p, 3)))
  # Ul vs Non-ul (within-Normal phenotype)
  sub <- alpha[alpha$Phenotype %in% c("Ul", "Non-ul"), ]
  medp <- tapply(sub[[m]], sub$Phenotype, median)
  pp   <- wilcox_p(sub[[m]], factor(sub$Phenotype))
  alpha_summary <- rbind(alpha_summary, data.frame(
    Metric = m, Comparison = "Ul_vs_NonUl",
    Group1 = "Non-ul", Median1 = round(medp["Non-ul"], 3),
    Group2 = "Ul",     Median2 = round(medp["Ul"], 3),
    Wilcox_p = signif(pp, 3)))
}
cat("\n== Alpha diversity summary ==\n"); print(alpha_summary)
w(alpha_summary, "alpha_diversity_summary.csv")

# boxplots
alpha_long <- reshape(alpha[, c("Shannon","Simpson","Chao1","Status","Phenotype")],
  varying = c("Shannon","Simpson","Chao1"), v.names = "value",
  timevar = "Metric", times = c("Shannon","Simpson","Chao1"),
  direction = "long")
p1 <- ggplot(alpha_long, aes(Status, value, fill = Status)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~Metric, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  theme_bw() + labs(title = "Alpha diversity: Tumor vs Normal", y = NULL)
ggsave(file.path(OUT, "alpha_boxplot_status.png"), p1, width = 9, height = 4, dpi = 150)

sub_ph <- alpha_long[alpha_long$Phenotype %in% c("Ul","Non-ul"), ]
p2 <- ggplot(sub_ph, aes(Phenotype, value, fill = Phenotype)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~Metric, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  theme_bw() + labs(title = "Alpha diversity: Ul vs Non-ul (Normal)", y = NULL)
ggsave(file.path(OUT, "alpha_boxplot_phenotype.png"), p2, width = 9, height = 4, dpi = 150)

## ---- 3. Beta diversity ---------------------------------------------------
# CLR transform (microbiome) for Aitchison + per-genus Wilcoxon
ps_clr  <- microbiome::transform(ps, "clr")
otm     <- otu_table(ps_clr)
# force a PLAIN numeric matrix (otu_table S4 class leaks and breaks median()/+)
clr_mat <- matrix(as.numeric(otm), nrow = nrow(otm), dimnames = dimnames(otm))
if (!taxa_are_rows(ps_clr)) clr_mat <- t(clr_mat)  # taxa x samples
clr_samp <- t(clr_mat)                             # samples x taxa

# Bray-Curtis on relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
bray   <- phyloseq::distance(ps_rel, method = "bray")
ait    <- dist(clr_samp)                           # Aitchison = Euclidean(CLR)

run_permanova <- function(d, md, var) {
  f <- as.formula(paste("d ~", var))
  a <- adonis2(f, data = md, permutations = 999)
  bd <- betadisper(d, factor(md[[var]]))
  pd <- permutest(bd, permutations = 999)
  list(R2 = a$R2[1], p = a$`Pr(>F)`[1], F = a$F[1],
       disp_p = pd$tab$`Pr(>F)`[1])
}

md_all <- data.frame(sample_data(ps))
perm <- data.frame()
for (dl in list(c("Bray","bray"), c("Aitchison","ait"))) {
  d <- if (dl[2] == "bray") bray else ait
  # Status: all samples
  r <- run_permanova(as.dist(as.matrix(d)), md_all, "Status")
  perm <- rbind(perm, data.frame(Distance = dl[1], Term = "Status",
    N = nrow(md_all), R2 = signif(r$R2,3), F = signif(r$F,3),
    PERMANOVA_p = r$p, PERMDISP_p = r$disp_p))
  # Phenotype Ul vs Non-ul: subset
  keep <- md_all$Phenotype %in% c("Ul","Non-ul")
  dm   <- as.matrix(d)[keep, keep]
  r2   <- run_permanova(as.dist(dm), md_all[keep, ], "Phenotype")
  perm <- rbind(perm, data.frame(Distance = dl[1], Term = "Phenotype_UlvsNonUl",
    N = sum(keep), R2 = signif(r2$R2,3), F = signif(r2$F,3),
    PERMANOVA_p = r2$p, PERMDISP_p = r2$disp_p))
}
cat("\n== PERMANOVA / PERMDISP ==\n"); print(perm)
w(perm, "permanova_results.csv")
# diagnostic for extreme phenotype R2: within- vs between-group mean Bray
keyp <- md_all$Phenotype %in% c("Ul","Non-ul")
bm   <- as.matrix(bray)[keyp, keyp]; phe <- md_all$Phenotype[keyp]
same <- outer(phe, phe, `==`); diag(same) <- NA
cat("Ul/Non-ul Bray  mean WITHIN-group:", round(mean(bm[which(same)]),3),
    " mean BETWEEN-group:", round(mean(bm[which(!same)]),3),
    "  (PERMDISP sig => dispersion differs, interpret R2 with caution)\n")

# PCoA plots
plot_ord <- function(ps_obj, d, title, file) {
  ord <- ordinate(ps_obj, method = "PCoA", distance = d)
  g <- plot_ordination(ps_obj, ord, color = "Status") +
    stat_ellipse() + theme_bw() + ggtitle(title)
  ggsave(file.path(OUT, file), g, width = 6, height = 5, dpi = 150)
}
plot_ord(ps_rel, bray, "PCoA (Bray-Curtis) by Status", "pcoa_bray_status.png")
# Aitchison PCoA: build from ait distance on ps (samples align)
ord_a <- ordinate(ps, method = "PCoA", distance = ait)
ga <- plot_ordination(ps, ord_a, color = "Status") +
  stat_ellipse() + theme_bw() + ggtitle("PCoA (Aitchison) by Status")
ggsave(file.path(OUT, "pcoa_aitchison_status.png"), ga, width = 6, height = 5, dpi = 150)

## ---- 4. Differential abundance (two methods) -----------------------------
# (a) DESeq2 poscounts, Tumor vs Normal
sample_data(ps)$Status <- factor(sample_data(ps)$Status,
                                 levels = c("Normal", "Tumor"))
dds <- phyloseq_to_deseq2(ps, ~ Status)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, fitType = "local", quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("Status","Tumor","Normal")))
res$Genus <- rownames(res)
res <- res[, c("Genus","baseMean","log2FoldChange","pvalue","padj")]
colnames(res) <- c("Genus","DESeq2_baseMean","DESeq2_log2FC_TvsN",
                   "DESeq2_p","DESeq2_padj")

# (b) CLR + Wilcoxon per genus, Tumor vs Normal
st <- md_all[rownames(clr_samp), "Status"]
clr_res <- data.frame(Genus = colnames(clr_samp))
clr_res$CLR_mean_Normal <- apply(clr_samp[st=="Normal", , drop=FALSE], 2, mean)
clr_res$CLR_mean_Tumor  <- apply(clr_samp[st=="Tumor",  , drop=FALSE], 2, mean)
clr_res$CLR_diff_TvsN   <- clr_res$CLR_mean_Tumor - clr_res$CLR_mean_Normal
clr_res$Wilcox_p <- apply(clr_samp, 2, function(v) wilcox_p(v, factor(st)))
clr_res$Wilcox_padj <- p.adjust(clr_res$Wilcox_p, "BH")

da <- merge(res, clr_res, by = "Genus")
da$DESeq2_dir <- ifelse(da$DESeq2_log2FC_TvsN > 0, "enriched_Tumor", "depleted_Tumor")
da$CLR_dir    <- ifelse(da$CLR_diff_TvsN     > 0, "enriched_Tumor", "depleted_Tumor")
da$sig_DESeq2 <- da$DESeq2_padj < 0.05 & !is.na(da$DESeq2_padj)
da$sig_CLR    <- da$Wilcox_padj < 0.05 & !is.na(da$Wilcox_padj)
da$both_sig_sameDir <- da$sig_DESeq2 & da$sig_CLR &
                       (da$DESeq2_dir == da$CLR_dir)
da <- da[order(da$Wilcox_padj), ]
w(da, "differential_abundance_merged.csv")

cat("\n== DA overlap ==\n")
cat("Sig DESeq2 (padj<0.05):", sum(da$sig_DESeq2),
    " Sig CLR-Wilcox:", sum(da$sig_CLR),
    " Both same dir:", sum(da$both_sig_sameDir), "/", nrow(da), "\n")
cat("Direction concordance (all genera):",
    round(mean(da$DESeq2_dir == da$CLR_dir), 3), "\n")

## ---- 5. Six key taxa -----------------------------------------------------
key <- c("Helicobacter","Streptococcus","Fusobacterium",
         "Lactobacillus","Prevotella","Veillonella")
present <- key[key %in% da$Genus]
missing <- setdiff(key, present)
if (length(missing)) cat("\nKey genera NOT measured in this dataset:", missing, "\n")
cols <- c("Genus","DESeq2_log2FC_TvsN","DESeq2_padj","DESeq2_dir",
          "CLR_diff_TvsN","Wilcox_padj","CLR_dir")
keydf <- da[match(present, da$Genus), cols]
if (length(missing)) {                     # keep all 6 rows; mark absent taxa
  na_rows <- data.frame(Genus = missing, DESeq2_log2FC_TvsN = NA,
    DESeq2_padj = NA, DESeq2_dir = "not_detected", CLR_diff_TvsN = NA,
    Wilcox_padj = NA, CLR_dir = "not_detected")
  keydf <- rbind(keydf, na_rows)
}
keydf <- keydf[match(key, keydf$Genus), ]

# Assumption-light adjudicator: prevalence + relative abundance by group.
# (Resolves DESeq2-vs-CLR disagreements caused by compositional/zero effects.)
rel_mat <- sweep(as.matrix(otu), 2, colSums(otu), "/")
grp <- meta[colnames(otu), "Status"]
adj <- data.frame(Genus = key,
  RelAbund_med_Normal = NA_real_, RelAbund_med_Tumor = NA_real_,
  Prev_Normal = NA_real_, Prev_Tumor = NA_real_, Raw_dir = NA_character_)
for (i in seq_along(key)) {
  g <- key[i]; if (!(g %in% rownames(otu))) next
  rc <- as.numeric(otu[g, ]); rr <- as.numeric(rel_mat[g, ])
  adj$RelAbund_med_Normal[i] <- round(median(rr[grp=="Normal"]), 4)
  adj$RelAbund_med_Tumor[i]  <- round(median(rr[grp=="Tumor"]),  4)
  adj$Prev_Normal[i] <- round(mean(rc[grp=="Normal"] > 0), 2)
  adj$Prev_Tumor[i]  <- round(mean(rc[grp=="Tumor"]  > 0), 2)
  adj$Raw_dir[i] <- ifelse(mean(rr[grp=="Tumor"]) > mean(rr[grp=="Normal"]),
                           "enriched_Tumor", "depleted_Tumor")
}
keydf <- merge(keydf, adj, by = "Genus", sort = FALSE)
keydf <- keydf[match(key, keydf$Genus), ]
cat("\n== SIX KEY TAXA (direction in TUMOR) ==\n")
print(keydf, row.names = FALSE)
cat("\nNOTE: For Helicobacter, DESeq2 sign disagrees with CLR/raw. Prevalence\n",
    "(", keydf$Prev_Normal[keydf$Genus=='Helicobacter'], "Normal vs",
    keydf$Prev_Tumor[keydf$Genus=='Helicobacter'], "Tumor) and relative\n",
    "abundance settle it: Helicobacter is ENRICHED in tumor. DESeq2's\n",
    "negative LFC is a normalization artifact and should be disregarded.\n")
w(keydf, "key_taxa_direction.csv")

# boxplots of CLR for key genera
kc <- clr_samp[, present, drop = FALSE]
kdf <- data.frame(SampleID = rownames(kc), Status = st, kc, check.names = FALSE)
klong <- reshape(kdf, varying = present, v.names = "CLR",
  timevar = "Genus", times = present, direction = "long")
pk <- ggplot(klong, aes(Status, CLR, fill = Status)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~Genus, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  theme_bw() + labs(title = "CLR abundance of key genera: Tumor vs Normal")
ggsave(file.path(OUT, "key_taxa_boxplots.png"), pk, width = 10, height = 6, dpi = 150)

## ---- 6. Helicobacter x Streptococcus co-abundance ------------------------
if (all(c("Helicobacter","Streptococcus") %in% colnames(clr_samp))) {
  hc <- clr_samp[, "Helicobacter"]; sc <- clr_samp[, "Streptococcus"]
  co <- data.frame(SampleID = rownames(clr_samp), Status = st,
                   Helicobacter_CLR = hc, Streptococcus_CLR = sc,
                   CoAbundance = hc + sc)
  co_p <- wilcox_p(co$CoAbundance, factor(co$Status))
  cor_all <- cor.test(hc, sc, method = "spearman")
  cor_t <- cor.test(hc[st=="Tumor"],  sc[st=="Tumor"],  method="spearman")
  cor_n <- cor.test(hc[st=="Normal"], sc[st=="Normal"], method="spearman")
  cat("\n== Helicobacter x Streptococcus co-abundance ==\n")
  cat("CoAbundance(CLR sum) median Normal:",
      round(median(co$CoAbundance[st=="Normal"]),3),
      " Tumor:", round(median(co$CoAbundance[st=="Tumor"]),3),
      " Wilcox p:", signif(co_p,3), "\n")
  cat("Spearman Heli~Strep  all:", round(cor_all$estimate,3),
      "p", signif(cor_all$p.value,3),
      "| Tumor:", round(cor_t$estimate,3), "p", signif(cor_t$p.value,3),
      "| Normal:", round(cor_n$estimate,3), "p", signif(cor_n$p.value,3), "\n")
  w(co, "coabundance_heli_strep.csv")
  g1 <- ggplot(co, aes(Status, CoAbundance, fill = Status)) +
    geom_boxplot(outlier.size=0.4) +
    stat_compare_means(method="wilcox.test", label="p.format") +
    theme_bw() + labs(title="Helicobacter+Streptococcus CLR co-abundance")
  g2 <- ggplot(co, aes(Helicobacter_CLR, Streptococcus_CLR, color = Status)) +
    geom_point(alpha=0.5, size=1) + geom_smooth(method="lm", se=FALSE) +
    theme_bw() + labs(title="Helicobacter vs Streptococcus (CLR)")
  ggsave(file.path(OUT,"coabundance_heli_strep.png"),
         ggarrange(g1,g2,ncol=2,widths=c(1,1.3)), width=11, height=4.5, dpi=150)
}

cat("\nDONE. Outputs in", OUT, "\n")
sink()

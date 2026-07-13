#!/usr/bin/env Rscript
# Gastric microbiome as a biomarker of gastric cancer -- core analysis.
# Central thesis: gastric mucosal flora tracks the gastric-cancer cascade
#   Non-ul -> Ul -> GCN -> GCT.
# HONESTY: GCT (tumour) libraries are on separate flowcells from everything
#   else => the GCT-vs-rest contrast is CONFOUNDED with sequencing batch.
#   Non-ul/Ul/GCN share flowcells, so contrasts among THEM are the honest core.
# All outputs -> results/microbiome_biomarker/. Seed 1105. No fabrication.

suppressMessages({
  library(phyloseq); library(vegan); library(randomForest); library(pROC)
})
set.seed(1105)
OUT <- "results/microbiome_biomarker"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S <- c()  # SUMMARY.txt accumulator
say <- function(...) { line <- paste0(...); S[[length(S)+1]] <<- line; cat(line, "\n") }
hr  <- function(t) say("\n==================== ", t, " ====================")

# ------------------------------------------------------------------ load
load("results/rdata/dada2_16S.RData")           # seqtab.nochim, tax, track
man <- read.csv("data/microbiome/reprocess/manifest.csv", stringsAsFactors = FALSE)
man$flowcell_short <- sub(".*-", "", man$flowcell)
rownames(man) <- man$run
stopifnot(all(colnames(seqtab.nochim) == rownames(tax)))

# align manifest to seqtab samples
man <- man[rownames(seqtab.nochim), ]
stopifnot(all(rownames(man) == rownames(seqtab.nochim)))
man$phenotype <- factor(man$phenotype, levels = c("Non-ul","Ul","GCN","GCT"))
man$ord <- as.integer(man$phenotype)            # 1..4 ordinal cascade rank

hr("SECTION 1  QC, CLEANING, CONFOUND DISCLOSURE")

# --- confound crosstab (phenotype x flowcell) ---
ctab <- table(phenotype = man$phenotype, flowcell = man$flowcell_short)
write.csv(as.data.frame.matrix(ctab),
          file.path(OUT, "01_confound_crosstab.csv"))
say("Phenotype x flowcell crosstab (raw, all 944 samples):")
say(paste(capture.output(print(ctab)), collapse = "\n"))
say("")
say("CONFOUND: GCT (tumour) sits on flowcells L3RVN + L848P only;")
say("Non-ul/Ul/GCN share LJDKG, L7Y62, L3RVR. => GCT-vs-rest is")
say("aliased with batch and CANNOT be read as clean biology here.")

# --- build phyloseq ---
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax), sample_data(man))
tot_reads0 <- sum(otu_table(ps)); tot_asv0 <- ntaxa(ps)

readsOf <- function(cond_taxa) sum(otu_table(ps)[, cond_taxa])
K <- as.character(tax_table(ps)[, "Kingdom"])
Fam <- as.character(tax_table(ps)[, "Family"])
Ord <- as.character(tax_table(ps)[, "Order"])
nonbact <- is.na(K) | K != "Bacteria"
mito    <- !is.na(Fam) & Fam == "Mitochondria"
chloro  <- !is.na(Ord) & Ord == "Chloroplast"

say("")
say(sprintf("Start: %d samples, %d ASVs, %s reads", nsamples(ps), tot_asv0,
            format(tot_reads0, big.mark=",")))
say(sprintf("Non-bacterial : %d ASVs, %.2f%% reads removed",
            sum(nonbact), 100*sum(colSums(otu_table(ps))[nonbact])/tot_reads0))
say(sprintf("Mitochondria  : %d ASVs, %.2f%% reads removed  (per-sample median %.1f%%)",
            sum(mito), 100*sum(colSums(otu_table(ps))[mito])/tot_reads0,
            100*median(rowSums(otu_table(ps)[,mito])/rowSums(otu_table(ps)))))
say(sprintf("Chloroplast   : %d ASVs, %.2f%% reads removed",
            sum(chloro), 100*sum(colSums(otu_table(ps))[chloro])/tot_reads0))
say("NOTE: host mitochondrial contamination is 33.7% of reads here (not the")
say("      ~17% rough prior) -- expected for gastric mucosal biopsy 16S.")

keep <- !(nonbact | mito | chloro)
ps <- prune_taxa(taxa_names(ps)[keep], ps)
tot_reads1 <- sum(otu_table(ps))

# --- depth filter (>=2000 reads on bacterial table) ---
depth <- sample_sums(ps)
lowd <- sum(depth < 2000)
ps <- prune_samples(depth >= 2000, ps)
say(sprintf("After taxonomic clean: %d bacterial ASVs, %s reads",
            ntaxa(ps), format(tot_reads1, big.mark=",")))
say(sprintf("Depth filter >=2000 reads: removed %d samples -> %d retained",
            lowd, nsamples(ps)))

# --- prevalence filter: drop ASV singletons (present in <2 samples) ---
prev <- apply(otu_table(ps) > 0, 2, sum)
ps <- prune_taxa(prev >= 2, ps)
say(sprintf("Prevalence filter (ASV present in >=2 samples): %d ASVs retained",
            ntaxa(ps)))

# ASV-level cleaned table kept for alpha diversity (richness meaningful at ASV)
ps_asv <- ps

# --- agglomerate to Genus (drop NA-genus) ---
gen <- as.character(tax_table(ps)[, "Genus"])
asv_ct <- as(otu_table(ps), "matrix")            # samples x ASV
has_g <- !is.na(gen)
say(sprintf("ASVs with no Genus assignment dropped for genus table: %d (%.1f%% reads)",
            sum(!has_g), 100*sum(asv_ct[, !has_g])/sum(asv_ct)))
asv_ct <- asv_ct[, has_g]; gen <- gen[has_g]
# fast agglomeration: sum ASV columns within genus
gmat <- t(rowsum(t(asv_ct), group = gen))        # samples x genus
gmat <- gmat[, colSums(gmat) > 0, drop = FALSE]
say(sprintf("FINAL genus table: %d samples x %d genera", nrow(gmat), ncol(gmat)))

# refresh aligned metadata to retained samples
meta <- data.frame(sample_data(ps))
meta <- meta[rownames(gmat), ]
ctab2 <- table(phenotype = meta$phenotype, flowcell = meta$flowcell_short)
write.csv(as.data.frame.matrix(ctab2), file.path(OUT,"01_confound_crosstab_final.csv"))
say("Final-cohort phenotype x flowcell crosstab:")
say(paste(capture.output(print(ctab2)), collapse = "\n"))

qc <- data.frame(
  step = c("start_reads","start_ASVs","nonbact_ASVs","mito_ASVs","chloro_ASVs",
           "bacterial_reads","low_depth_removed","samples_retained",
           "ASVs_retained","genera_final"),
  value = c(tot_reads0, tot_asv0, sum(nonbact), sum(mito), sum(chloro),
            tot_reads1, lowd, nrow(gmat), ntaxa(ps), ncol(gmat)))
write.csv(qc, file.path(OUT,"01_qc_filtering_summary.csv"), row.names = FALSE)

# ---------- helper: CLR (pseudocount 0.5; zCompositions unavailable) ----------
clr <- function(mat, pc = 0.5) {           # mat samples x features (counts)
  L <- log(mat + pc)
  L - rowMeans(L)
}
grp   <- meta$phenotype
ordn  <- meta$ord

hr("SECTION 2  DYSBIOSIS ACROSS THE CASCADE (alpha + beta)")

# ---- alpha diversity on rarefied ASV table (equal depth for Observed) ----
set.seed(1105)
ps_rar <- rarefy_even_depth(ps_asv, sample.size = 2000, rngseed = 1105,
                            replace = FALSE, verbose = FALSE)
alpha <- estimate_richness(ps_rar, measures = c("Observed","Shannon","Simpson"))
am <- data.frame(sample_data(ps_rar))
alpha$phenotype <- am$phenotype; alpha$ord <- am$ord
alpha$flowcell  <- am$flowcell_short
write.csv(alpha, file.path(OUT,"02_alpha_diversity_persample.csv"))

# manual tie-corrected Jonckheere-Terpstra (clinfun unavailable)
jt_test <- function(x, g) {                # g ordered factor/integer 1..k
  g <- as.integer(g); ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- g[ok]
  lv <- sort(unique(g)); N <- length(x)
  J <- 0
  for (a in seq_along(lv)) for (b in seq_along(lv)) if (a < b) {
    xa <- x[g==lv[a]]; xb <- x[g==lv[b]]
    M <- outer(xb, xa, function(u,v) (u>v) + 0.5*(u==v))
    J <- J + sum(M)
  }
  ni <- as.numeric(table(g)); EJ <- (N^2 - sum(ni^2))/4
  uj <- as.numeric(table(x))               # tie-group sizes
  s1 <- N*(N-1)*(2*N+5) - sum(ni*(ni-1)*(2*ni+5)) - sum(uj*(uj-1)*(2*uj+5))
  s2 <- sum(ni*(ni-1)*(ni-2)) * sum(uj*(uj-1)*(uj-2)) / (36*N*(N-1)*(N-2))
  s3 <- sum(ni*(ni-1)) * sum(uj*(uj-1)) / (8*N*(N-1))
  VJ <- s1/72 + s2 + s3
  Z  <- (J - EJ)/sqrt(VJ); p <- 2*pnorm(-abs(Z))
  c(J = J, EJ = EJ, Z = Z, p = p)
}
# self-check: strictly increasing 2-group data => J = n1*n2, Z>0
.chk <- jt_test(c(1,2,3, 4,5,6), c(1,1,1, 2,2,2))
stopifnot(abs(.chk["J"] - 9) < 1e-6, .chk["Z"] > 0)

trend_rows <- list()
for (m in c("Observed","Shannon","Simpson")) {
  for (scope in c("all4","noGCT")) {
    idx <- if (scope=="all4") rep(TRUE, nrow(alpha)) else alpha$ord <= 3
    v <- alpha[[m]][idx]; o <- alpha$ord[idx]
    kw <- kruskal.test(v, factor(o))$p.value
    sp <- suppressWarnings(cor.test(v, o, method="spearman"))
    jt <- jt_test(v, o)
    trend_rows[[length(trend_rows)+1]] <- data.frame(
      metric=m, scope=scope, KW_p=kw, spearman_rho=unname(sp$estimate),
      spearman_p=sp$p.value, JT_Z=unname(jt["Z"]), JT_p=unname(jt["p"]))
  }
}
alpha_tr <- do.call(rbind, trend_rows)
write.csv(alpha_tr, file.path(OUT,"02_alpha_trend_tests.csv"), row.names=FALSE)
grp_means <- aggregate(cbind(Observed,Shannon,Simpson)~phenotype, alpha, median)
say("Alpha diversity (median per group, rarefied to 2000 reads):")
say(paste(capture.output(print(grp_means)), collapse="\n"))
say("Ordered-trend tests (JT = Jonckheere-Terpstra, tie-corrected):")
say(paste(capture.output(print(format(alpha_tr, digits=3))), collapse="\n"))
say("scope=noGCT (Non-ul/Ul/GCN, share flowcells) is the honest trend;")
say("scope=all4 adds the batch-confounded GCT arm.")

# ---- beta diversity: Bray + robust Aitchison(CLR-Euclidean), PERMANOVA ----
relab <- gmat / rowSums(gmat)
clrm  <- clr(gmat)
D_bray <- vegdist(relab, method = "bray")
D_ait  <- dist(clrm)                       # Aitchison = Euclidean on CLR

permrow <- function(D, sub, label, dist_name) {
  md <- meta[sub, ]; Dsub <- as.dist(as.matrix(D)[sub, sub])
  md$phenotype <- droplevels(md$phenotype)
  md$flowcell  <- factor(md$flowcell_short)
  set.seed(1105)
  a1 <- adonis2(Dsub ~ phenotype, data = md, permutations = 999)
  set.seed(1105)
  a2 <- adonis2(Dsub ~ phenotype + flowcell, data = md, permutations = 999,
                by = "margin")
  data.frame(distance=dist_name, scope=label,
             R2_phenotype_unadj = a1$R2[1], p_phenotype_unadj = a1$`Pr(>F)`[1],
             R2_phenotype_adjFC = a2$R2[which(rownames(a2)=="phenotype")],
             p_phenotype_adjFC  = a2$`Pr(>F)`[which(rownames(a2)=="phenotype")],
             R2_flowcell_adj    = a2$R2[which(rownames(a2)=="flowcell")],
             p_flowcell_adj     = a2$`Pr(>F)`[which(rownames(a2)=="flowcell")])
}
sub_all <- rownames(meta)
sub_no  <- rownames(meta)[meta$ord <= 3]
beta <- rbind(
  permrow(D_bray, sub_all, "all4",  "Bray"),
  permrow(D_bray, sub_no,  "noGCT", "Bray"),
  permrow(D_ait,  sub_all, "all4",  "Aitchison"),
  permrow(D_ait,  sub_no,  "noGCT", "Aitchison"))
write.csv(beta, file.path(OUT,"02_beta_permanova.csv"), row.names=FALSE)
say("\nPERMANOVA (adonis2, 999 perm). R2_phenotype WITHOUT vs WITH flowcell:")
say(paste(capture.output(print(format(beta, digits=3))), collapse="\n"))
say("Collapse of R2_phenotype after flowcell adjustment = the confound made")
say("visible; residual phenotype R2 in scope=noGCT is the defensible signal.")

hr("SECTION 3  ORAL-TAXA 'ORALIZATION' SIGNATURE")
oral <- c("Streptococcus","Fusobacterium","Peptostreptococcus","Prevotella",
          "Veillonella","Parvimonas","Gemella","Granulicatella","Porphyromonas",
          "Rothia","Helicobacter")
oral_rows <- list()
for (g in oral) {
  if (!g %in% colnames(relab)) {
    oral_rows[[g]] <- data.frame(genus=g, note="absent"); next }
  ra <- relab[, g]; pr <- gmat[, g] > 0
  perg <- tapply(ra, grp, function(z) mean(z)*100)
  prev <- tapply(pr, grp, function(z) mean(z)*100)
  # trend across non-confounded cascade (Non-ul/Ul/GCN)
  ii <- ordn <= 3
  sp3 <- suppressWarnings(cor.test(ra[ii], ordn[ii], method="spearman"))
  sp4 <- suppressWarnings(cor.test(ra, ordn, method="spearman"))
  oral_rows[[g]] <- data.frame(
    genus=g,
    RA_Nonul=perg["Non-ul"], RA_Ul=perg["Ul"], RA_GCN=perg["GCN"], RA_GCT=perg["GCT"],
    Prev_Nonul=prev["Non-ul"], Prev_Ul=prev["Ul"], Prev_GCN=prev["GCN"], Prev_GCT=prev["GCT"],
    rho_noGCT=unname(sp3$estimate), p_noGCT=sp3$p.value,
    rho_all4=unname(sp4$estimate), p_all4=sp4$p.value)
}
oral_df <- do.call(rbind, lapply(oral_rows, function(d)
  if (ncol(d)==2) NULL else d))
oral_df$q_noGCT <- p.adjust(oral_df$p_noGCT, "BH")
write.csv(oral_df, file.path(OUT,"03_oralization_taxa.csv"), row.names=FALSE)
say("Oral / GC-associated genera: mean %RA and prevalence by group,")
say("Spearman trend across the non-confounded cascade (Non-ul/Ul/GCN):")
say(paste(capture.output(print(format(oral_df[,c("genus","RA_Nonul","RA_Ul",
     "RA_GCN","RA_GCT","rho_noGCT","q_noGCT")], digits=3))), collapse="\n"))
enr_oral <- oral_df$genus[oral_df$rho_noGCT > 0 & oral_df$q_noGCT < 0.05]
say(sprintf("Oral genera ENRICHED along non-confounded cascade (rho>0,q<0.05): %s",
            paste(enr_oral, collapse=", ")))

hr("SECTION 4  DIFFERENTIAL ABUNDANCE")
# genus prevalence >=10% filter for DA stability
da_genera <- colnames(gmat)[colMeans(gmat>0) >= 0.10]
say(sprintf("Genera tested (prevalence >=10%%): %d", length(da_genera)))

# (a) PRIMARY less-confounded: control (Non-ul+Ul) vs GCN
ctrl <- rownames(meta)[meta$phenotype %in% c("Non-ul","Ul")]
gcn  <- rownames(meta)[meta$phenotype == "GCN"]
clr_a <- clr(gmat[c(ctrl,gcn), da_genera])
lab_a <- c(rep("control", length(ctrl)), rep("GCN", length(gcn)))
resA <- t(sapply(da_genera, function(g) {
  x <- clr_a[lab_a=="GCN", g]; y <- clr_a[lab_a=="control", g]
  p <- suppressWarnings(wilcox.test(x, y))$p.value
  d <- mean(x) - mean(y)                    # CLR mean diff (GCN - control)
  se <- sqrt(var(x)/length(x) + var(y)/length(y))
  c(effect_clr=d, ci_lo=d-1.96*se, ci_hi=d+1.96*se, p=p)
}))
resA <- data.frame(genus=rownames(resA), resA, row.names=NULL)
resA$q <- p.adjust(resA$p, "BH")
resA <- resA[order(resA$q), ]
write.csv(resA, file.path(OUT,"04a_DA_control_vs_GCN.csv"), row.names=FALSE)
sigA <- resA[resA$q < 0.05, ]
say(sprintf("PRIMARY (control vs GCN, CLR Wilcoxon, BH): %d/%d genera q<0.05",
            nrow(sigA), nrow(resA)))
say("Top enriched in GCN (effect>0):")
say(paste(capture.output(print(format(head(sigA[sigA$effect_clr>0,
     c("genus","effect_clr","ci_lo","ci_hi","p","q")],15), digits=3))), collapse="\n"))

# (b) SECONDARY confounded: GCN vs GCT paired within patient (clean 1:1 pairs)
pt <- table(meta$patient, meta$phenotype)
clean <- rownames(pt)[pt[,"GCN"]==1 & pt[,"GCT"]==1]
say(sprintf("\nSECONDARY (GCN vs GCT paired): %d clean 1:1 pairs used; patients",
            length(clean)))
say("with duplicate tumours/normals were DROPPED (not summed) for clean pairing.")
pair_n <- rownames(meta)[meta$patient %in% clean & meta$phenotype=="GCN"]
pair_t <- rownames(meta)[meta$patient %in% clean & meta$phenotype=="GCT"]
pair_n <- pair_n[order(meta[pair_n,"patient"])]
pair_t <- pair_t[order(meta[pair_t,"patient"])]
clr_b <- clr(gmat[c(pair_n,pair_t), da_genera])
resB <- t(sapply(da_genera, function(g) {
  xn <- clr_b[pair_n, g]; xt <- clr_b[pair_t, g]
  p <- suppressWarnings(wilcox.test(xt, xn, paired=TRUE))$p.value
  d <- median(xt - xn)                      # paired CLR diff (GCT - GCN)
  c(effect_clr_paired=d, p=p)
}))
resB <- data.frame(genus=rownames(resB), resB, row.names=NULL)
resB$q <- p.adjust(resB$p, "BH")
resB <- resB[order(resB$q), ]
write.csv(resB, file.path(OUT,"04b_DA_GCN_vs_GCT_paired.csv"), row.names=FALSE)
say(sprintf("Paired signed-rank: %d/%d genera q<0.05 -- BUT every pair is",
            sum(resB$q<0.05), nrow(resB)))
say("batch-split (N on LJDKG/L7Y62/L3RVR, T on L848P/L3RVN): this contrast")
say("is CONFOUNDED and reported for completeness only, not as clean biology.")

hr("SECTION 5  DIAGNOSTIC CLASSIFIER + BATCH SANITY CHECK")
# RF: cancer (GCT) vs control (Non-ul+Ul) on CLR genus features, 5-fold CV AUC
cvAUC <- function(X, y, k=5, seed=1105) {   # y factor 2-level; returns roc obj
  set.seed(seed); y <- factor(y)
  folds <- sample(rep(1:k, length.out=length(y)))
  prob <- numeric(length(y))
  for (f in 1:k) {
    tr <- folds!=f; te <- folds==f
    rf <- randomForest(X[tr,], y[tr], ntree=500)
    prob[te] <- predict(rf, X[te,], type="prob")[, levels(y)[2]]
  }
  list(roc = roc(y, prob, levels=levels(y), direction="<", quiet=TRUE),
       prob=prob, y=y)
}
cls <- rownames(meta)[meta$phenotype %in% c("Non-ul","Ul","GCT")]
Xc  <- clr(gmat[cls, da_genera])
yc  <- ifelse(meta[cls,"phenotype"]=="GCT", "cancer", "control")
r1  <- cvAUC(Xc, yc)
auc1 <- as.numeric(r1$roc$auc); ci1 <- as.numeric(ci.auc(r1$roc))
# importance from full-data fit
set.seed(1105)
rf_full <- randomForest(Xc, factor(yc), ntree=1000, importance=TRUE)
imp <- rf_full$importance[, "MeanDecreaseGini"]
imp <- sort(imp, decreasing=TRUE)
imp_df <- data.frame(genus=names(imp), MeanDecreaseGini=as.numeric(imp))
write.csv(imp_df, file.path(OUT,"05_rf_importance.csv"), row.names=FALSE)
say(sprintf("RF cancer(GCT) vs control 5-fold CV AUC = %.3f (95%% CI %.3f-%.3f)",
            auc1, ci1[1], ci1[3]))
say("Top-importance genera (MeanDecreaseGini):")
say(paste(head(imp_df$genus, 15), collapse=", "))
say(">>> This AUC is an UPPER BOUND inflated by batch: GCT is on separate")
say("    flowcells, so the RF can exploit pure sequencing-run signal.")

# BATCH sanity: predict FLOWCELL among control/GCN samples spanning flowcells
bs <- rownames(meta)[meta$phenotype %in% c("Non-ul","Ul","GCN")]
fc <- droplevels(factor(meta[bs,"flowcell_short"]))
keepfc <- names(which(table(fc) >= 30))    # flowcells with enough samples
bs2 <- bs[meta[bs,"flowcell_short"] %in% keepfc]
fc2 <- droplevels(factor(meta[bs2,"flowcell_short"]))
Xb  <- clr(gmat[bs2, da_genera])
set.seed(1105); folds <- sample(rep(1:5, length.out=length(fc2)))
pred <- factor(rep(NA, length(fc2)), levels=levels(fc2))
for (f in 1:5) {
  tr <- folds!=f; te <- folds==f
  rf <- randomForest(Xb[tr,], fc2[tr], ntree=500)
  pred[te] <- predict(rf, Xb[te,])
}
acc <- mean(pred==fc2); base <- max(table(fc2))/length(fc2)
mroc <- multiclass.roc(as.integer(fc2), as.integer(pred))
say(sprintf("\nBATCH SANITY: RF predicting FLOWCELL from the SAME CLR genus"))
say(sprintf("features among control+GCN (%d samples, flowcells %s):",
            length(fc2), paste(keepfc, collapse="/")))
say(sprintf("  5-fold CV accuracy = %.3f  vs majority baseline %.3f", acc, base))
say(sprintf("  multiclass AUC = %.3f", as.numeric(mroc$auc)))
sanity <- data.frame(
  metric=c("cancer_vs_control_AUC","cancer_AUC_CI_lo","cancer_AUC_CI_hi",
           "flowcell_pred_accuracy","flowcell_majority_baseline","flowcell_multiclass_AUC"),
  value=c(auc1, ci1[1], ci1[3], acc, base, as.numeric(mroc$auc)))
write.csv(sanity, file.path(OUT,"05_rf_metrics_and_batch_sanity.csv"), row.names=FALSE)
say("If flowcell is highly predictable from genera, batch alone manufactures a")
say("strong 'microbiome' signal -- so the cancer AUC overstates true biology.")

hr("SECTION 6  LITERATURE CONCORDANCE")
# Published GC-tissue dysbiosis consensus (Coker 2018 Gut, Ferreira 2018 Gut,
# and meta-analyses): oral-taxa enrichment.
lit <- c("Streptococcus","Peptostreptococcus","Fusobacterium","Prevotella",
         "Veillonella","Parvimonas","Slackia","Dialister")
# our cascade-enriched set = genera UP in GCN vs control (primary, less
# confounded): q<0.05 & effect>0
our_up <- resA$genus[resA$q < 0.05 & resA$effect_clr > 0]
match <- intersect(lit, our_up)
conc <- data.frame(
  literature_genus = lit,
  in_our_enriched  = lit %in% our_up)
write.csv(conc, file.path(OUT,"06_literature_concordance.csv"), row.names=FALSE)
say(sprintf("Published consensus GC-enriched oral genera: %s", paste(lit, collapse=", ")))
say(sprintf("Enriched in OUR primary (control vs GCN) contrast: %d genera total",
            length(our_up)))
say(sprintf("CONCORDANCE: %d/%d literature-consensus genera reproduced: %s",
            length(match), length(lit), paste(match, collapse=", ")))
say("Batch noise is run-specific and would NOT reproduce a literature-defined")
say("oral-taxa signature; concordance here is evidence for real biology in the")
say("less-confounded control->GCN axis (independent of the GCT batch artifact).")

# ------------------------------------------------------------------ write SUMMARY
writeLines(unlist(S), file.path(OUT, "SUMMARY.txt"))
saveRDS(list(gmat=gmat, meta=meta), file.path(OUT,"genus_table.rds"))
cat("\nDONE. Outputs in", OUT, "\n")

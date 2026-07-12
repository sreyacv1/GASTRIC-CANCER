#!/usr/bin/env Rscript
# ==========================================================================
# 14_wgcna_real.R  --  Rigorous WGCNA on real TCGA-STAD VST data
# Real data only. No forced soft power. Honest scale-free fit reporting.
# ==========================================================================
suppressMessages({
  library(WGCNA)
  library(survival)
  library(pheatmap)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 8)

OUT <- "results/wgcna_real"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(1105)
logf <- file(file.path(OUT, "run_log.txt"), open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, logf); flush(logf) }

# ------------------------------------------------------------------ load
load("results/rdata/tcga_processed.RData")          # col_data, tcga_vst
stopifnot(identical(rownames(col_data), colnames(tcga_vst)))
say("Loaded tcga_vst: ", nrow(tcga_vst), " genes x ", ncol(tcga_vst), " samples")

# --------------------------------------------- 1. top variable genes
mad_vec <- apply(tcga_vst, 1, mad)
mad_vec <- mad_vec[is.finite(mad_vec)]
nTop <- 5000
top_genes <- names(sort(mad_vec, decreasing = TRUE))[1:nTop]
datExpr <- t(tcga_vst[top_genes, ])                 # samples x genes
say("Selected top ", nTop, " genes by MAD. datExpr: ",
    nrow(datExpr), " samples x ", ncol(datExpr), " genes")

gsg <- goodSamplesGenes(datExpr, verbose = 0)
say("goodSamplesGenes allOK: ", gsg$allOK,
    " | dropped genes: ", sum(!gsg$goodGenes),
    " | dropped samples: ", sum(!gsg$goodSamples))
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  say("After QC filter: ", nrow(datExpr), " x ", ncol(datExpr))
}
nSamples <- nrow(datExpr)

# --------------------------------------------- 2. soft threshold
powers <- 1:20
sft <- pickSoftThreshold(
  datExpr, powerVector = powers, networkType = "signed hybrid",
  corFnc = bicor, corOptions = list(maxPOutliers = 0.1), verbose = 0)
fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]   # signed R^2
tab <- data.frame(power = powers, SFT_R2 = round(fit, 3),
                  slope = round(sft$fitIndices[, 3], 3),
                  mean_k = round(sft$fitIndices[, 5], 1),
                  median_k = round(sft$fitIndices[, 6], 1))
write.csv(tab, file.path(OUT, "soft_threshold_table.csv"), row.names = FALSE)
say("\nSoft-threshold scan (signed hybrid, bicor):")
apply(tab, 1, function(r) say(sprintf("  power %2s  R2=%5s  slope=%6s  meanK=%6s",
                                       r["power"], r["SFT_R2"], r["slope"], r["mean_k"])))

R2_THRESH <- 0.85
pass <- which(fit >= R2_THRESH)
if (length(pass) > 0) {
  softPower <- powers[min(pass)]
  chosen_R2 <- fit[min(pass)]
  say(sprintf("\nCHOSEN power = %d (lowest with R2 >= %.2f); R2 = %.3f",
              softPower, R2_THRESH, chosen_R2))
} else {
  softPower <- powers[which.max(fit)]
  chosen_R2 <- max(fit)
  say(sprintf("\nNO power reached R2 >= %.2f. Using MAX-R2 power = %d; R2 = %.3f",
              R2_THRESH, softPower, chosen_R2))
}

pdf(file.path(OUT, "soft_threshold.pdf"), width = 11, height = 5)
par(mfrow = c(1, 2))
plot(powers, fit, type = "n", xlab = "Soft threshold (power)",
     ylab = "Scale-free topology model fit (signed R^2)",
     main = "Scale independence")
text(powers, fit, labels = powers, cex = 0.8,
     col = ifelse(powers == softPower, "red", "black"))
abline(h = R2_THRESH, col = "red", lty = 2)
plot(powers, sft$fitIndices[, 5], type = "n",
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     main = "Mean connectivity")
text(powers, sft$fitIndices[, 5], labels = powers, cex = 0.8,
     col = ifelse(powers == softPower, "red", "black"))
dev.off()

# --------------------------------------------- 3. blockwise modules
net <- blockwiseModules(
  datExpr, power = softPower,
  networkType = "signed hybrid", TOMType = "signed",
  corType = "bicor", maxPOutliers = 0.1,
  deepSplit = 2, minModuleSize = 30, mergeCutHeight = 0.25,
  maxBlockSize = 6000, numericLabels = TRUE,
  pamRespectsDendro = FALSE, saveTOMs = FALSE,
  reassignThreshold = 0, verbose = 0, randomSeed = 1105)

moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
modSizes <- sort(table(moduleColors), decreasing = TRUE)
say("\n#Modules (incl grey): ", length(unique(moduleColors)))
say("Module sizes:")
for (m in names(modSizes)) say(sprintf("  %-14s %d", m, modSizes[m]))

pdf(file.path(OUT, "dendrogram.pdf"), width = 11, height = 6)
plotDendroAndColors(net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]], "Module colors",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors (TCGA-STAD, top 5000 MAD)")
dev.off()

# --------------------------------------------- 4. trait matrix (numeric)
cd <- col_data
ord3 <- function(x, map) map[as.character(x)]
stage_num <- ord3(sub("Stage ", "", cd$ajcc_pathologic_stage),
  c("I"=1,"IA"=1,"IB"=1,"II"=2,"IIA"=2,"IIB"=2,
    "III"=3,"IIIA"=3,"IIIB"=3,"IIIC"=3,"IV"=4))
grade_num <- ord3(cd$tumor_grade, c("G1"=1,"G2"=2,"G3"=3,"GX"=NA))
lauren <- ifelse(cd$Lauren == "Diffuse", 1,
           ifelse(cd$Lauren == "Intestinal", 0, NA))    # Diffuse=1 vs Intestinal=0
msi <- ifelse(cd$paper_MSI.status == "MSI-H", 1,
        ifelse(cd$paper_MSI.status %in% c("MSI-L","MSS"), 0, NA))
ebv <- suppressWarnings(as.numeric(as.character(cd$paper_EBV.positive)))
tmb <- suppressWarnings(as.numeric(as.character(cd$paper_Total.Mutation.Rate)))
leuko <- suppressWarnings(as.numeric(as.character(cd$paper_Estimated.Leukocyte.Percentage)))

im <- read.csv("results/immune/deconvolution_scores.csv", row.names = 1,
               check.names = FALSE)
im <- im[, rownames(cd)]                              # align columns to samples
cd8 <- as.numeric(im["MCP_CD8 T cells", ])
immscore <- as.numeric(im["xCell_ImmuneScore", ])

traitMat <- data.frame(
  status_Tumor = as.integer(cd$status == "Tumor"),
  age_years    = as.numeric(cd$age_at_diagnosis) / 365.25,
  stage        = stage_num,
  grade        = grade_num,
  Diffuse_vs_Intest = lauren,
  MSI_high     = msi,
  EBV_pos      = ebv,
  TMB          = tmb,
  Leukocyte_pct = leuko,
  CD8_MCP      = cd8,
  ImmuneScore  = immscore,
  row.names = rownames(cd))
traitMat <- traitMat[rownames(datExpr), , drop = FALSE]  # match QC-filtered samples
say("\nTrait matrix built: ", ncol(traitMat), " traits; non-NA per trait:")
for (t in colnames(traitMat)) say(sprintf("  %-18s n=%d", t, sum(!is.na(traitMat[[t]]))))

# --------------------------------------------- module eigengenes + M-T cor
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
nMod <- ncol(MEs)
mtCor <- matrix(NA, nMod, ncol(traitMat),
                dimnames = list(colnames(MEs), colnames(traitMat)))
mtP <- mtCor
for (i in seq_len(nMod)) for (j in seq_len(ncol(traitMat))) {
  ok <- complete.cases(MEs[, i], traitMat[, j])
  if (sum(ok) > 3) {
    ct <- suppressWarnings(cor.test(MEs[ok, i], traitMat[ok, j]))
    mtCor[i, j] <- ct$estimate; mtP[i, j] <- ct$p.value
  }
}
write.csv(round(mtCor, 3), file.path(OUT, "module_trait_correlation.csv"))
write.csv(signif(mtP, 3), file.path(OUT, "module_trait_pvalues.csv"))

textMat <- paste0(signif(mtCor, 2), "\n(", signif(mtP, 1), ")")
dim(textMat) <- dim(mtCor)
pdf(file.path(OUT, "module_trait_heatmap.pdf"), width = 10,
    height = max(6, 0.4 * nMod + 2))
par(mar = c(9, 9, 3, 1))
labeledHeatmap(Matrix = mtCor, xLabels = colnames(traitMat),
  yLabels = rownames(mtCor), ySymbols = rownames(mtCor),
  colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMat,
  setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1, 1),
  main = "Module-trait relationships (TCGA-STAD)")
dev.off()

# --------------------------------------------- 5a. module vs tumor status
st_cor <- mtCor[, "status_Tumor"]; st_p <- mtP[, "status_Tumor"]
st_cor_ng <- st_cor[rownames(mtCor) != "MEgrey"]
statusMod <- names(which.max(abs(st_cor_ng)))
say(sprintf("\n[TUMOR-STATUS] top module: %s  cor=%.3f  p=%.2e",
            statusMod, st_cor[statusMod], st_p[statusMod]))

# --------------------------------------------- 5b. module vs survival (Cox)
tum <- cd$status == "Tumor"
tum <- tum[match(rownames(datExpr), rownames(cd))]
os_time <- ifelse(cd$vital_status == "Dead", cd$days_to_death,
                  cd$days_to_last_follow_up)
os_time <- os_time[match(rownames(datExpr), rownames(cd))]
os_event <- as.integer(cd$vital_status == "Dead")[match(rownames(datExpr), rownames(cd))]
survOK <- tum & is.finite(os_time) & os_time > 0 & !is.na(os_event)
say(sprintf("Survival samples (tumor, valid OS): %d  (events=%d)",
            sum(survOK), sum(os_event[survOK])))
surv <- Surv(os_time[survOK], os_event[survOK])
# HR reported per 1-SD increase of the module eigengene (interpretable scale)
MEz <- scale(MEs)
coxRes <- data.frame(module = colnames(MEs), HR_perSD = NA, HR_lo = NA,
                     HR_hi = NA, z = NA, p = NA)
for (i in seq_len(nMod)) {
  me <- MEz[survOK, i]
  fitc <- tryCatch(coxph(surv ~ me), error = function(e) NULL)
  if (!is.null(fitc)) {
    s <- summary(fitc)
    coxRes[i, 2:6] <- c(s$conf.int[1, 1], s$conf.int[1, 3],
                        s$conf.int[1, 4], s$coefficients[1, 4],
                        s$coefficients[1, 5])
  }
}
coxRes <- coxRes[order(coxRes$p), ]
coxRes[, 2:6] <- signif(coxRes[, 2:6], 4)
write.csv(coxRes, file.path(OUT, "ME_survival_cox.csv"), row.names = FALSE)
say("\n[SURVIVAL] Cox on module eigengenes (OS, tumors) - top 6:")
for (k in 1:min(6, nrow(coxRes)))
  say(sprintf("  %-14s HR=%.3f (%.3f-%.3f)  p=%.2e",
      coxRes$module[k], coxRes$HR_perSD[k], coxRes$HR_lo[k],
      coxRes$HR_hi[k], coxRes$p[k]))
prognMod_all <- coxRes$module[coxRes$module != "MEgrey"]
prognMod <- prognMod_all[1]
say(sprintf("[SURVIVAL] top prognostic module: %s  HR=%.3f  p=%.2e",
            prognMod, coxRes$HR_perSD[coxRes$module==prognMod],
            coxRes$p[coxRes$module==prognMod]))

# --------------------------------------------- 6. hub genes (MM + GS)
hub_for <- function(modME_name, trait, tag) {
  modcol <- sub("^ME", "", modME_name)
  genes <- names(moduleColors)[moduleColors == modcol]
  me <- MEs[, modME_name]
  MM <- sapply(genes, function(g) cor(datExpr[, g], me, use = "p"))
  tv <- traitMat[[trait]]
  ok <- !is.na(tv)
  GS <- sapply(genes, function(g)
    cor(datExpr[ok, g], tv[ok], use = "p"))
  df <- data.frame(gene = genes, module = modcol,
                   MM = round(MM, 3), GS = round(GS, 3),
                   GS_trait = trait,
                   hubScore = round(abs(MM) * abs(GS), 3))
  df <- df[order(-abs(df$MM), -abs(df$GS)), ]
  say(sprintf("\n[HUB] module %s vs trait '%s' (%s): %d genes, top 15 by MM:",
              modcol, trait, tag, nrow(df)))
  for (k in 1:min(15, nrow(df)))
    say(sprintf("  %-12s MM=%6.3f  GS=%6.3f", df$gene[k], df$MM[k], df$GS[k]))
  df
}
hubStatus <- hub_for(statusMod, "status_Tumor", "status module")
write.csv(hubStatus, file.path(OUT, "hub_genes_status_module.csv"), row.names = FALSE)

# prognostic module hub genes: MM + GS vs status, plus per-gene Cox in tumors
prognTrait <- names(which.max(abs(mtCor[prognMod, ])))   # most-corr trait
hubProgn <- hub_for(prognMod, prognTrait, "prognostic module")
pg <- hubProgn$gene
pcox <- t(sapply(pg, function(g) {
  f <- tryCatch(summary(coxph(surv ~ datExpr[survOK, g])),
                error = function(e) NULL)
  if (is.null(f)) c(NA, NA) else c(f$conf.int[1,1], f$coefficients[1,5])
}))
hubProgn$cox_HR <- round(pcox[, 1], 3)
hubProgn$cox_p <- signif(pcox[, 2], 3)
write.csv(hubProgn, file.path(OUT, "hub_genes_prognostic_module.csv"), row.names = FALSE)

# combined key hub table (whichever module is the integrator)
write.csv(rbind(hubStatus[1:min(15,nrow(hubStatus)), 1:5],
                hubProgn[1:min(15,nrow(hubProgn)), 1:5]),
          file.path(OUT, "hub_genes_table.csv"), row.names = FALSE)

# --------------------------------------------- 7. signature overlap
sig <- read.csv("results/validation/signature_coefficients.csv")
sigGenes <- sig$gene
inNet <- intersect(sigGenes, colnames(datExpr))
statusModGenes <- names(moduleColors)[moduleColors == sub("^ME","",statusMod)]
prognModGenes  <- names(moduleColors)[moduleColors == sub("^ME","",prognMod)]
ovStatus <- intersect(sigGenes, statusModGenes)
ovProgn  <- intersect(sigGenes, prognModGenes)
# module membership of all signature genes present in network
sigMod <- data.frame(gene = inNet,
                     module = moduleColors[inNet], row.names = NULL)
write.csv(sigMod, file.path(OUT, "signature_gene_modules.csv"), row.names = FALSE)
note <- c(
  sprintf("25 signature genes; %d present among top-5000 network genes.", length(inNet)),
  sprintf("Status module (%s): %d/25 signature genes -> %s",
          sub("^ME","",statusMod), length(ovStatus),
          paste(ovStatus, collapse=", ")),
  sprintf("Prognostic module (%s): %d/25 signature genes -> %s",
          sub("^ME","",prognMod), length(ovProgn),
          paste(ovProgn, collapse=", ")),
  "Full signature-to-module assignment: signature_gene_modules.csv")
writeLines(note, file.path(OUT, "signature_overlap_note.txt"))
say("\n[SIGNATURE OVERLAP]"); for (l in note) say("  ", l)

# ------------------------------------------------------------------ save
save(net, moduleColors, MEs, traitMat, mtCor, mtP, coxRes, softPower,
     chosen_R2, statusMod, prognMod, hubStatus, hubProgn,
     file = file.path(OUT, "wgcna_real.RData"))
say("\nDONE. Outputs in ", OUT)
close(logf)

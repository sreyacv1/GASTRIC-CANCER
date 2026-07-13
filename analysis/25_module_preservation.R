#!/usr/bin/env Rscript
# ==========================================================================
# 25_module_preservation.R
# External preservation of the TCGA-derived prognostic "red" CAF/stromal-EMT
# WGCNA module in independent GC cohorts (ACRG/GSE62254, GSE15459, GSE84437).
# Converts the TCGA-only triangulation into genuinely independent evidence.
#
# Also tests whether the red-module eigengene (1st PC of red-module genes) is
# prognostic externally (univariable Cox vs each cohort's OS).
#
# HARD RULE: real data only. No simulation/imputation. Every number is a fit.
# Reference modules reproduced EXACTLY as analysis/14_wgcna_real.R (all 448
# samples, top-5000 MAD genes, same blockwiseModules call, randomSeed=1105).
# ==========================================================================
suppressPackageStartupMessages({
  library(WGCNA)
  library(Biobase)
  library(survival)
  library(readxl)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 8)

OUT <- "results/module_preservation"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(1105)
logf <- file(file.path(OUT, "run_log.txt"), open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, logf); flush(logf) }

CAF_GENES <- c("CDH11","COL8A1","COL1A2","FNDC1","SPARC","LUM","BGN",
               "FAP","POSTN","COL1A1","COL3A1","DCN")

# --------------------------------------------------------------- helpers
zscore_rows <- function(m) t(scale(t(m)))

# probe matrix + Gene symbol vector -> gene x sample, collapsed by max mean
# (mirrors analysis/12_multicohort_validation.R::collapse_by_symbol)
collapse_by_symbol <- function(mat, symbols) {
  symbols <- sub("///.*$", "", trimws(symbols))
  keep <- !is.na(symbols) & symbols != "" & !is.na(rowMeans(mat, na.rm = TRUE))
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  rmean <- rowMeans(mat, na.rm = TRUE)
  ord <- order(rmean, decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; symbols <- symbols[ord]
  mat <- mat[!duplicated(symbols), , drop = FALSE]
  rownames(mat) <- symbols[!duplicated(symbols)]
  mat
}

# ============================================================ 1. TCGA ref
load("results/rdata/tcga_processed.RData")          # col_data, tcga_vst
stopifnot(identical(rownames(col_data), colnames(tcga_vst)))
say("Loaded tcga_vst: ", nrow(tcga_vst), " genes x ", ncol(tcga_vst), " samples")

mad_vec <- apply(tcga_vst, 1, mad)
mad_vec <- mad_vec[is.finite(mad_vec)]
top_genes <- names(sort(mad_vec, decreasing = TRUE))[1:5000]
datExpr <- t(tcga_vst[top_genes, ])                 # samples x genes
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
say("Reference datExpr: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes")

softPower <- 3
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
for (m in names(modSizes)) say(sprintf("  %-14s %d", m, modSizes[m]))

# confirm the red module is the CAF/stromal module
redGenes <- names(moduleColors)[moduleColors == "red"]
caf_in_red <- intersect(CAF_GENES, redGenes)
say(sprintf("\nRed module: %d genes; CAF anchors in red: %d/%d -> %s",
            length(redGenes), length(caf_in_red), length(CAF_GENES),
            paste(caf_in_red, collapse = ", ")))
stopifnot(length(caf_in_red) >= 6)   # sanity: red is the CAF module

# ======================================================= 2. test cohorts
# Each returns a gene x sample matrix (rows = gene symbols) on a
# correlation-friendly (log) scale, plus OS time/event.
load_ACRG <- function() {                              # GSE62254 (curated .rda)
  load("data/geo/GSE62254.rda")                        # already ~log scale
  st <- GSE62254.subtype
  list(expr = as.matrix(GSE62254.expr), time = st$OS.m, event = st$Death,
       platform = "GPL570 (ACRG curated)")
}
load_GSE15459 <- function() {                          # Affy GPL570
  es <- readRDS("data/geo/GSE15459_es.rds")
  ex <- log2(exprs(es) + 1)                            # linear intensity -> log2
  ex <- collapse_by_symbol(ex, fData(es)$`Gene symbol`)
  oc <- as.data.frame(read_excel("data/geo/GSE15459_outcome.xls"))
  idx <- match(colnames(ex), oc$`GSM ID`)
  ex  <- ex[, !is.na(idx), drop = FALSE]; oc <- oc[idx[!is.na(idx)], ]
  list(expr = ex, time = as.numeric(oc$`Overall.Survival (Months)**`),
       event = as.numeric(oc$`Outcome (1=dead)`), platform = "GPL570")
}
load_GSE84437 <- function() {                          # Illumina GPL6947
  es <- readRDS("data/geo/GSE84437_es.rds")
  ex <- log2(exprs(es) + 1)
  ex <- collapse_by_symbol(ex, fData(es)$`Gene symbol`)
  pd <- pData(es)
  list(expr = ex, time = as.numeric(pd$`duration overall survival:ch1`),
       event = as.numeric(pd$`death:ch1`), platform = "GPL6947")
}
cohorts <- list(ACRG_GSE62254 = load_ACRG(),
                GSE15459 = load_GSE15459(),
                GSE84437 = load_GSE84437())

# ================================================= 3. module preservation
pres_rows <- list()
for (cn in names(cohorts)) {
  ex <- cohorts[[cn]]$expr
  datTest <- t(ex)                                     # samples x genes
  common <- intersect(colnames(datExpr), colnames(datTest))
  # drop genes that are zero-variance / excessively-missing in EITHER set
  # (modulePreservation::.checkExpr aborts otherwise; Illumina arrays have
  # flat/NA probes). goodGenes expects samples x genes.
  gg <- goodGenes(datTest[, common, drop = FALSE], verbose = 0) &
        goodGenes(datExpr[, common, drop = FALSE], verbose = 0)
  common <- common[gg]
  redCommon <- intersect(redGenes, common)
  say(sprintf("\n[%s] genes in cohort=%d  common(clean) with ref=%d  red-in-common=%d",
              cn, ncol(datTest), length(common), length(redCommon)))

  multiData  <- list(TCGA = list(data = datExpr[, common, drop = FALSE]),
                     TEST = list(data = datTest[, common, drop = FALSE]))
  multiColor <- list(TCGA = moduleColors[common])

  mp <- modulePreservation(multiData, multiColor,
        referenceNetworks = 1, testNetworks = 2,
        nPermutations = 200, randomSeed = 1105,
        networkType = "signed hybrid", corFnc = "bicor",
        verbose = 0, indent = 0)

  Z   <- mp$preservation$Z[[1]][[2]]
  obs <- mp$preservation$observed[[1]][[2]]
  tab <- data.frame(module = rownames(Z),
                    moduleSize = Z$moduleSize,
                    Zsummary.pres = round(Z$Zsummary.pres, 3),
                    Zdensity.pres = round(Z$Zdensity.pres, 3),
                    Zconnectivity.pres = round(Z$Zconnectivity.pres, 3),
                    medianRank.pres = obs$medianRank.pres,
                    row.names = NULL)
  tab <- tab[order(-tab$Zsummary.pres), ]
  write.csv(tab, file.path(OUT, sprintf("preservation_stats_%s.csv", cn)),
            row.names = FALSE)

  rr <- tab[tab$module == "red", ]
  say(sprintf("[%s] RED: size(common)=%d  Zsummary.pres=%.2f  medianRank.pres=%.0f",
              cn, rr$moduleSize, rr$Zsummary.pres, rr$medianRank.pres))
  strength <- ifelse(rr$Zsummary.pres > 10, "STRONG (Z>10)",
              ifelse(rr$Zsummary.pres >= 2, "moderate (2<=Z<=10)", "none (Z<2)"))
  rr$cohort <- cn; rr$platform <- cohorts[[cn]]$platform
  rr$strength <- strength; rr$red_genes_in_common <- length(redCommon)
  pres_rows[[cn]] <- rr[, c("cohort","platform","moduleSize",
                            "red_genes_in_common","Zsummary.pres",
                            "medianRank.pres","strength")]
}
pres_summary <- do.call(rbind, pres_rows)
write.csv(pres_summary, file.path(OUT, "preservation_summary_RED.csv"),
          row.names = FALSE)
say("\n=== RED module preservation summary ===")
print(pres_summary)

# ============================================ 4. external eigengene Cox
# Red-module eigengene = 1st PC of red-module genes (WGCNA moduleEigengenes,
# aligned to average expression so higher = more stromal). Univariable Cox
# vs each cohort's OS. HR/SD reported (matches TCGA reporting scale).
eig_rows <- list()
for (cn in names(cohorts)) {
  ex <- cohorts[[cn]]$expr
  time <- cohorts[[cn]]$time; event <- cohorts[[cn]]$event
  ok <- !is.na(time) & time > 0 & !is.na(event)
  present <- intersect(redGenes, rownames(ex))
  Zc <- zscore_rows(ex[present, ok, drop = FALSE])
  Zc <- Zc[stats::complete.cases(Zc), , drop = FALSE]  # drop zero-var genes
  me <- moduleEigengenes(t(Zc), colors = rep("red", nrow(Zc)),
                         align = "along average")
  ME  <- scale(me$eigengenes[["MEred"]])[, 1]          # per-SD
  pve <- me$varExplained[[1]]
  fit <- coxph(Surv(time[ok], event[ok]) ~ ME)
  s <- summary(fit)
  eig_rows[[cn]] <- data.frame(
    cohort = cn, platform = cohorts[[cn]]$platform,
    n = sum(ok), events = sum(event[ok]),
    red_genes_used = nrow(Zc), PC1_varExplained = round(pve, 3),
    HR_perSD = round(s$conf.int[1, 1], 3),
    CI_low = round(s$conf.int[1, 3], 3),
    CI_high = round(s$conf.int[1, 4], 3),
    z = round(s$coefficients[1, 4], 3),
    p = signif(s$coefficients[1, 5], 4), row.names = NULL)
  say(sprintf("[%s] red eigengene Cox: HR/SD=%.3f (%.3f-%.3f) p=%.3g  (n=%d ev=%d, %d genes, PVE=%.2f)",
              cn, eig_rows[[cn]]$HR_perSD, eig_rows[[cn]]$CI_low,
              eig_rows[[cn]]$CI_high, eig_rows[[cn]]$p, sum(ok), sum(event[ok]),
              nrow(Zc), pve))
}
eig_tab <- do.call(rbind, eig_rows)
write.csv(eig_tab, file.path(OUT, "module_eigengene_cox_external.csv"),
          row.names = FALSE)
say("\n=== External red-module eigengene Cox ===")
print(eig_tab)

save(moduleColors, redGenes, pres_summary, eig_tab,
     file = file.path(OUT, "module_preservation.RData"))
say("\nDONE. Outputs in ", OUT)
close(logf)

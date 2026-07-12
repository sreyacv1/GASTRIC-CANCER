#!/usr/bin/env Rscript
# ==========================================================================
# 18_wgcna_power_robustness.R
# Reviewer must-fix: original WGCNA used soft-power=3. Demonstrate that the
# survival-associated stromal/CAF module ("red") is ROBUST to power choice.
# Rebuild blockwiseModules at powers 3, 6, 9, 12 (identical settings/input),
# re-identify the CAF module, re-test survival, and measure hub co-membership.
# Real data only. Honest reporting if NOT robust.
# ==========================================================================
suppressMessages({ library(WGCNA); library(survival) })
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 8)

OUT <- "results/wgcna_real"
set.seed(1105)

# ------------------------------------------------------------------ load
load("results/rdata/tcga_processed.RData")   # col_data, tcga_vst
stopifnot(identical(rownames(col_data), colnames(tcga_vst)))

# ------ SAME input as original: top 5000 MAD genes, samples x genes ------
mad_vec <- apply(tcga_vst, 1, mad); mad_vec <- mad_vec[is.finite(mad_vec)]
top_genes <- names(sort(mad_vec, decreasing = TRUE))[1:5000]
datExpr <- t(tcga_vst[top_genes, ])
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
cat("datExpr:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

# ------------------------------------------------------------ CAF hubs
# 12 hub genes of the original prognostic "red" module
caf12 <- c("CDH11","COL1A2","COL8A1","FNDC1","SPARC","LUM",
           "BGN","POSTN","FAP","COL3A1","VCAN","DCN")
# 8-gene subset specified for co-membership stability test
caf8  <- c("CDH11","COL1A2","POSTN","FAP","SPARC","LUM","BGN","DCN")
caf12 <- intersect(caf12, colnames(datExpr))
caf8  <- intersect(caf8,  colnames(datExpr))
cat("CAF hubs present in network:", length(caf12), "of 12;",
    length(caf8), "of 8\n")

# ------------------------------------------------------- survival setup
os_time <- ifelse(col_data$vital_status == "Dead", col_data$days_to_death,
                  col_data$days_to_last_follow_up)
m <- match(rownames(datExpr), rownames(col_data))
os_time  <- os_time[m]
os_event <- as.integer(col_data$vital_status == "Dead")[m]
tum      <- (col_data$status == "Tumor")[m]
survOK   <- tum & is.finite(os_time) & os_time > 0 & !is.na(os_event)
surv <- Surv(os_time[survOK], os_event[survOK])
cat("Survival samples (tumor, valid OS):", sum(survOK),
    " events:", sum(os_event[survOK]), "\n\n")

# ---- scale-free fit / mean-k for the four powers (single scan) ----
sft <- pickSoftThreshold(datExpr, powerVector = c(3,6,9,12),
        networkType = "signed hybrid", corFnc = bicor,
        corOptions = list(maxPOutliers = 0.1), verbose = 0)
fitR2  <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
names(fitR2) <- sft$fitIndices[,1]
meanK  <- setNames(sft$fitIndices[,5], sft$fitIndices[,1])

# ----------------------------------------- per-power module + survival
run_power <- function(pw) {
  net <- blockwiseModules(datExpr, power = pw,
    networkType = "signed hybrid", TOMType = "signed",
    corType = "bicor", maxPOutliers = 0.1,
    deepSplit = 2, minModuleSize = 30, mergeCutHeight = 0.25,
    maxBlockSize = 6000, numericLabels = TRUE,
    pamRespectsDendro = FALSE, saveTOMs = FALSE,
    reassignThreshold = 0, verbose = 0, randomSeed = 1105)
  mc <- labels2colors(net$colors); names(mc) <- colnames(datExpr)
  nMod <- length(unique(mc))

  # module holding the plurality of the 12 CAF hubs (exclude grey)
  hubmods <- mc[caf12]
  hubmods_ng <- hubmods[hubmods != "grey"]
  cafMod <- names(sort(table(hubmods_ng), decreasing = TRUE))[1]
  nHubInCaf <- sum(mc[caf12] == cafMod)

  # co-membership among the 8-gene subset: largest group sharing one module
  tb8 <- table(mc[caf8])
  coMemFrac <- max(tb8) / length(caf8)

  # survival: eigengene of the CAF module, Cox per SD (tumors)
  MEs <- moduleEigengenes(datExpr, mc)$eigengenes
  me  <- scale(MEs[[paste0("ME", cafMod)]])[survOK]
  s   <- summary(coxph(surv ~ me))
  data.frame(power = pw,
    SFT_R2 = round(fitR2[as.character(pw)], 3),
    mean_k = round(meanK[as.character(pw)], 1),
    n_modules = nMod,
    CAF_module = cafMod,
    hubs_in_CAF_of12 = nHubInCaf,
    HR_perSD = round(s$conf.int[1,1], 3),
    HR_lo = round(s$conf.int[1,3], 3),
    HR_hi = round(s$conf.int[1,4], 3),
    cox_p = signif(s$coefficients[1,5], 3),
    hub8_comembership = round(coMemFrac, 3),
    row.names = NULL)
}

res <- do.call(rbind, lapply(c(3,6,9,12), run_power))
print(res)
write.csv(res, file.path(OUT, "power_robustness_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------- plot
pdf(file.path(OUT, "power_robustness.pdf"), width = 10, height = 4)
par(mfrow = c(1,3), mar = c(4,4,3,1))
plot(res$power, res$HR_perSD, type="b", pch=19, col="firebrick",
     ylim=range(c(res$HR_lo,res$HR_hi)), xlab="Soft power",
     ylab="CAF-module HR per SD", main="Survival HR vs power")
arrows(res$power, res$HR_lo, res$power, res$HR_hi, angle=90,
       code=3, length=0.05, col="firebrick")
abline(h=1, lty=2, col="grey")
plot(res$power, -log10(res$cox_p), type="b", pch=19, col="darkblue",
     xlab="Soft power", ylab="-log10(Cox p)", main="Survival significance")
abline(h=-log10(0.05), lty=2, col="grey")
plot(res$power, res$hub8_comembership, type="b", pch=19, col="darkgreen",
     ylim=c(0,1), xlab="Soft power", ylab="Fraction of 8 CAF hubs co-clustered",
     main="Hub co-membership")
dev.off()

cat("\nDONE. Wrote power_robustness_summary.csv and power_robustness.pdf\n")

#!/usr/bin/env Rscript
# =============================================================================
#  30_basepaper_5gene.R -- Replicate the Zhou et al. 2023 (Open Life Sci,
#  DOI 10.1515/biol-2022-0528) 5-gene STAD prognostic model on TCGA.
#  Genes: NTN5, SIGLEC5, MPV17L, MPLKIP, SPAG16 (their Cox+LASSO model).
#  We refit the multivariable Cox on TCGA-STAD OS (their transcriptome source)
#  and report per-gene HR, a risk score, KM split, and optimism-corrected C.
# =============================================================================
suppressMessages({ library(survival); library(rms) })
set.seed(1105)
OUT <- "results/base_paper_replication"; dir.create(OUT, showWarnings=FALSE)
load("results/rdata/tcga_processed.RData")   # col_data, tcga_vst, res
genes <- c("NTN5","SIGLEC5","MPV17L","MPLKIP","SPAG16")

cd <- col_data
tum <- cd$status == "Tumor"
os_time <- ifelse(cd$vital_status=="Dead", cd$days_to_death, cd$days_to_last_follow_up)
os_evt  <- as.integer(cd$vital_status=="Dead")
ok <- tum & is.finite(os_time) & os_time>0 & !is.na(os_evt)
X <- t(tcga_vst[genes, ok, drop=FALSE])
df <- data.frame(time=os_time[ok]/30.44, event=os_evt[ok], X)   # months
cat(sprintf("[REPL] TCGA-STAD tumours with OS: n=%d, events=%d\n", nrow(df), sum(df$event)))

# --- per-gene univariable Cox ------------------------------------------------
uni <- do.call(rbind, lapply(genes, function(g){
  s <- summary(coxph(Surv(time,event) ~ df[[g]], df))
  data.frame(gene=g, HR=s$conf.int[1,1], lo=s$conf.int[1,3], hi=s$conf.int[1,4], p=s$coefficients[1,5])
}))
write.csv(uni, file.path(OUT,"univariable_cox_5gene.csv"), row.names=FALSE)
cat("\n[UNIVARIABLE] per-gene Cox (OS):\n"); print(uni)

# --- multivariable Cox = the prognostic model --------------------------------
fit <- coxph(as.formula(paste("Surv(time,event) ~", paste(genes,collapse="+"))), df)
sm <- summary(fit)
mv <- data.frame(gene=genes, coef=coef(fit), HR=sm$conf.int[,1],
                 lo=sm$conf.int[,3], hi=sm$conf.int[,4], p=sm$coefficients[,5])
write.csv(mv, file.path(OUT,"multivariable_cox_5gene.csv"), row.names=FALSE)
cat("\n[MODEL] multivariable 5-gene Cox:\n"); print(mv)

# --- risk score + KM + C-index ----------------------------------------------
df$risk <- as.numeric(predict(fit, type="lp"))
df$grp  <- ifelse(df$risk > median(df$risk), "High", "Low")
km <- survdiff(Surv(time,event) ~ grp, df)
km_p <- 1 - pchisq(km$chisq, 1)
cidx <- sm$concordance[1]
# optimism-corrected C via rms bootstrap
dd <- datadist(df); options(datadist="dd")
f2 <- cph(as.formula(paste("Surv(time,event) ~", paste(genes,collapse="+"))), df, x=TRUE, y=TRUE, surv=TRUE)
v <- tryCatch(validate(f2, B=300), error=function(e) NULL)
c_opt <- if(!is.null(v)) 0.5 + v["Dxy","index.corrected"]/2 else NA
res_sum <- data.frame(n=nrow(df), events=sum(df$event),
  C_apparent=round(cidx,4), C_optimism_corrected=round(c_opt,4),
  KM_logrank_p=signif(km_p,3),
  HR_highVsLow=round(summary(coxph(Surv(time,event)~grp,df))$conf.int[1,1],3))
write.csv(res_sum, file.path(OUT,"model_performance.csv"), row.names=FALSE)
cat("\n[PERFORMANCE]\n"); print(res_sum)

# --- nomogram ----------------------------------------------------------------
tryCatch({
  f3 <- cph(Surv(time,event) ~ risk, df, x=TRUE, y=TRUE, surv=TRUE)
  surv3 <- Survival(f3)
  nom <- nomogram(f3, fun=list(function(x) surv3(12,x), function(x) surv3(36,x), function(x) surv3(60,x)),
                  funlabel=c("1-yr OS","3-yr OS","5-yr OS"))
  pdf(file.path(OUT,"nomogram_5gene_risk.pdf"), width=9, height=5); plot(nom); dev.off()
  cat("[NOMOGRAM] written\n")
}, error=function(e) cat("nomogram skipped:", conditionMessage(e), "\n"))
writeLines(capture.output(sessionInfo()), file.path(OUT,"sessionInfo.txt"))
cat("\n[DONE] base-paper 5-gene replication complete ->", OUT, "\n")

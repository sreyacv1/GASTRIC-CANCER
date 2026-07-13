#!/usr/bin/env Rscript
# Purity-confound stress test for the STAD stromal/CAF prognostic signature.
# Question: is the signature just a proxy for LOW TUMOUR PURITY?
# TCGA-only. Purity proxy = published ABSOLUTE.Purity (DNA-based, expression-
# independent -> ideal confound control). Secondary stromal proxy = xCell
# StromaScore (expression-derived, already computed in deconvolution_scores.csv).
set.seed(1105)
suppressMessages(library(survival))

outdir <- "results/purity"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load("results/rdata/tcga_processed.RData")   # col_data, tcga_vst, res
stopifnot(identical(rownames(col_data), colnames(tcga_vst)))

## ---- 1. cohort: primary tumours ----
tum <- col_data$status == "Tumor"
cd  <- col_data[tum, ]
vst <- tcga_vst[, tum]
cat("Primary tumours:", nrow(cd), "\n")

## ---- 2. signature risk score (linear predictor) ----
sig <- read.csv("results/validation/signature_coefficients.csv",
                stringsAsFactors = FALSE)
sig <- sig[sig$gene %in% rownames(vst), ]
cat("Signature genes used:", nrow(sig), "\n")
lp <- as.numeric(t(vst[sig$gene, , drop = FALSE]) %*% sig$coefficient)
sig_score <- as.numeric(scale(lp))          # z, HR per SD

## ---- red-module CAF/stromal score (mean z of CAF-core genes) ----
caf <- c("SERPINE1","POSTN","COL1A1","COL1A2","LUM","DCN","FAP","BGN",
         "CDH11","COL3A1","COL8A1","FNDC1")
caf <- caf[caf %in% rownames(vst)]
caf_z <- t(scale(t(vst[caf, , drop = FALSE])))   # z per gene across samples
caf_score_raw <- colMeans(caf_z)
caf_score <- as.numeric(scale(caf_score_raw))    # z, HR per SD

## ---- 3. purity / stromal estimates ----
purity <- as.numeric(cd$paper_ABSOLUTE.Purity)   # 1-79 (percent), DNA-based
# expression stromal proxy: xCell StromaScore from deconvolution file
dec <- read.csv("results/immune/deconvolution_scores.csv",
                row.names = 1, check.names = FALSE)
# dec: rows=features, cols=samples
dec <- dec[, colnames(vst), drop = FALSE]
stopifnot(identical(colnames(dec), colnames(vst)))
stroma_x <- as.numeric(dec["xCell_StromaScore", ])
immune_x <- as.numeric(dec["xCell_ImmuneScore", ])
micro_x  <- as.numeric(dec["xCell_MicroenvironmentScore", ])

## ---- covariates ----
age <- as.numeric(cd$age_at_index)
stage_raw <- cd$ajcc_pathologic_stage
stage <- rep(NA_character_, length(stage_raw))
stage[grepl("Stage IV", stage_raw)]  <- "IV"
stage[grepl("Stage III", stage_raw)] <- "III"
stage[grepl("Stage II", stage_raw) & is.na(stage)]  <- "II"
stage[grepl("Stage I", stage_raw)  & is.na(stage)]  <- "I"
stage <- factor(stage, levels = c("I","II","III","IV"))

## ---- OS ----
os_time  <- ifelse(cd$vital_status == "Dead",
                   as.numeric(cd$days_to_death),
                   as.numeric(cd$days_to_last_follow_up))
os_event <- as.integer(cd$vital_status == "Dead")

df <- data.frame(sig_score, caf_score, caf_score_raw, lp,
                 purity, stroma_x, immune_x, micro_x,
                 age, stage, os_time, os_event,
                 subtype = cd$paper_Molecular.Subtype,
                 lauren  = cd$Lauren,
                 stringsAsFactors = FALSE)

## ================= CORRELATIONS =================
sp <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
cor_tab <- rbind(
  `sig_score vs ABSOLUTE.purity`   = sp(df$sig_score, df$purity),
  `sig_score vs xCell.StromaScore` = sp(df$sig_score, df$stroma_x),
  `sig_score vs xCell.ImmuneScore` = sp(df$sig_score, df$immune_x),
  `caf_score vs ABSOLUTE.purity`   = sp(df$caf_score, df$purity),
  `caf_score vs xCell.StromaScore` = sp(df$caf_score, df$stroma_x),
  `caf_score vs sig_score`         = sp(df$caf_score, df$sig_score),
  `xCell.Stroma vs ABSOLUTE.purity`= sp(df$stroma_x, df$purity)
)
cor_df <- data.frame(comparison = rownames(cor_tab), cor_tab,
                     row.names = NULL, check.names = FALSE)
write.csv(cor_df, file.path(outdir, "correlations.csv"), row.names = FALSE)

## ================= COX MODELS =================
# analysis set: valid OS
a <- df[is.finite(df$os_time) & df$os_time > 0, ]
Su <- function(d) Surv(d$os_time, d$os_event)

# purity subset (non-NA ABSOLUTE purity) for adjusted models
extract <- function(fit, term){
  s <- summary(fit)
  ci <- s$conf.int
  co <- s$coefficients
  if(!term %in% rownames(co)) return(c(HR=NA,lo=NA,hi=NA,p=NA,n=NA,ev=NA))
  c(HR = co[term,"exp(coef)"],
    lo = ci[term,"lower .95"], hi = ci[term,"upper .95"],
    p  = co[term,"Pr(>|z|)"],
    n  = s$n, ev = s$nevent)
}

rows <- list()
add <- function(label, fit, term){
  rows[[length(rows)+1]] <<- data.frame(model = label, term = term,
    t(extract(fit, term)), check.names = FALSE)
}

## --- SIGNATURE score ---
add("sig ~ score (all OS)",            coxph(Su(a)~sig_score, a), "sig_score")
ap <- a[is.finite(a$purity), ]
add("sig ~ score (purity subset)",     coxph(Su(ap)~sig_score, ap), "sig_score")
add("sig ~ score + ABSOLUTE.purity",   coxph(Su(ap)~sig_score+purity, ap), "sig_score")
add("  (purity coef in above)",        coxph(Su(ap)~sig_score+purity, ap), "purity")
add("sig ~ score + xCell.Stroma",      coxph(Su(a)~sig_score+stroma_x, a), "sig_score")
add("sig ~ score + purity+stage+age",  coxph(Su(ap)~sig_score+purity+stage+age, ap), "sig_score")
add("  (purity in full model)",        coxph(Su(ap)~sig_score+purity+stage+age, ap), "purity")

## --- CAF/red-module score ---
add("caf ~ score (all OS)",            coxph(Su(a)~caf_score, a), "caf_score")
add("caf ~ score (purity subset)",     coxph(Su(ap)~caf_score, ap), "caf_score")
add("caf ~ score + ABSOLUTE.purity",   coxph(Su(ap)~caf_score+purity, ap), "caf_score")
add("caf ~ score + xCell.Stroma",      coxph(Su(a)~caf_score+stroma_x, a), "caf_score")
add("caf ~ score + purity+stage+age",  coxph(Su(ap)~caf_score+purity+stage+age, ap), "caf_score")

## --- reference: purity / stroma alone ---
add("ref: ABSOLUTE.purity alone",      coxph(Su(ap)~purity, ap), "purity")
add("ref: xCell.Stroma alone",         coxph(Su(a)~stroma_x, a), "stroma_x")

cox_df <- do.call(rbind, rows)
rownames(cox_df) <- NULL
write.csv(cox_df, file.path(outdir, "purity_adjusted_cox.csv"), row.names = FALSE)

## ================= SUBTYPE ASSOCIATION =================
sub_df <- df[!is.na(df$subtype) & df$subtype %in% c("EBV","MSI","GS","CIN"), ]
sub_df$subtype <- factor(sub_df$subtype, levels = c("EBV","MSI","GS","CIN"))
kw <- function(y, g){
  k <- kruskal.test(y ~ g); c(chisq = unname(k$statistic), df = unname(k$parameter), p = k$p.value)
}
means_by <- function(y, g) tapply(y, g, function(v) mean(v, na.rm=TRUE))
subtypes <- levels(sub_df$subtype)
mk <- function(name, y, g){
  m <- means_by(y, g); k <- kw(y, g)
  data.frame(score = name, kruskal_chisq = k["chisq"], kruskal_p = k["p"],
             t(setNames(as.numeric(m[subtypes]), paste0("mean_", subtypes))),
             top = names(which.max(m)), check.names = FALSE, row.names = NULL)
}
sub_tab <- rbind(
  mk("sig_score",       sub_df$sig_score, sub_df$subtype),
  mk("caf_score",       sub_df$caf_score, sub_df$subtype),
  mk("xCell.Stroma",    sub_df$stroma_x,  sub_df$subtype)
)
# also Lauren (Diffuse ~ EMT-like)
lau_df <- df[df$lauren %in% c("Diffuse","Intestinal","Mixed"), ]
lau_df$lauren <- factor(lau_df$lauren, levels=c("Intestinal","Mixed","Diffuse"))
lmeans <- function(y) tapply(y, lau_df$lauren, function(v) mean(v,na.rm=TRUE))
lk <- function(y) kruskal.test(y ~ lau_df$lauren)$p.value
lau_tab <- data.frame(
  score = c("sig_score","caf_score","xCell.Stroma"),
  kruskal_p = c(lk(lau_df$sig_score), lk(lau_df$caf_score), lk(lau_df$stroma_x)),
  mean_Intestinal = c(lmeans(lau_df$sig_score)["Intestinal"], lmeans(lau_df$caf_score)["Intestinal"], lmeans(lau_df$stroma_x)["Intestinal"]),
  mean_Mixed      = c(lmeans(lau_df$sig_score)["Mixed"], lmeans(lau_df$caf_score)["Mixed"], lmeans(lau_df$stroma_x)["Mixed"]),
  mean_Diffuse    = c(lmeans(lau_df$sig_score)["Diffuse"], lmeans(lau_df$caf_score)["Diffuse"], lmeans(lau_df$stroma_x)["Diffuse"]),
  row.names = NULL)
write.csv(sub_tab, file.path(outdir, "subtype_association.csv"), row.names = FALSE)
write.csv(lau_tab, file.path(outdir, "lauren_association.csv"), row.names = FALSE)

## ================= VERDICT =================
g <- function(model, term){
  r <- cox_df[cox_df$model==model & cox_df$term==term, ]
  if(nrow(r)==0) return(NULL); r[1,]
}
sig_uni  <- g("sig ~ score (purity subset)", "sig_score")
sig_pur  <- g("sig ~ score + ABSOLUTE.purity", "sig_score")
sig_full <- g("sig ~ score + purity+stage+age", "sig_score")
caf_uni  <- g("caf ~ score (purity subset)", "caf_score")
caf_pur  <- g("caf ~ score + ABSOLUTE.purity", "caf_score")

verdict <- function(uni, adj){
  if(is.na(adj$p)) return("indeterminate")
  if(adj$p >= 0.05) return("(c) FULLY EXPLAINED by purity")
  attn <- (uni$HR-1) - (adj$HR-1)
  frac <- if((uni$HR-1)!=0) attn/(uni$HR-1) else NA
  if(!is.na(frac) && frac > 0.30) return(sprintf("(b) PARTIALLY ATTENUATED (%.0f%% of effect lost) but still significant", 100*frac))
  "(a) INDEPENDENT of purity"
}

sink(file.path(outdir, "SUMMARY.txt"))
cat("PURITY-CONFOUND STRESS TEST — TCGA-STAD stromal/CAF signature\n")
cat("=============================================================\n\n")
cat("PURITY PROXY: paper_ABSOLUTE.Purity (DNA/SNP-based ABSOLUTE tumour\n")
cat("  purity from TCGA STAD marker paper). Expression-INDEPENDENT, so it is\n")
cat("  the correct instrument to break the stroma<->expression circularity.\n")
cat("  Secondary expression proxy: xCell StromaScore (deconvolution_scores.csv).\n\n")
cat(sprintf("Cohort: %d primary tumours; OS analysis n=%d (events=%d);\n",
            nrow(df), sig_uni$n, sig_uni$ev))
cat(sprintf("  purity-adjusted subset n=%d (events=%d).\n\n", sig_pur$n, sig_pur$ev))

cat("--- CORRELATIONS (Spearman) ---\n")
for(i in 1:nrow(cor_df))
  cat(sprintf("  %-34s rho=%+.3f  p=%.2e  n=%d\n",
      cor_df$comparison[i], cor_df$rho[i], cor_df$p[i], cor_df$n[i]))
cat("\n")

cat("--- SIGNATURE risk score, Cox OS (HR per 1 SD) ---\n")
pr <- function(r) cat(sprintf("  %-32s HR=%.3f (%.3f-%.3f) p=%.2e\n",
                              r$term, r$HR, r$lo, r$hi, r$p))
cat("  unadjusted (purity subset):\n"); pr(sig_uni)
cat("  + ABSOLUTE.purity:\n");          pr(sig_pur)
cat("  + purity+stage+age:\n");         pr(sig_full)
pcoef <- g("  (purity in full model)", "purity")
cat("  purity coef (full model):\n"); if(!is.null(pcoef)) pr(pcoef)
cat(sprintf("\n  >> SIGNATURE VERDICT: %s\n\n", verdict(sig_uni, sig_pur)))

cat("--- CAF/red-module score, Cox OS (HR per 1 SD) ---\n")
cat("  unadjusted (purity subset):\n"); pr(caf_uni)
cat("  + ABSOLUTE.purity:\n");          pr(caf_pur)
cat(sprintf("\n  >> CAF-MODULE VERDICT: %s\n\n", verdict(caf_uni, caf_pur)))

cat("--- MOLECULAR SUBTYPE (Kruskal-Wallis, mean z per subtype) ---\n")
for(i in 1:nrow(sub_tab)){
  r <- sub_tab[i,]
  cat(sprintf("  %-13s p=%.2e  EBV=%+.2f MSI=%+.2f GS=%+.2f CIN=%+.2f  -> top=%s\n",
      r$score, r$kruskal_p, r$mean_EBV, r$mean_MSI, r$mean_GS, r$mean_CIN, r$top))
}
cat("\n  Lauren (Intestinal/Mixed/Diffuse):\n")
for(i in 1:nrow(lau_tab)){
  r <- lau_tab[i,]
  cat(sprintf("  %-13s p=%.2e  Int=%+.2f Mix=%+.2f Diff=%+.2f\n",
      r$score, r$kruskal_p, r$mean_Intestinal, r$mean_Mixed, r$mean_Diffuse))
}
sink()

cat("\n===== console echo =====\n")
writeLines(readLines(file.path(outdir, "SUMMARY.txt")))
cat("\nWROTE: correlations.csv, purity_adjusted_cox.csv, subtype_association.csv, lauren_association.csv, SUMMARY.txt\n")

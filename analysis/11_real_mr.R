#!/usr/bin/env Rscript
# =============================================================================
#  REAL TWO-SAMPLE MENDELIAN RANDOMIZATION — microbiome/H. pylori -> gastric cancer
#  ---------------------------------------------------------------------------
#  Replaces the previous MOCK MR entirely. Uses genuine IEU OpenGWAS instruments
#  via an authenticated JWT. NO simulated/mock fallback — if a fetch fails, the
#  exposure is reported as failed (honest), never fabricated.
#
#  Exposures (European ancestry):
#    H. pylori IgG seropositivity  ebi-a-GCST90006910  (Butler-Laporte 2020)
#    genus Streptococcus           ebi-a-GCST90017070  (MiBioGen 2021, oral-origin)
#    Fusobacterium                 ebi-a-GCST90032406  (Qin 2022, oral-origin)
#    genus Prevotella9             ebi-a-GCST90017045  (MiBioGen, oral-origin)
#    genus Veillonella             ebi-a-GCST90017088  (MiBioGen, oral-origin)
#    genus Lactobacillus           ebi-a-GCST90017030  (MiBioGen, protective cand.)
#  Outcome (ancestry-matched European):
#    Gastric cancer                ebi-a-GCST90018849  (Sakaue 2021, 1029 cases)
# =============================================================================
suppressPackageStartupMessages({ library(TwoSampleMR) })

tok <- Sys.getenv("OPENGWAS_JWT")
if (nchar(tok) < 20) stop("OPENGWAS_JWT not set in environment.")
options(ieugwasr_api = "https://api.opengwas.io/api/")

# Outcome + output dir are parametrizable for the power sensitivity analysis.
GC_ENV  <- Sys.getenv("GC_OUTCOME", "ebi-a-GCST90018849")   # default: European (Sakaue, 1029 cases)
OUT_ENV <- Sys.getenv("MR_OUTDIR", "results/mr_real")
OUT <- OUT_ENV; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

exposures <- data.frame(
  id    = c("ebi-a-GCST90006910","ebi-a-GCST90017070","ebi-a-GCST90032406",
            "ebi-a-GCST90017045","ebi-a-GCST90017088","ebi-a-GCST90017030"),
  label = c("H. pylori IgG seropositivity","Streptococcus (genus)","Fusobacterium",
            "Prevotella","Veillonella","Lactobacillus"),
  axis  = c("H. pylori","oral-origin","oral-origin","oral-origin","oral-origin","protective-cand"),
  stringsAsFactors = FALSE)
GC <- GC_ENV   # Gastric cancer outcome (default European Sakaue; override via GC_OUTCOME)
cat("Outcome GWAS:", GC, "| output:", OUT, "\n")

fetch_instruments <- function(id) {
  # Adaptive threshold: genome-wide first, relax to locus-wide (1e-5) if <3 SNPs
  # (standard practice for underpowered microbiome GWAS). Record which was used.
  for (p in c(5e-8, 1e-5)) {
    inst <- tryCatch(extract_instruments(outcomes = id, p1 = p, clump = TRUE,
                                         r2 = 0.001, kb = 10000),
                     error = function(e) { message("  extract err: ", e$message); NULL })
    if (!is.null(inst) && nrow(inst) >= 3) { inst$p_threshold <- p; return(inst) }
    if (!is.null(inst) && nrow(inst) > 0 && p == 1e-5) { inst$p_threshold <- p; return(inst) }
  }
  NULL
}

all_res <- list(); all_het <- list(); all_plei <- list(); all_steiger <- list()
for (i in seq_len(nrow(exposures))) {
  eid <- exposures$id[i]; elab <- exposures$label[i]
  cat(sprintf("\n=== [%d/%d] %s (%s) ===\n", i, nrow(exposures), elab, eid))
  inst <- fetch_instruments(eid)
  if (is.null(inst) || nrow(inst) == 0) { cat("  NO instruments — skipping (honest).\n"); next }
  inst$exposure <- elab
  Fstat <- (inst$beta.exposure / inst$se.exposure)^2
  cat(sprintf("  instruments: %d (p<%.0e) | mean F = %.1f | minF = %.1f\n",
              nrow(inst), unique(inst$p_threshold)[1], mean(Fstat), min(Fstat)))

  out <- tryCatch(extract_outcome_data(snps = inst$SNP, outcomes = GC),
                  error = function(e) { message("  outcome err: ", e$message); NULL })
  if (is.null(out) || nrow(out) == 0) { cat("  NO outcome SNPs — skipping.\n"); next }
  harm <- harmonise_data(inst, out, action = 2)
  harm <- harm[harm$mr_keep, , drop = FALSE]
  if (nrow(harm) < 2) { cat("  <2 SNPs after harmonise — skipping.\n"); next }
  cat(sprintf("  harmonised SNPs used: %d\n", nrow(harm)))

  res <- mr(harm)                                   # IVW, Egger, WM, mode, Wald
  res$exposure <- elab; res$axis <- exposures$axis[i]
  res$mean_F <- mean(Fstat); res$n_instruments <- nrow(harm)
  oo <- generate_odds_ratios(res)
  print(oo[, c("method","nsnp","b","se","pval","or","or_lci95","or_uci95")])
  all_res[[elab]] <- oo

  het  <- tryCatch(mr_heterogeneity(harm),  error=function(e) NULL)
  plei <- tryCatch(mr_pleiotropy_test(harm), error=function(e) NULL)
  st   <- tryCatch(directionality_test(harm), error=function(e) NULL)
  if (!is.null(het))  { het$exposure  <- elab; all_het[[elab]]  <- het }
  if (!is.null(plei)) { plei$exposure <- elab; all_plei[[elab]] <- plei }
  if (!is.null(st))   { st$exposure   <- elab; all_steiger[[elab]] <- st }

  # per-exposure plots + leave-one-out
  tryCatch({
    p1 <- mr_scatter_plot(res, harm)
    ggplot2::ggsave(file.path(OUT, sprintf("scatter_%s.png", gsub("[^A-Za-z0-9]","_",elab))),
                    p1[[1]], width=6, height=5, dpi=150)
    loo <- mr_leaveoneout(harm)
    p2 <- mr_leaveoneout_plot(loo)
    ggplot2::ggsave(file.path(OUT, sprintf("loo_%s.png", gsub("[^A-Za-z0-9]","_",elab))),
                    p2[[1]], width=6, height=5, dpi=150)
  }, error = function(e) message("  plot err: ", e$message))

  # MR-PRESSO if enough SNPs and package available
  if (nrow(harm) >= 4 && requireNamespace("MRPRESSO", quietly = TRUE)) {
    pr <- tryCatch(run_mr_presso(harm, NbDistribution = 1000), error=function(e) NULL)
    if (!is.null(pr)) saveRDS(pr, file.path(OUT, sprintf("presso_%s.rds", gsub("[^A-Za-z0-9]","_",elab))))
  }
}

if (length(all_res)) {
  summ <- do.call(rbind, all_res)
  write.csv(summ, file.path(OUT, "MR_results_all_methods_REAL.csv"), row.names = FALSE)
  if (length(all_het))  write.csv(do.call(rbind, all_het),  file.path(OUT, "MR_heterogeneity_REAL.csv"), row.names=FALSE)
  if (length(all_plei)) write.csv(do.call(rbind, all_plei), file.path(OUT, "MR_pleiotropy_REAL.csv"), row.names=FALSE)
  if (length(all_steiger)) write.csv(do.call(rbind, all_steiger), file.path(OUT, "MR_steiger_REAL.csv"), row.names=FALSE)
  cat("\n================ REAL MR SUMMARY (IVW rows) ================\n")
  ivw <- summ[summ$method %in% c("Inverse variance weighted","Wald ratio"), ]
  print(ivw[, c("exposure","axis","method","nsnp","or","or_lci95","or_uci95","pval","mean_F")])
  cat("\nSaved to ", OUT, "\n")
} else {
  cat("\nNo exposures yielded usable instruments. Nothing written.\n")
}

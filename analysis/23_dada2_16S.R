#!/usr/bin/env Rscript
# =============================================================================
#  23_dada2_16S.R  --  REAL 16S ASV inference from raw DDBJ PRJDB20660 FASTQ
#
#  Replaces the previously fabricated (multinomial-sampled) OTU table with a
#  genuine DADA2 ASV table denoised from the raw 2x301bp MiSeq V3-V4 reads.
#
#  Region : V3-V4, primers 341F (CCTACGGGNGGCWGCAG, 17nt) / 805R
#           (GACTACHVGGGTATCTAATCC, 21nt) -> trimLeft = c(17, 21)
#  Input  : data/microbiome/raw/DRR*_{1,2}.fastq.bz2  (944 samples)
#  Manifest (run -> phenotype/patient): data/microbiome/reprocess/manifest.csv
#  Validation: per-sample read tracking is compared to the published
#              Supplementary Table 3 (input/filtered/denoised/merged/non-chim);
#              a high correlation proves both the pipeline and the DDBJ<->supp
#              sample mapping are correct.
#
#  Env vars:
#    PILOT_N   : if set (e.g. 16), process only this many samples spanning all
#                phenotypes (for parameter validation before the full run).
#    TRUNC_F / TRUNC_R : override forward/reverse truncation lengths.
#    DADA_THREADS : threads (default: detectCores()-2)
# =============================================================================
suppressMessages({ library(dada2) })
set.seed(1105)

RAW  <- "data/microbiome/raw"
OUT  <- "data/microbiome/dada2"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
FILT <- file.path(OUT, "filtered"); dir.create(FILT, showWarnings=FALSE)
RD   <- "results/rdata"; dir.create(RD, showWarnings=FALSE)
man  <- read.csv("data/microbiome/reprocess/manifest.csv", stringsAsFactors=FALSE)

THREADS <- as.integer(Sys.getenv("DADA_THREADS", unset=NA))
if (is.na(THREADS)) THREADS <- max(1, parallel::detectCores() - 2)
truncF <- as.integer(Sys.getenv("TRUNC_F", unset="260"))
truncR <- as.integer(Sys.getenv("TRUNC_R", unset="220"))

# --- discover FASTQ pairs actually present ----------------------------------
fnFs <- sort(list.files(RAW, pattern="_1\\.fastq\\.bz2$", full.names=TRUE))
fnRs <- sort(list.files(RAW, pattern="_2\\.fastq\\.bz2$", full.names=TRUE))
runF <- sub("_1\\.fastq\\.bz2$","",basename(fnFs))
runR <- sub("_2\\.fastq\\.bz2$","",basename(fnRs))
common <- intersect(runF, runR)
fnFs <- fnFs[match(common, runF)]; fnRs <- fnRs[match(common, runR)]
runs <- common
cat(sprintf("[DISCOVER] %d paired samples on disk (manifest lists %d)\n",
            length(runs), nrow(man)))

# --- optional pilot subset spanning phenotypes ------------------------------
pilotN <- suppressWarnings(as.integer(Sys.getenv("PILOT_N", unset="")))
if (!is.na(pilotN) && pilotN > 0) {
  ph <- man$phenotype[match(runs, man$run)]
  keep <- unlist(lapply(split(seq_along(runs), ph),
                        function(ix) head(ix, ceiling(pilotN/length(unique(ph))))))
  keep <- sort(keep)[seq_len(min(pilotN, length(keep)))]
  runs <- runs[keep]; fnFs <- fnFs[keep]; fnRs <- fnRs[keep]
  cat(sprintf("[PILOT] processing %d samples: %s\n", length(runs),
              paste(table(ph[keep]), names(table(ph[keep])), collapse=", ")))
  OUT_TRACK <- file.path(OUT, sprintf("track_pilot_%d.csv", pilotN))
} else {
  OUT_TRACK <- file.path(OUT, "track_reads.csv")
}

filtFs <- file.path(FILT, paste0(runs, "_F_filt.fastq.gz"))
filtRs <- file.path(FILT, paste0(runs, "_R_filt.fastq.gz"))
names(filtFs) <- runs; names(filtRs) <- runs

# --- 1. filter & trim (remove primers via trimLeft) -------------------------
cat(sprintf("[FILTER] trimLeft=c(17,21) truncLen=c(%d,%d) maxEE=c(2,2)\n", truncF, truncR))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     trimLeft=c(17,21), truncLen=c(truncF,truncR),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=THREADS, verbose=TRUE)
rownames(out) <- runs
# drop samples that lost all reads
ok <- file.exists(filtFs)
filtFs <- filtFs[ok]; filtRs <- filtRs[ok]; runsOk <- runs[ok]
cat(sprintf("[FILTER] %d/%d samples retained after filtering\n", sum(ok), length(runs)))

# --- 2. learn errors --------------------------------------------------------
cat("[ERRORS] learnErrors (forward, reverse)...\n")
errF <- learnErrors(filtFs, multithread=THREADS, nbases=1e8, verbose=1)
errR <- learnErrors(filtRs, multithread=THREADS, nbases=1e8, verbose=1)

# --- 3. denoise (pseudo-pooling for sensitivity at scale) -------------------
cat("[DADA] denoising with pool='pseudo'...\n")
ddF <- dada(filtFs, err=errF, multithread=THREADS, pool="pseudo", verbose=0)
ddR <- dada(filtRs, err=errR, multithread=THREADS, pool="pseudo", verbose=0)

# --- 4. merge + sequence table + chimera removal ----------------------------
cat("[MERGE] mergePairs (minOverlap=12)...\n")
mg <- mergePairs(ddF, filtFs, ddR, filtRs, minOverlap=12, maxMismatch=0, verbose=TRUE)
seqtab <- makeSequenceTable(mg)
cat(sprintf("[SEQTAB] %d samples x %d ASVs; length dist:\n", nrow(seqtab), ncol(seqtab)))
print(table(nchar(getSequences(seqtab))))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=THREADS, verbose=TRUE)
cat(sprintf("[CHIMERA] %d ASVs kept (%.1f%% of reads)\n",
            ncol(seqtab.nochim), 100*sum(seqtab.nochim)/sum(seqtab)))

# --- 5. track reads through the pipeline ------------------------------------
getN <- function(x) sum(getUniques(x))
track <- data.frame(
  run      = runsOk,
  input    = out[runsOk, 1],
  filtered = out[runsOk, 2],
  denoisedF= sapply(ddF, getN),
  denoisedR= sapply(ddR, getN),
  merged   = sapply(mg, getN),
  nonchim  = rowSums(seqtab.nochim)[runsOk],
  stringsAsFactors=FALSE)
track$phenotype <- man$phenotype[match(track$run, man$run)]
track$patient   <- man$patient[match(track$run, man$run)]
write.csv(track, OUT_TRACK, row.names=FALSE)
cat(sprintf("[TRACK] wrote %s\n", OUT_TRACK))
print(head(track))

# --- 6. assign taxonomy (SILVA) ---------------------------------------------
silva <- "data/microbiome/ref/silva_nr99_v138.1_train_set.fa.gz"
if (file.exists(silva)) {
  cat("[TAXONOMY] assignTaxonomy against SILVA v138.1...\n")
  tax <- assignTaxonomy(seqtab.nochim, silva, multithread=THREADS,
                        tryRC=TRUE, verbose=TRUE)
  sp <- "data/microbiome/ref/silva_species_assignment_v138.1.fa.gz"
  if (file.exists(sp)) tax <- addSpecies(tax, sp)
} else {
  cat("[TAXONOMY] SILVA ref not found at", silva, "- skipping taxonomy this run.\n")
  tax <- NULL
}

# --- 7. save ----------------------------------------------------------------
save(seqtab.nochim, tax, track, out, file=file.path(RD, "dada2_16S.RData"))
# ASV table with short IDs
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))
otu <- t(seqtab.nochim); rownames(otu) <- asv_ids
write.csv(data.frame(ASV=asv_ids, sequence=asv_seqs),
          file.path(OUT, "asv_sequences.csv"), row.names=FALSE)
write.csv(otu, file.path(OUT, "asv_table.csv"))
if (!is.null(tax)) { rownames(tax) <- asv_ids
  write.csv(tax, file.path(OUT, "asv_taxonomy.csv")) }
cat("[DONE] DADA2 pipeline complete.\n")

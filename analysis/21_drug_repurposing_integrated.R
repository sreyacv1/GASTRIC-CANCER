#!/usr/bin/env Rscript
# 21_drug_repurposing_integrated.R
# Connectivity-map style drug repurposing for gastric cancer, using the
# INTEGRATED TCGA+GTEx tumour-vs-normal DEG signature (script 16 used TCGA-only).
# Logic identical to 16: a candidate reverser is a compound whose OWN perturbation
# signature is ANTI-correlated with the GC tumour signature, i.e. it DOWN-regulates
# the tumour-UP genes and/or UP-regulates the tumour-DOWN genes.
# Real Enrichr API results only. Hypothesis-generating, in-silico, no validation.

suppressMessages({library(enrichR); library(ggplot2)})
set.seed(1)

OUT <- "results/drug_repurposing_integrated"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- file.path(OUT, "raw"); dir.create(RAW, showWarnings = FALSE)

# ---- 1. Define tumour signature from INTEGRATED DEG ------------------------
# Integrated matrix is scale-standardised, so logFC is compressed: rank by the
# moderated t-statistic instead. UP = top 150 by most-positive t (adj.P<0.05),
# DOWN = top 150 by most-negative t (adj.P<0.05).
deg <- read.csv("results/tables/DEG_integrated_TCGA_GTEx.csv", stringsAsFactors = FALSE)
deg <- deg[!is.na(deg$adj.P.Val) & !is.na(deg$t) & nzchar(deg$gene), ]

up <- deg[deg$adj.P.Val < 0.05 & deg$t > 0, ]
dn <- deg[deg$adj.P.Val < 0.05 & deg$t < 0, ]
up <- head(up[order(-up$t), ], 150)
dn <- head(dn[order( dn$t), ], 150)
up_genes <- unique(up$gene)
dn_genes <- unique(dn$gene)
cat(sprintf("Integrated signature: %d UP genes, %d DOWN genes\n",
            length(up_genes), length(dn_genes)))

setEnrichrSite("Enrichr")

# ---- drug-name parsers (per library term format) --------------------------
parse_lincs <- function(t) {
  x <- strsplit(t, "-", fixed = TRUE)[[1]]
  if (length(x) < 3) return(tolower(trimws(x[length(x)])))
  tolower(trimws(paste(x[2:(length(x) - 1)], collapse = "-")))
}
parse_consensus <- function(t) tolower(trimws(sub("\\s+(Up|Down)$", "", t)))
parse_first    <- function(t) tolower(trimws(sub("\\s.*$", "", t)))

# ---- 2/3. Query reverser arms ---------------------------------------------
run <- function(genes, db, arm, parser, direction) {
  r <- tryCatch(enrichr(genes, db)[[1]], error = function(e) {message("  FAIL ", db, ": ", conditionMessage(e)); NULL})
  if (is.null(r) || nrow(r) == 0) return(NULL)
  write.csv(r, file.path(RAW, sprintf("%s__arm-%s.csv", db, arm)), row.names = FALSE)
  r$library <- db; r$arm <- arm; r$drug_regulation <- direction
  r$drug <- vapply(r$Term, parser, character(1))
  r[, c("drug","library","arm","drug_regulation","Term","Overlap",
        "P.value","Adjusted.P.value","Combined.Score","Genes")]
}

jobs <- list(
  list(up_genes, "LINCS_L1000_Chem_Pert_down",             "UP",   parse_lincs,     "Down"),
  list(up_genes, "Drug_Perturbations_from_GEO_down",       "UP",   parse_first,     "Down"),
  list(up_genes, "LINCS_L1000_Chem_Pert_Consensus_Sigs",   "UP",   parse_consensus, "consensus"),
  list(dn_genes, "LINCS_L1000_Chem_Pert_up",               "DOWN", parse_lincs,     "Up"),
  list(dn_genes, "Drug_Perturbations_from_GEO_up",         "DOWN", parse_first,     "Up"),
  list(dn_genes, "LINCS_L1000_Chem_Pert_Consensus_Sigs",   "DOWN", parse_consensus, "consensus")
)

res <- do.call(rbind, lapply(jobs, function(j) run(j[[1]], j[[2]], j[[3]], j[[4]], j[[5]])))

keep_consensus <- with(res, drug_regulation != "consensus" |
  (arm == "UP"   & grepl("\\bDown$", Term)) |
  (arm == "DOWN" & grepl("\\bUp$",   Term)))
res <- res[keep_consensus, ]
res$drug_regulation[res$drug_regulation == "consensus"] <-
  ifelse(res$arm[res$drug_regulation == "consensus"] == "UP", "Down", "Up")

sig <- res[res$Adjusted.P.value < 0.05 & nzchar(res$drug), ]
write.csv(res, file.path(OUT, "all_reverser_enrichment.csv"), row.names = FALSE)

ds <- rbind(run(up_genes, "DSigDB", "UP",   parse_first, "assoc"),
            run(dn_genes, "DSigDB", "DOWN", parse_first, "assoc"))
if (!is.null(ds)) write.csv(ds, file.path(OUT, "DSigDB_supporting.csv"), row.names = FALSE)

# ---- 4. Rank candidate drugs ----------------------------------------------
agg <- do.call(rbind, lapply(split(sig, sig$drug), function(d) {
  data.frame(
    drug        = d$drug[1],
    n_arms      = length(unique(d$arm)),
    both_arms   = length(unique(d$arm)) == 2,
    n_signatures= nrow(d),
    min_adj_p   = min(d$Adjusted.P.value),
    arms        = paste(sort(unique(d$arm)), collapse = "+"),
    libraries   = paste(sort(unique(d$library)), collapse = "; "),
    max_combined= max(d$Combined.Score),
    stringsAsFactors = FALSE)
}))
agg <- agg[order(-agg$n_arms, agg$min_adj_p, -agg$max_combined), ]
agg$rank <- seq_len(nrow(agg))
write.csv(agg, file.path(OUT, "candidate_drugs_ranked.csv"), row.names = FALSE)

byarm <- function(a, n) {
  d <- agg[grepl(a, agg$arms), ]
  head(d[order(d$min_adj_p), ], n)
}
top <- rbind(byarm("UP", 8), byarm("DOWN", 8))
top <- top[!duplicated(top$drug), ]

cat("\n=== TOP CANDIDATE REVERSING COMPOUNDS (integrated DEG, balanced per arm) ===\n")
print(top[order(top$arms, top$min_adj_p),
          c("drug","arms","n_signatures","min_adj_p","libraries")], row.names = FALSE)
cat("\nBoth-arm (n_arms==2) candidates:", sum(agg$n_arms == 2),
    "  |  UP-arm sig drugs:", sum(grepl("UP", agg$arms)),
    "  DOWN-arm sig drugs:", sum(grepl("DOWN", agg$arms)), "\n")

# ---- barplot --------------------------------------------------------------
top$label <- factor(top$drug, levels = rev(top$drug[order(top$arms, top$min_adj_p)]))
p <- ggplot(top, aes(label, -log10(min_adj_p), fill = arms)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("UP" = "#4C72B0", "DOWN" = "#DD8452",
                               "DOWN+UP" = "#55A868"),
                    name = "Reverser arm",
                    labels = c("UP" = "represses tumour-UP genes",
                               "DOWN" = "induces tumour-DOWN genes",
                               "DOWN+UP" = "both arms")) +
  labs(title = "Candidate GC-signature-reversing compounds (integrated TCGA+GTEx DEG)",
       subtitle = "In-silico, hypothesis-generating; no experimental validation",
       x = NULL, y = expression(-log[10]~"(best adjusted p)")) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "top_candidate_drugs.png"), p, width = 9, height = 6, dpi = 150)

cat("\nOutputs written to", OUT, "\n")

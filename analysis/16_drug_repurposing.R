#!/usr/bin/env Rscript
# 16_drug_repurposing.R
# Connectivity-map style drug repurposing for gastric cancer.
# Logic: a candidate reverser is a compound whose OWN perturbation signature is
# ANTI-correlated with the GC tumour signature, i.e. it DOWN-regulates the
# tumour-UP genes and/or UP-regulates the tumour-DOWN genes.
# Real Enrichr API results only. Hypothesis-generating, in-silico, no validation.

suppressMessages({library(enrichR); library(ggplot2)})
set.seed(1)

OUT <- "results/drug_repurposing"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- file.path(OUT, "raw"); dir.create(RAW, showWarnings = FALSE)

# ---- 1. Define tumour signature -------------------------------------------
deg <- read.csv("results/tables/TCGA_DEG_results_symbols.csv", stringsAsFactors = FALSE)
deg <- deg[!is.na(deg$padj) & !is.na(deg$log2FoldChange) & nzchar(deg$gene_symbol), ]

up   <- deg[deg$padj < 0.05 & deg$log2FoldChange >  1, ]
dn   <- deg[deg$padj < 0.05 & deg$log2FoldChange < -1, ]
up   <- head(up[order(-up$log2FoldChange), ], 150)
dn   <- head(dn[order( dn$log2FoldChange), ], 150)
up_genes <- unique(up$gene_symbol)
dn_genes <- unique(dn$gene_symbol)
cat(sprintf("Signature: %d UP genes, %d DOWN genes\n", length(up_genes), length(dn_genes)))

setEnrichrSite("Enrichr")

# ---- drug-name parsers (per library term format) --------------------------
# ponytail: heuristic name extraction, tuned to the 4 term formats probed live.
# Ceiling: drug names containing hyphens in the per-signature LINCS lib parse
# imperfectly; consensus + GEO libs give clean names and carry the ranking.
parse_lincs <- function(t) {           # "<batch cell time>-<drug>-<dose>", drug may hold hyphens
  x <- strsplit(t, "-", fixed = TRUE)[[1]]
  if (length(x) < 3) return(tolower(trimws(x[length(x)])))
  tolower(trimws(paste(x[2:(length(x) - 1)], collapse = "-")))  # drop cell prefix + dose suffix
}
parse_consensus <- function(t) tolower(trimws(sub("\\s+(Up|Down)$", "", t)))  # "Palbociclib Down"
parse_first    <- function(t) tolower(trimws(sub("\\s.*$", "", t)))           # GEO / DSigDB first token

# ---- 2/3. Query reverser arms ---------------------------------------------
# arm "UP"   = tumour-UP genes vs drug-DOWN libraries  (drug represses tumour-UP)
# arm "DOWN" = tumour-DOWN genes vs drug-UP libraries  (drug induces tumour-DOWN)
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
  # tumour-UP genes vs drug-DOWN signatures
  list(up_genes, "LINCS_L1000_Chem_Pert_down",             "UP",   parse_lincs,     "Down"),
  list(up_genes, "Drug_Perturbations_from_GEO_down",       "UP",   parse_first,     "Down"),
  list(up_genes, "LINCS_L1000_Chem_Pert_Consensus_Sigs",   "UP",   parse_consensus, "consensus"),
  # tumour-DOWN genes vs drug-UP signatures
  list(dn_genes, "LINCS_L1000_Chem_Pert_up",               "DOWN", parse_lincs,     "Up"),
  list(dn_genes, "Drug_Perturbations_from_GEO_up",         "DOWN", parse_first,     "Up"),
  list(dn_genes, "LINCS_L1000_Chem_Pert_Consensus_Sigs",   "DOWN", parse_consensus, "consensus")
)

res <- do.call(rbind, lapply(jobs, function(j) run(j[[1]], j[[2]], j[[3]], j[[4]], j[[5]])))

# Consensus lib is directional inside the Term; keep only reverser-direction rows.
# arm UP wants drug "Down" consensus sigs; arm DOWN wants drug "Up" consensus sigs.
keep_consensus <- with(res, drug_regulation != "consensus" |
  (arm == "UP"   & grepl("\\bDown$", Term)) |
  (arm == "DOWN" & grepl("\\bUp$",   Term)))
res <- res[keep_consensus, ]
res$drug_regulation[res$drug_regulation == "consensus"] <-
  ifelse(res$arm[res$drug_regulation == "consensus"] == "UP", "Down", "Up")

# significance filter for the ranking pool
sig <- res[res$Adjusted.P.value < 0.05 & nzchar(res$drug), ]
write.csv(res, file.path(OUT, "all_reverser_enrichment.csv"), row.names = FALSE)

# DSigDB: non-directional supporting evidence (drug-associated gene sets), both lists
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
# Prioritise both-arm reversers, then strongest adjusted p (pooled ranking).
agg <- agg[order(-agg$n_arms, agg$min_adj_p, -agg$max_combined), ]
agg$rank <- seq_len(nrow(agg))
write.csv(agg, file.path(OUT, "candidate_drugs_ranked.csv"), row.names = FALSE)

# The two arms sit on very different p-value scales and surface different
# biology, so pooled ranking lets one arm dominate. Present a BALANCED top set:
# strongest candidates from each arm, so both are represented honestly.
byarm <- function(a, n) {
  d <- agg[grepl(a, agg$arms), ]
  head(d[order(d$min_adj_p), ], n)
}
top <- rbind(byarm("UP", 8), byarm("DOWN", 8))
top <- top[!duplicated(top$drug), ]

cat("\n=== TOP CANDIDATE REVERSING COMPOUNDS (balanced per arm) ===\n")
print(top[order(top$arms, top$min_adj_p),
          c("drug","arms","n_signatures","min_adj_p","libraries")], row.names = FALSE)
cat("\nBoth-arm (n_arms==2) candidates:", sum(agg$n_arms == 2),
    "  |  UP-arm sig drugs:", sum(grepl("UP", agg$arms)),
    "  DOWN-arm sig drugs:", sum(grepl("DOWN", agg$arms)), "\n")

# ---- barplot --------------------------------------------------------------
top$label <- factor(top$drug, levels = rev(top$drug[order(top$arms, top$min_adj_p)]))
p <- ggplot(top, aes(label, -log10(min_adj_p), fill = arms)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("UP" = "#4C72B0", "DOWN" = "#DD8452"),
                    name = "Reverser arm",
                    labels = c("UP" = "represses tumour-UP genes",
                               "DOWN" = "induces tumour-DOWN genes")) +
  labs(title = "Candidate GC-signature-reversing compounds (Enrichr)",
       subtitle = "In-silico, hypothesis-generating; no experimental validation",
       x = NULL, y = expression(-log[10]~"(best adjusted p)")) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "top_candidate_drugs.png"), p, width = 9, height = 6, dpi = 150)

cat("\nOutputs written to", OUT, "\n")

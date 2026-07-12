###############################################################################
##                                                                           ##
##   GASTRIC CANCER MULTI-OMICS PIPELINE  v2.0                             ##
##   Microbiome Dysbiosis → Host Gene Expression → Causal Inference (MR)   ##
##                                                                           ##
##  NOVEL BIOLOGY:                                                           ##
##   • H. pylori + Streptococcus co-abundance synergy scoring               ##
##   • Lauren subtype-specific EMT drivers (Diffuse vs Intestinal)          ##
##   • Fanconi anemia / DNA repair pathway as microbiome-modulated hub      ##
##   • Calcium signaling as Streptococcus-responsive axis                   ##
##   • Immune microenvironment deconvolution (ESTIMATE)                     ##
##   • Bidirectional & Multivariable MR for causal architecture             ##
##                                                                           ##
##  DATASETS (ALL AUTO-DOWNLOADED):                                          ##
##   Transcriptome : TCGA-STAD, GTEx Stomach, GSE27342, GSE63089, GSE62254  ##
##   Microbiome    : DDBJ PRJDB20660, NCBI SRA PRJNA830774                 ##
##   MR GWAS       : MiBioGen taxa GWAS + FinnGen/UKBB gastric cancer       ##
##                                                                           ##
##  PIPELINE STEPS:                                                          ##
##   00. Configuration, Paths, Helpers                                       ##
##   01. Package Management                                                   ##
##   02. TCGA-STAD RNA-seq (Auto-Download + DEG + Subtypes)                 ##
##   03. GTEx Normal Stomach (Auto-Download)                                 ##
##   04. GEO Datasets (GSE27342, GSE63089, GSE62254 — Auto-Download)        ##
##   05. Multi-Cohort Harmonization & ComBat Batch Correction               ##
##   [continued in Part 2: Steps 06-11]                                     ##
##                                                                           ##
##  USAGE: source("gastric_cancer_multiomics_v2_part1.R")                   ##
##         then source("gastric_cancer_multiomics_v2_part2.R")              ##
##                                                                           ##
##  SYSTEM REQUIREMENTS:                                                     ##
##   R >= 4.3.0 | RAM >= 32 GB recommended | Disk >= 20 GB for raw data     ##
###############################################################################


###############################################################################
## 00. CONFIGURATION, PATHS & UTILITY HELPERS
###############################################################################

# ── 0a. Working directory ─────────────────────────────────────────────────────
# Change this path to your project root
# BASE_DIR <- "gastric_cancer_project"
# setwd(file.path(path.expand("~"), BASE_DIR))

# ── 0b. Directory scaffold ────────────────────────────────────────────────────
dirs <- c(
  "data/tcga", "data/gtex", "data/geo", "data/microbiome",
  "data/gwas", "data/processed",
  "results/plots/qc", "results/plots/transcriptome",
  "results/plots/microbiome", "results/plots/integration",
  "results/plots/wgcna", "results/plots/mr",
  "results/tables", "results/rdata", "logs"
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ── 0c. Logging ───────────────────────────────────────────────────────────────
LOG_FILE <- file.path("logs", paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M"), ".log"))

log_msg <- function(msg, level = "INFO") {
  ts      <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  full    <- sprintf("%s [%s] %s", ts, level, msg)
  cat(full, "\n")
  cat(full, "\n", file = LOG_FILE, append = TRUE)
}

# ── 0d. Checkpoint system: skip completed steps on re-run ─────────────────────
# Each major step saves a .done sentinel file
checkpoint_done <- function(step_name) {
  file.exists(file.path("results/rdata", paste0(step_name, ".done")))
}
mark_done <- function(step_name) {
  writeLines(as.character(Sys.time()),
             file.path("results/rdata", paste0(step_name, ".done")))
}

# ── 0e. Global settings ───────────────────────────────────────────────────────
options(
  stringsAsFactors  = FALSE,
  timeout           = 600,      # 10-min download timeout
  warn              = 1,        # print warnings immediately
  scipen             = 999
)
set.seed(42)

# ── 0f. Parallel processing setup ──────────────────────────────────────────────
library(BiocParallel)
N_CORES <- max(4, min(64, parallel::detectCores() - 4))
register(MulticoreParam(workers = N_CORES))
log_msg(sprintf("Parallel processing enabled with %d workers.", N_CORES))

# ── 0g. Publication-quality ggplot2 theme ─────────────────────────────────────
# Applied uniformly to all plots for consistency
theme_gc <- function(base_size = 12) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey40", size = base_size - 1),
      axis.title       = element_text(face = "bold"),
      legend.position  = "right",
      legend.key.size  = unit(0.45, "cm"),
      legend.text      = element_text(size = base_size - 2),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      panel.grid.major = element_line(colour = "grey94"),
      plot.margin      = margin(8, 8, 8, 8)
    )
}

# ── 0g. Standard colour palettes ──────────────────────────────────────────────
COL_STATUS  <- c(Tumor   = "#D7191C", Normal  = "#2C7BB6")
COL_DATASET <- c(TCGA    = "#1B7837", GTEx    = "#762A83",
                 GSE27342 = "#E08214", GSE63089 = "#4DAC26",
                 GSE62254 = "#D01C8B")
COL_LAUREN  <- c(Intestinal = "#F4A582", Diffuse = "#92C5DE",
                 Mixed      = "#A6DBA0", Unknown = "grey70")
COL_SUBTYPE <- c(EBV  = "#8073AC", MSI  = "#E08214",
                 GS   = "#1B7837", CIN  = "#D7191C")

log_msg("Configuration complete. Pipeline v2.0 initialized.")


###############################################################################
## 01. PACKAGE MANAGEMENT
###############################################################################
log_msg("=== STEP 01: Package Management ===")

# ── 1a-1d. Bypassed installation loop
log_msg("✓ Package setup assumed complete (bypassing loop).\n")

# ── 1e. Load ONLY essential packages first (save memory) ──────────────────────
pkgs_to_load <- c(
  "SummarizedExperiment", "DESeq2", "edgeR", "limma", "sva", "GSVA", "ashr",
  "biomaRt", "tidyverse", "data.table", "ggrepel", "patchwork", "pheatmap", "GEOquery"
)

suppressPackageStartupMessages({
  loaded <- sapply(pkgs_to_load, function(p) {
    tryCatch({ library(p, character.only = TRUE); TRUE },
             error = function(e) { log_msg(sprintf("Failed to load %s", p), "WARN"); FALSE })
  })
})

log_msg(sprintf("Packages loaded: %d/%d", sum(loaded), length(loaded)))
log_msg("✓ Package setup complete.\n")


###############################################################################
## 02. TCGA-STAD RNA-seq: DOWNLOAD, PREPROCESS, DEG, LAUREN SUBTYPES
###############################################################################
log_msg("=== STEP 02: TCGA-STAD ===")

TCGA_RDATA <- "data/processed/TCGA_STAD_processed.RData"

if (!checkpoint_done("step02_tcga")) {

  # ── 2a. Use local RDS (bypass GDC query) ───────────────────────────────────
  log_msg("  Loading TCGA-STAD from local RDS file…")
  rse_file <- "data/host/tcga_stad_rse.rds"
  if (file.exists(rse_file)) {
    tcga_se <- readRDS(rse_file)
  } else {
    stop("TCGA RDS file not found in data/host/!")
  }

  tcga_raw   <- assay(tcga_se, "unstranded")
  tcga_meta  <- as.data.frame(colData(tcga_se))

  # ── 2b. Sample classification ───────────────────────────────────────────────
  tcga_meta$sample_type_simple <- dplyr::case_when(
    grepl("Primary Tumor",       tcga_meta$sample_type) ~ "Tumor",
    grepl("Solid Tissue Normal", tcga_meta$sample_type) ~ "Normal",
    TRUE ~ "Other"
  )
  keep_idx  <- tcga_meta$sample_type_simple %in% c("Tumor", "Normal")
  tcga_raw  <- tcga_raw[, keep_idx]
  tcga_meta <- tcga_meta[keep_idx, ]

  # ── 2c. Lauren subtype extraction ───────────────────────────────────────────
  # TCGA STAD has Lauren classification in paper_Lauren.Class
  lauren_col <- grep("Lauren|lauren", colnames(tcga_meta), value = TRUE, ignore.case = TRUE)
  if (length(lauren_col) > 0) {
    tcga_meta$Lauren <- tcga_meta[[lauren_col[1]]]
    # Normalize Lauren values
    tcga_meta$Lauren <- dplyr::case_when(
      grepl("Diffuse",    tcga_meta$Lauren, ignore.case = TRUE) ~ "Diffuse",
      grepl("Intestinal", tcga_meta$Lauren, ignore.case = TRUE) ~ "Intestinal",
      grepl("Mixed",      tcga_meta$Lauren, ignore.case = TRUE) ~ "Mixed",
      TRUE ~ "Unknown"
    )
    tcga_meta$Lauren[is.na(tcga_meta$Lauren)] <- "Unknown"
  } else {
    tcga_meta$Lauren <- "Unknown"
  }

  # TCGA molecular subtypes (EBV/MSI/GS/CIN)
  subtype_col <- grep("Subtype|subtype|molecular", colnames(tcga_meta), value = TRUE, ignore.case = TRUE)
  if (length(subtype_col) > 0) {
    tcga_meta$TCGA_subtype <- tcga_meta[[subtype_col[1]]]
  } else {
    tcga_meta$TCGA_subtype <- "Unknown"
  }

  log_msg(sprintf("  TCGA-STAD: %d genes × %d samples (%d Tumor, %d Normal)",
                  nrow(tcga_raw), ncol(tcga_raw),
                  sum(tcga_meta$sample_type_simple == "Tumor"),
                  sum(tcga_meta$sample_type_simple == "Normal")))

  # ── 2d. Limma-Voom normalization (optimized for 8GB RAM) ────────────────────
  log_msg("  Running limma-voom normalization…")
  
  # Filter low counts (≥10 in ≥5% of samples)
  dge_tcga <- edgeR::DGEList(counts = round(tcga_raw), samples = tcga_meta)
  min_samps <- max(5, round(0.05 * ncol(dge_tcga)))
  keep_genes <- rowSums(edgeR::cpm(dge_tcga) >= 1) >= min_samps
  dge_tcga <- dge_tcga[keep_genes, , keep.lib.sizes = FALSE]
  dge_tcga <- edgeR::calcNormFactors(dge_tcga)
  
  design_tcga <- model.matrix(~ 0 + sample_type_simple, data = tcga_meta)
  colnames(design_tcga) <- gsub("sample_type_simple", "", colnames(design_tcga))
  
  # Voom transformation
  vst_tcga  <- limma::voom(dge_tcga, design_tcga, plot = FALSE)
  tcga_expr <- vst_tcga$E   # log2-CPM values
  
  # ── 2e. DEG analysis (Tumor vs Normal) using limma ──────────────────────────
  log_msg("  Extracting DEGs (Tumor vs Normal)…")
  fit_tcga  <- limma::lmFit(vst_tcga, design_tcga)
  cont_tcga <- limma::makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design_tcga)
  fit_tcga  <- limma::contrasts.fit(fit_tcga, cont_tcga)
  fit_tcga  <- limma::eBayes(fit_tcga)
  
  res_tcga <- limma::topTable(fit_tcga, coef = "Tumor_vs_Normal", number = Inf, sort.by = "P")
  
  res_tcga_df <- as.data.frame(res_tcga) %>%
    rownames_to_column("gene_id") %>%
    dplyr::rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val) %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(
      gene_id_clean = gsub("\\..*", "", gene_id),
      sig = dplyr::case_when(
        padj < 0.05 & log2FoldChange >  1 ~ "Up",
        padj < 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE                               ~ "NS"
      )
    )

  # ── 2f. Ensembl → Gene Symbol mapping ───────────────────────────────────────
  log_msg("  Mapping Ensembl IDs to gene symbols locally via org.Hs.eg.db…")
  gene_map <- NULL
  
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    ensembl_ids <- unique(res_tcga_df$gene_id_clean)
    sym_map <- tryCatch({
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys      = ensembl_ids,
        columns   = c("SYMBOL", "GENENAME"),
        keytype   = "ENSEMBL"
      )
    }, error = function(e) {
      log_msg(paste("  Local mapping error:", e$message), "WARN")
      NULL
    })
    
    if (!is.null(sym_map)) {
      gene_map <- sym_map %>%
        dplyr::rename(ensembl_gene_id = ENSEMBL, hgnc_symbol = SYMBOL, gene_biotype = GENENAME)
    }
  }
  
  # Fallback to BioMart only if local mapping failed or yielded nothing
  if (is.null(gene_map) || nrow(gene_map) == 0) {
    log_msg("  Local mapping failed. Querying BioMart as fallback…", "WARN")
    tryCatch({
      mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                      host = "https://www.ensembl.org")
      gene_map <- getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
        filters    = "ensembl_gene_id",
        values     = unique(res_tcga_df$gene_id_clean),
        mart       = mart
      )
    }, error = function(e) {
      log_msg("  BioMart primary failed, trying mirror…", "WARN")
      tryCatch({
        mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host = "https://useast.ensembl.org")
        gene_map <- getBM(
          attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
          filters    = "ensembl_gene_id",
          values     = unique(res_tcga_df$gene_id_clean),
          mart       = mart
        )
      }, error = function(e2) {
        log_msg("  BioMart failed entirely. Using Ensembl IDs as symbols.", "WARN")
      })
    })
  }

  if (is.null(gene_map) || nrow(gene_map) == 0) {
    gene_map <- data.frame(
      ensembl_gene_id = unique(res_tcga_df$gene_id_clean),
      hgnc_symbol     = unique(res_tcga_df$gene_id_clean),
      gene_biotype    = "unknown"
    )
  }

  gene_map <- gene_map %>%
    dplyr::filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
  
  res_tcga_df <- res_tcga_df %>%
    dplyr::left_join(gene_map, by = c("gene_id_clean" = "ensembl_gene_id"))
  
  write.csv(res_tcga_df, "results/tables/TCGA_DEG_results.csv", row.names = FALSE)

  # ── 2g. Convert TCGA expression to gene symbols ───────────────────────────
  rownames(tcga_expr) <- gsub("\\..*", "", rownames(tcga_expr))
  tcga_expr_sym <- tcga_expr[rownames(tcga_expr) %in% gene_map$ensembl_gene_id, ]
  id2sym <- setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)
  rownames(tcga_expr_sym) <- id2sym[rownames(tcga_expr_sym)]
  tcga_expr_sym <- tcga_expr_sym[!is.na(rownames(tcga_expr_sym)) &
                                   !duplicated(rownames(tcga_expr_sym)), ]

  # ── 2h. DEG count summary ────────────────────────────────────────────────────
  n_up   <- sum(res_tcga_df$sig == "Up",   na.rm = TRUE)
  n_down <- sum(res_tcga_df$sig == "Down", na.rm = TRUE)
  log_msg(sprintf("  DEGs: %d upregulated, %d downregulated (|LFC|>1, padj<0.05)",
                  n_up, n_down))

  # ── 2i. Volcano plot ──────────────────────────────────────────────────────────
  top_labels <- res_tcga_df %>%
    dplyr::filter(sig != "NS", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
    dplyr::group_by(sig) %>%
    dplyr::slice_min(padj, n = 12) %>%
    dplyr::ungroup()

  p_volcano <- ggplot(res_tcga_df, aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(alpha = 0.35, size = 0.7, stroke = 0) +
    geom_text_repel(data = top_labels, aes(label = hgnc_symbol),
                    size = 2.6, max.overlaps = 25, min.segment.length = 0.2,
                    fontface = "italic", colour = "black") +
    scale_colour_manual(values = c(Up = "#D7191C", Down = "#2C7BB6", NS = "grey75"),
                        labels = c(Up = sprintf("Up (%d)", n_up),
                                   Down = sprintf("Down (%d)", n_down),
                                   NS = "NS")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "black", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black", linewidth = 0.3) +
    scale_x_continuous(limits = c(-8, 8), oob = scales::squish) +
    labs(title    = "TCGA-STAD: Tumor vs Normal",
         subtitle = sprintf("n = %d genes tested | limma-voom with eBayes moderated t-test", nrow(res_tcga_df)),
         x        = expression(log[2]~Fold~Change),
         y        = expression(-log[10]~(adj.~p)),
         colour   = NULL) +
    theme_gc() +
    theme(legend.position = c(0.15, 0.90))

  ggsave("results/plots/transcriptome/TCGA_volcano.pdf", p_volcano,
         width = 8, height = 6, useDingbats = FALSE)

  # ── 2j. PCA coloured by sample type AND Lauren subtype ───────────────────────
  pca_tcga   <- prcomp(t(tcga_expr), scale. = TRUE)
  var_exp    <- round(summary(pca_tcga)$importance[2, 1:3] * 100, 1)
  pca_df     <- as.data.frame(pca_tcga$x[, 1:3]) %>%
    dplyr::mutate(sample_id = rownames(.)) %>%
    dplyr::left_join(
      tcga_meta %>% dplyr::mutate(sample_id = rownames(.)) %>%
        dplyr::select(sample_id, sample_type_simple, Lauren, TCGA_subtype),
      by = "sample_id"
    )

  p_pca_type <- ggplot(pca_df, aes(PC1, PC2, colour = sample_type_simple)) +
    geom_point(size = 1.8, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.5) +
    scale_colour_manual(values = COL_STATUS) +
    labs(title    = "TCGA-STAD PCA",
         subtitle = "Coloured by sample type",
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         colour = NULL) +
    theme_gc()

  p_pca_lauren <- ggplot(pca_df %>% dplyr::filter(sample_type_simple == "Tumor"),
                         aes(PC1, PC2, colour = Lauren)) +
    geom_point(size = 1.8, alpha = 0.8) +
    stat_ellipse(level = 0.90, linetype = "dashed", linewidth = 0.5, aes(group = Lauren)) +
    scale_colour_manual(values = COL_LAUREN) +
    labs(title    = "TCGA-STAD Tumors — Lauren Subtypes",
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         colour = "Lauren\nSubtype") +
    theme_gc()

  p_pca_subtype <- ggplot(pca_df %>% dplyr::filter(sample_type_simple == "Tumor",
                                                    TCGA_subtype != "Unknown"),
                          aes(PC1, PC2, colour = TCGA_subtype)) +
    geom_point(size = 1.8, alpha = 0.8) +
    stat_ellipse(level = 0.90, linetype = "dashed", linewidth = 0.5) +
    scale_colour_manual(values = COL_SUBTYPE) +
    labs(title    = "TCGA-STAD Tumors — Molecular Subtypes",
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         colour = "TCGA\nSubtype") +
    theme_gc()

  ggsave("results/plots/qc/TCGA_PCA_triptych.pdf",
         p_pca_type + p_pca_lauren + p_pca_subtype,
         width = 18, height = 5.5, useDingbats = FALSE)

  # ── 2k. Top-50 DEG heatmap (TCGA only) ───────────────────────────────────────
  top50_genes <- res_tcga_df %>%
    dplyr::filter(sig != "NS", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
    dplyr::slice_min(padj, n = 50) %>%
    dplyr::pull(hgnc_symbol)
  top50_genes <- intersect(top50_genes, rownames(tcga_expr_sym))

  ann_col_tcga <- data.frame(
    Status  = tcga_meta$sample_type_simple,
    Lauren  = tcga_meta$Lauren,
    Subtype = tcga_meta$TCGA_subtype,
    row.names = rownames(tcga_meta)
  )

  pdf("results/plots/transcriptome/TCGA_top50_DEG_heatmap.pdf", width = 14, height = 9)
  pheatmap(
    tcga_expr_sym[top50_genes, ],
    annotation_col    = ann_col_tcga,
    annotation_colors = list(
      Status  = COL_STATUS,
      Lauren  = COL_LAUREN,
      Subtype = COL_SUBTYPE
    ),
    show_colnames     = FALSE,
    scale             = "row",
    clustering_method = "ward.D2",
    color             = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
    main              = "Top 50 DEGs — TCGA-STAD (Tumor vs Normal)"
  )
  dev.off()

  # ── 2l. Lauren-subtype DEG (Diffuse vs Intestinal in tumors only) ─────────
  # Run Lauren subtype DEG analysis using limma-voom for methodological consistency
  log_msg("  Running Lauren subtype DEG analysis (limma-voom, within-tumor subgroup)…")
  tumor_meta_lauren <- tcga_meta %>%
    dplyr::filter(sample_type_simple == "Tumor",
                  Lauren %in% c("Diffuse", "Intestinal"))
  tumor_counts_lauren <- tcga_raw[, rownames(tumor_meta_lauren)]

  if (nrow(tumor_meta_lauren) >= 20 &&
      all(c("Diffuse", "Intestinal") %in% tumor_meta_lauren$Lauren)) {
    tumor_meta_lauren$Lauren <- factor(tumor_meta_lauren$Lauren, levels = c("Intestinal", "Diffuse"))
    design_lauren <- model.matrix(~ Lauren, data = tumor_meta_lauren)
    dge_lauren <- DGEList(counts = tumor_counts_lauren)
    dge_lauren <- calcNormFactors(dge_lauren)
    v_lauren <- voom(dge_lauren, design_lauren, plot = FALSE)
    fit_lauren <- lmFit(v_lauren, design_lauren)
    fit_lauren <- eBayes(fit_lauren)
    res_lauren_df <- topTable(fit_lauren, coef = 2, number = Inf) %>%
      rownames_to_column("gene_id") %>%
      dplyr::mutate(gene_id_clean = gsub("\\..*", "", gene_id)) %>%
      dplyr::left_join(gene_map, by = c("gene_id_clean" = "ensembl_gene_id")) %>%
      dplyr::filter(!is.na(adj.P.Val)) %>%
      dplyr::arrange(adj.P.Val)
    
    # Map column names to standard
    res_lauren_df$log2FoldChange <- res_lauren_df$logFC
    res_lauren_df$padj <- res_lauren_df$adj.P.Val

    write.csv(res_lauren_df, "results/tables/TCGA_Lauren_DEG_Diffuse_vs_Intestinal.csv", row.names = FALSE)
    log_msg(sprintf("  Lauren DEGs: %d up (Diffuse), %d down (Intestinal)",
                    sum(res_lauren_df$padj < 0.05 & res_lauren_df$log2FoldChange > 1, na.rm = TRUE),
                    sum(res_lauren_df$padj < 0.05 & res_lauren_df$log2FoldChange < -1, na.rm = TRUE)))
  } else {
    res_lauren_df <- NULL
    log_msg("  Insufficient samples for Lauren DEG analysis", "WARN")
  }

  # ── 2m. Save ─────────────────────────────────────────────────────────────────
  save(vst_tcga, tcga_expr, tcga_expr_sym, tcga_meta, tcga_raw,
       res_tcga_df, res_lauren_df, gene_map,
       file = TCGA_RDATA)
  mark_done("step02_tcga")
  log_msg("  ✓ TCGA-STAD processing complete.\n")

} else {
  log_msg("  Loading cached TCGA data…")
  load(TCGA_RDATA)
  log_msg("  ✓ TCGA cache loaded.\n")
}


###############################################################################
## 03. GTEx NORMAL STOMACH — AUTO-DOWNLOAD v10
###############################################################################
log_msg("=== STEP 03: GTEx Normal Stomach ===")

GTEX_RDATA <- "data/processed/GTEx_stomach.RData"

if (!checkpoint_done("step03_gtex")) {

  # ── 3a. Download GTEx v10 TPM ─────────────────────────────────────────────
  # GTEx v10 TPM matrix (hg38) — requires no login
  gtex_tpm_url  <- paste0(
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/",
    "rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"
  )
  gtex_attr_url <- paste0(
    "https://storage.googleapis.com/adult-gtex/annotations/v10/",
    "metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
  )

  gtex_tpm_path  <- "data/host/GTEx_v10_tpm.gct.gz"
  gtex_attr_path <- "data/host/GTEx_v10_attrs.txt"

  log_msg("  Using pre-existing GTEx v10 files from data/host/")

  if (file.exists(gtex_tpm_path) && file.exists(gtex_attr_path)) {

    log_msg("  Parsing GTEx attributes → identifying stomach samples…")
    gtex_attr <- fread(gtex_attr_path, sep = "\t", quote = "")
    stomach_ids <- gtex_attr[SMTSD == "Stomach", SAMPID]
    log_msg(sprintf("  GTEx stomach samples found: %d", length(stomach_ids)))

    log_msg("  Reading GTEx TPM matrix (this may take a few minutes)…")
    # GCT v1.3: line 1 = version, line 2 = dims, line 3 = header
    gtex_raw <- fread(cmd = sprintf("zcat %s", gtex_tpm_path),
                      skip = 2, header = TRUE, sep = "\t", data.table = FALSE)

    rownames(gtex_raw) <- gtex_raw$Name
    gene_ids_gtex      <- gtex_raw$Name
    gtex_raw           <- gtex_raw[, -(1:2)]

    keep_cols    <- intersect(stomach_ids, colnames(gtex_raw))
    gtex_stomach <- as.matrix(gtex_raw[, keep_cols])
    rownames(gtex_stomach) <- gene_ids_gtex
    gtex_stomach <- log2(gtex_stomach + 1)  # log2(TPM + 1)

    # Clean Ensembl IDs
    rownames(gtex_stomach) <- gsub("\\..*", "", rownames(gtex_stomach))

    # Convert to gene symbols
    gtex_sym_map <- gene_map %>%
      dplyr::filter(ensembl_gene_id %in% rownames(gtex_stomach)) %>%
      dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

    gtex_stomach_sym <- gtex_stomach[gtex_sym_map$ensembl_gene_id, ]
    rownames(gtex_stomach_sym) <- gtex_sym_map$hgnc_symbol
    gtex_stomach_sym <- gtex_stomach_sym[!duplicated(rownames(gtex_stomach_sym)), ]

    log_msg(sprintf("  GTEx Stomach: %d genes × %d samples", nrow(gtex_stomach_sym), ncol(gtex_stomach_sym)))

    # GTEx metadata for batch correction
    gtex_meta_stomach <- data.frame(
      sample_id          = colnames(gtex_stomach_sym),
      sample_type_simple = "Normal",
      dataset            = "GTEx",
      Lauren             = "Unknown",
      TCGA_subtype       = "Normal",
      row.names          = colnames(gtex_stomach_sym),
      stringsAsFactors   = FALSE
    )

    save(gtex_stomach_sym, gtex_meta_stomach,
         file = GTEX_RDATA)
    mark_done("step03_gtex")
    log_msg("  ✓ GTEx stomach data saved.\n")

  } else {
    log_msg("  GTEx files not available; creating placeholder.", "WARN")
    gtex_stomach_sym  <- NULL
    gtex_meta_stomach <- NULL
    save(gtex_stomach_sym, gtex_meta_stomach, file = GTEX_RDATA)
    mark_done("step03_gtex")
  }

} else {
  log_msg("  Loading cached GTEx data…")
  load(GTEX_RDATA)
  log_msg("  ✓ GTEx cache loaded.\n")
}


###############################################################################
## 04. GEO DATASETS — AUTO-DOWNLOAD: GSE27342, GSE63089, GSE62254
###############################################################################
log_msg("=== STEP 04: GEO Datasets ===")

GEO_RDATA <- "data/processed/GEO_processed.RData"

if (!checkpoint_done("step04_geo")) {

  # Helper: process an ExpressionSet to a gene-symbol matrix
  process_gse <- function(gse_id, probe2gene_fn = NULL) {
    # Check if a series matrix file exists locally in data/host/ first
    local_series_file <- file.path("data/host", paste0(gse_id, "_series_matrix.txt.gz"))
    if (file.exists(local_series_file)) {
      log_msg(sprintf("  Loading %s from local series matrix file…", gse_id))
      gse_list <- getGEO(filename = local_series_file)
      if (is.list(gse_list)) gse <- gse_list[[1]] else gse <- gse_list
    } else {
      # Use getGEO with destdir = "data/geo/". It will automatically load the cached
      # file from data/geo/ if it exists.
      log_msg(sprintf("  Loading %s from GEO (using cached file in data/geo/ if available)…", gse_id))
      gse_list <- tryCatch(
        getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE,
               destdir = "data/geo/"),
        error = function(e) { log_msg(sprintf("  %s download/load failed: %s", gse_id, e$message), "WARN"); NULL }
      )
      if (is.null(gse_list)) return(NULL)
      gse <- gse_list[[1]]
    }

    expr_mat <- exprs(gse)
    pheno    <- pData(gse)
    feat     <- fData(gse)

    # ── Map probes to gene symbols ────────────────────────────────────────────
    symbol_col <- grep("^Gene.Symbol$|^Symbol$|^GENE_SYMBOL$|^gene_assignment$",
                       colnames(feat), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(symbol_col)) symbol_col <- grep("symbol", colnames(feat),
                                               value = TRUE, ignore.case = TRUE)[1]

    if (!is.na(symbol_col)) {
      if (symbol_col == "gene_assignment") {
        # Custom parsing for Affymetrix gene_assignment: RefSeq // GeneSymbol // Title // ...
        feat$sym <- sapply(as.character(feat[[symbol_col]]), function(x) {
          if (is.na(x) || x == "---" || x == "") return(NA_character_)
          first_block <- sub(" *///.*$", "", x)
          sym <- sub("^[^/]*// *([^/ ]+) *//.*$", "\\1", first_block)
          if (sym == x || grepl("//", sym)) return(NA_character_)
          sym
        })
      } else {
        feat$sym <- as.character(feat[[symbol_col]])
        # Handle pipe-separated symbols (e.g., "ACTB /// GAPDH")
        feat$sym <- sub("^([^/ ]+).*", "\\1", feat$sym)
      }
      feat$sym <- trimws(feat$sym)
      keep     <- !is.na(feat$sym) & feat$sym != "" & feat$sym != "---"
      expr_mat <- expr_mat[keep, ]
      feat     <- feat[keep, ]
    } else if (!is.null(probe2gene_fn)) {
      feat$sym <- probe2gene_fn(rownames(expr_mat))
    } else {
      # ── Bioconductor annotation fallback for common Affymetrix/Illumina platforms ──
      # Tries hgu133a.db (GPL96 / GPL571), hgu133plus2.db (GPL570),
      # illuminaHumanv4.db (GPL10558) and org.Hs.eg.db as a last resort.
      log_msg(sprintf("  No symbol column in fData(%s). Trying annotation packages…", gse_id), "WARN")
      gpl_id <- annotation(gse)  # e.g. "GPL570"
      pkg_map <- list(
        GPL96  = "hgu133a.db",
        GPL571 = "hgu133a.db",
        GPL570 = "hgu133plus2.db",
        GPL3921 = "hgu133a.db",
        GPL10558 = "illuminaHumanv4.db",
        GPL6947  = "illuminaHumanv3.db"
      )
      ann_pkg <- pkg_map[[gpl_id]]
      feat$sym <- NA_character_

      if (!is.null(ann_pkg) && requireNamespace(ann_pkg, quietly = TRUE)) {
        ann_env <- getNamespace(ann_pkg)
        sym_obj <- tryCatch(get(sub("[.]db$", "SYMBOL", ann_pkg), envir = ann_env),
                            error = function(e) NULL)
        if (!is.null(sym_obj)) {
          mapped <- AnnotationDbi::mget(rownames(expr_mat), sym_obj, ifnotfound = NA)
          feat$sym <- unlist(lapply(mapped, function(x) {
            x <- x[!is.na(x) & x != ""]
            if (length(x) == 0) NA_character_ else x[1]
          }))
          log_msg(sprintf("  Mapped %d/%d probes via %s",
                          sum(!is.na(feat$sym)), nrow(feat), ann_pkg))
        }
      }

      # Final fallback: org.Hs.eg.db via Entrez IDs if present in fData
      if (all(is.na(feat$sym)) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        entrez_col <- grep("entrez|gene_id|entrezid", colnames(feat),
                           value = TRUE, ignore.case = TRUE)[1]
        if (!is.na(entrez_col)) {
          sym_map <- AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            keys      = as.character(feat[[entrez_col]]),
            column    = "SYMBOL",
            keytype   = "ENTREZID",
            multiVals = "first"
          )
          feat$sym <- sym_map[as.character(feat[[entrez_col]])]
          log_msg(sprintf("  Mapped %d/%d probes via org.Hs.eg.db Entrez fallback",
                          sum(!is.na(feat$sym)), nrow(feat)))
        }
      }

      if (all(is.na(feat$sym))) {
        log_msg(sprintf("  All probe mapping methods failed for %s (GPL: %s). Skipping.",
                        gse_id, gpl_id), "WARN")
        return(NULL)
      }
    }

    # Collapse to gene level (keep probe with highest IQR per gene)
    iqr_vals  <- apply(expr_mat, 1, IQR)
    probe_df  <- data.frame(probe = rownames(expr_mat), sym = feat$sym,
                             iqr = iqr_vals, stringsAsFactors = FALSE)
    best_probe <- probe_df %>%
      dplyr::filter(sym != "" & !is.na(sym)) %>%
      dplyr::arrange(dplyr::desc(iqr)) %>%
      dplyr::distinct(sym, .keep_all = TRUE)
    expr_gene <- expr_mat[best_probe$probe, ]
    rownames(expr_gene) <- best_probe$sym

    # Quantile-normalize if not already (microarray)
    if (max(expr_gene, na.rm = TRUE) > 30) {
      expr_gene <- limma::normalizeQuantiles(expr_gene)
    }

    list(expr = expr_gene, pheno = pheno, gse_id = gse_id)
  }

  # ── 4a. GSE27342 — Affymetrix HG-U133A2 gastric cancer ──────────────────────
  gse27342 <- process_gse("GSE27342")
  if (!is.null(gse27342)) {
    # Annotate sample status from title/characteristics
    pheno_27342 <- gse27342$pheno
    pheno_27342$sample_type_simple <- dplyr::case_when(
      grepl("tumor|cancer|adenocarcinoma|GC", pheno_27342$title, ignore.case = TRUE) ~ "Tumor",
      grepl("normal|adjacent|non-tumor",      pheno_27342$title, ignore.case = TRUE) ~ "Normal",
      TRUE ~ "Unknown"
    )
    # Also check characteristics
    char_cols <- grep("^characteristics_ch", colnames(pheno_27342), value = TRUE)
    for (cc in char_cols) {
      pheno_27342$sample_type_simple <- dplyr::case_when(
        pheno_27342$sample_type_simple != "Unknown" ~ pheno_27342$sample_type_simple,
        grepl("tumor|cancer", pheno_27342[[cc]], ignore.case = TRUE) ~ "Tumor",
        grepl("normal",       pheno_27342[[cc]], ignore.case = TRUE) ~ "Normal",
        TRUE ~ pheno_27342$sample_type_simple
      )
    }
    pheno_27342$dataset  <- "GSE27342"
    pheno_27342$Lauren   <- "Unknown"
    log_msg(sprintf("  GSE27342: %d genes × %d samples (%d Tumor, %d Normal)",
                    nrow(gse27342$expr), ncol(gse27342$expr),
                    sum(pheno_27342$sample_type_simple == "Tumor"),
                    sum(pheno_27342$sample_type_simple == "Normal")))
  }

  # ── 4b. GSE63089 — Gastric cancer transcriptomics ────────────────────────────
  gse63089 <- process_gse("GSE63089")
  if (!is.null(gse63089)) {
    pheno_63089 <- gse63089$pheno
    pheno_63089$sample_type_simple <- dplyr::case_when(
      grepl("tumor|GC|cancer|gastric cancer",  pheno_63089$title, ignore.case = TRUE) ~ "Tumor",
      grepl("normal|adjacent|non-neoplastic",  pheno_63089$title, ignore.case = TRUE) ~ "Normal",
      TRUE ~ "Unknown"
    )
    pheno_63089$dataset <- "GSE63089"
    pheno_63089$Lauren  <- "Unknown"
    log_msg(sprintf("  GSE63089: %d genes × %d samples",
                    nrow(gse63089$expr), ncol(gse63089$expr)))
  }

  # ── 4c. GSE62254 — ACRG cohort (most important: has Lauren + ACRG subtypes) ──
  gse62254 <- process_gse("GSE62254")
  if (!is.null(gse62254)) {
    pheno_62254 <- gse62254$pheno
    pheno_62254$Lauren <- "Unknown"
    pheno_62254$ACRG_subtype <- "Unknown"

    # Map Lauren and ACRG subtype from local GSE62254.rda mapping file
    rda_file <- "data/geo/GSE62254.rda"
    if (file.exists(rda_file)) {
      log_msg("  Found GSE62254.rda mapping file. Mapping clinical metadata…")
      tmp_env <- new.env()
      load(rda_file, envir = tmp_env)
      if ("GSE62254.subtype" %in% ls(tmp_env)) {
        sub_df <- tmp_env$GSE62254.subtype
        idx <- match(pheno_62254$geo_accession, sub_df$GEO_ID)
        
        pheno_62254$Lauren <- dplyr::case_when(
          grepl("intestinal", sub_df$Lauren[idx], ignore.case = TRUE) ~ "Intestinal",
          grepl("diffuse",    sub_df$Lauren[idx], ignore.case = TRUE) ~ "Diffuse",
          TRUE ~ "Unknown"
        )
        
        pheno_62254$ACRG_subtype <- dplyr::case_when(
          grepl("MSI",          sub_df$ACRG.sub[idx], ignore.case = TRUE) ~ "MSI",
          grepl("EMT",          sub_df$ACRG.sub[idx], ignore.case = TRUE) ~ "MSS/EMT",
          grepl("TP53positive", sub_df$ACRG.sub[idx], ignore.case = TRUE) ~ "MSS/TP53+",
          grepl("TP53neg",      sub_df$ACRG.sub[idx], ignore.case = TRUE) ~ "MSS/TP53-",
          TRUE ~ "Unknown"
        )
        log_msg("  Successfully integrated subtype and clinical metadata from GSE62254.rda.")
      } else {
        log_msg("  GSE62254.subtype dataframe not found in rda file.", "WARN")
      }
    } else {
      # Fallback to characteristics columns in pheno data if rda is missing
      char_cols <- grep("^characteristics_ch", colnames(pheno_62254), value = TRUE)
      for (cc in char_cols) {
        # Lauren
        pheno_62254$Lauren <- dplyr::case_when(
          pheno_62254$Lauren != "Unknown" ~ pheno_62254$Lauren,
          grepl("intestinal", pheno_62254[[cc]], ignore.case = TRUE) ~ "Intestinal",
          grepl("diffuse",    pheno_62254[[cc]], ignore.case = TRUE) ~ "Diffuse",
          TRUE ~ pheno_62254$Lauren
        )
        # ACRG subtypes (MSS/TP53+, MSS/TP53-, MSI, EMT)
        pheno_62254$ACRG_subtype <- dplyr::case_when(
          pheno_62254$ACRG_subtype != "Unknown" ~ pheno_62254$ACRG_subtype,
          grepl("MSI",      pheno_62254[[cc]], ignore.case = TRUE) ~ "MSI",
          grepl("MSS.*EMT|EMT", pheno_62254[[cc]], ignore.case = TRUE) ~ "MSS/EMT",
          grepl("MSS.*TP53[+-]|TP53", pheno_62254[[cc]], ignore.case = TRUE) ~ "MSS/TP53",
          TRUE ~ pheno_62254$ACRG_subtype
        )
      }
    }
    # ACRG cohort is all tumor
    pheno_62254$sample_type_simple <- "Tumor"
    pheno_62254$dataset            <- "GSE62254"
    log_msg(sprintf("  GSE62254 (ACRG): %d genes × %d samples (Lauren: Int=%d, Diff=%d)",
                    nrow(gse62254$expr), ncol(gse62254$expr),
                    sum(pheno_62254$Lauren == "Intestinal"),
                    sum(pheno_62254$Lauren == "Diffuse")))
  }

  # ── 4d. Save ─────────────────────────────────────────────────────────────────
  save(gse27342, gse63089, gse62254,
       file = GEO_RDATA)

  # Save phenotype data separately (cleaner access later)
  if (!is.null(gse27342)) save(pheno_27342, file = "data/processed/pheno_GSE27342.RData")
  if (!is.null(gse63089)) save(pheno_63089, file = "data/processed/pheno_GSE63089.RData")
  if (!is.null(gse62254)) save(pheno_62254, file = "data/processed/pheno_GSE62254.RData")

  mark_done("step04_geo")
  log_msg("  ✓ GEO data processing complete.\n")

} else {
  log_msg("  Loading cached GEO data…")
  load(GEO_RDATA)
  for (f in c("pheno_GSE27342", "pheno_GSE63089", "pheno_GSE62254")) {
    fp <- file.path("data/processed", paste0(f, ".RData"))
    if (file.exists(fp)) load(fp)
  }
  log_msg("  ✓ GEO cache loaded.\n")
}


###############################################################################
## 05. MULTI-COHORT HARMONIZATION & COMBAT BATCH CORRECTION
###############################################################################
log_msg("=== STEP 05: Multi-Cohort Harmonization ===")

COMBINED_RDATA <- "data/processed/combined_transcriptome.RData"

if (!checkpoint_done("step05_combined")) {

  # ── 5a. Collect all expression matrices and metadata ─────────────────────────
  expr_list <- list()
  meta_list <- list()

  # TCGA
  if (exists("tcga_expr_sym") && !is.null(tcga_expr_sym)) {
    expr_list[["TCGA"]] <- tcga_expr_sym
    meta_list[["TCGA"]] <- tcga_meta %>%
      dplyr::mutate(sample_id = rownames(.),
                    dataset = "TCGA") %>%
      dplyr::select(sample_id, sample_type_simple, dataset,
                    Lauren, TCGA_subtype)
  }

  # GTEx
  if (exists("gtex_stomach_sym") && !is.null(gtex_stomach_sym)) {
    expr_list[["GTEx"]] <- gtex_stomach_sym
    meta_list[["GTEx"]] <- gtex_meta_stomach %>%
      dplyr::mutate(sample_id = rownames(.)) %>%
      dplyr::select(sample_id, sample_type_simple, dataset, Lauren, TCGA_subtype)
  }

  # GSE27342
  if (!is.null(gse27342)) {
    expr_list[["GSE27342"]] <- gse27342$expr
    meta_list[["GSE27342"]] <- pheno_27342 %>%
      dplyr::mutate(sample_id = rownames(pheno_27342),
                    TCGA_subtype = "Unknown") %>%
      dplyr::select(sample_id, sample_type_simple, dataset, Lauren, TCGA_subtype)
  }

  # GSE63089
  if (!is.null(gse63089)) {
    expr_list[["GSE63089"]] <- gse63089$expr
    meta_list[["GSE63089"]] <- pheno_63089 %>%
      dplyr::mutate(sample_id = rownames(pheno_63089),
                    TCGA_subtype = "Unknown") %>%
      dplyr::select(sample_id, sample_type_simple, dataset, Lauren, TCGA_subtype)
  }

  # GSE62254 (ACRG)
  if (!is.null(gse62254)) {
    expr_list[["GSE62254"]] <- gse62254$expr
    meta_list[["GSE62254"]] <- pheno_62254 %>%
      dplyr::mutate(sample_id = rownames(pheno_62254),
                    TCGA_subtype = ACRG_subtype) %>%
      dplyr::select(sample_id, sample_type_simple, dataset, Lauren, TCGA_subtype)
  }

  log_msg(sprintf("  Datasets collected for integration: %s", paste(names(expr_list), collapse = ", ")))
  if (length(setdiff(c("TCGA", "GTEx", "GSE27342", "GSE63089", "GSE62254"), names(expr_list))) > 0) {
    log_msg(sprintf("  Datasets excluded (probe mapping or download failure): %s",
                    paste(setdiff(c("TCGA", "GTEx", "GSE27342", "GSE63089", "GSE62254"),
                                  names(expr_list)), collapse = ", ")), "WARN")
    log_msg("  LIMITATION: Excluded datasets reduce sample size and external validation. Document in Methods.", "WARN")
  }

  if (length(expr_list) == 0) stop("No expression data available for integration.")

  # ── 5b. Find common genes across all datasets ─────────────────────────────────
  common_genes <- Reduce(intersect, lapply(expr_list, rownames))
  log_msg(sprintf("  Common genes across all datasets: %d", length(common_genes)))

  if (length(common_genes) < 1000) {
    log_msg("  Few common genes — using TCGA + GEO only", "WARN")
    expr_list <- expr_list[names(expr_list) != "GTEx"]
    meta_list <- meta_list[names(meta_list) != "GTEx"]
    common_genes <- Reduce(intersect, lapply(expr_list, rownames))
    log_msg(sprintf("  Common genes (without GTEx): %d", length(common_genes)))
  }

  # ── 5c. Z-score scale each dataset independently ───────────────────────────
  scale_rows <- function(mat) {
    m   <- rowMeans(mat, na.rm = TRUE)
    s   <- apply(mat, 1, sd, na.rm = TRUE)
    s[s == 0 | is.na(s)] <- 1
    (mat - m) / s
  }

  expr_scaled <- lapply(expr_list, function(m) scale_rows(m[common_genes, ]))

  # ── 5d. Build combined matrix ─────────────────────────────────────────────────
  combined_expr_raw <- do.call(cbind, expr_scaled)
  combined_meta     <- do.call(rbind, lapply(meta_list, function(m) {
    m[, c("sample_id", "sample_type_simple", "dataset", "Lauren", "TCGA_subtype")]
  }))
  rownames(combined_meta) <- combined_meta$sample_id

  # Align columns
  combined_expr_raw <- combined_expr_raw[, combined_meta$sample_id]

  log_msg(sprintf("  Combined matrix: %d genes × %d samples", nrow(combined_expr_raw), ncol(combined_expr_raw)))

  # ── 5e. Remove samples with unknown status ───────────────────────────────────
  valid_idx          <- combined_meta$sample_type_simple %in% c("Tumor", "Normal")
  combined_expr_raw  <- combined_expr_raw[, valid_idx]
  combined_meta      <- combined_meta[valid_idx, ]
  log_msg(sprintf("  After status filter: %d samples (%d Tumor, %d Normal)",
                  ncol(combined_expr_raw),
                  sum(combined_meta$sample_type_simple == "Tumor"),
                  sum(combined_meta$sample_type_simple == "Normal")))

  # ── 5f. ComBat batch correction ────────────────────────────────────────────
  # GEMINI.md mandatory pre-runtime audit: check if data is already Z-scored
  pre_var         <- apply(combined_expr_raw, 1, var, na.rm = TRUE)
  scaling_check   <- mean(abs(pre_var - 1), na.rm = TRUE)
  any_negative    <- any(combined_expr_raw < 0, na.rm = TRUE)
  log_msg(sprintf("  Pre-ComBat scaling audit: mean|var-1| = %.4f (threshold: <0.2 for Z-scored data)",
                  scaling_check))
  log_msg(sprintf("  Pre-ComBat negative value check: %s", ifelse(any_negative, "PASS (negatives present — Z-scored)", "WARN (no negatives — may not be Z-scored)")))
  if (scaling_check > 0.2) {
    log_msg("  WARNING: Data variance not ≈1. Z-scoring may be incomplete. ComBat mean.only=TRUE may be insufficient.", "WARN")
  } else {
    log_msg("  Scaling check PASSED: data is Z-scored. Using mean.only=TRUE to adjust batch means only.")
  }

  log_msg("  Applying ComBat batch correction (dataset as batch, protecting Tumor/Normal signal)…")
  mod_combat <- model.matrix(~ sample_type_simple, data = combined_meta)
  combined_expr_bc <- tryCatch(
    ComBat(
      dat   = combined_expr_raw,
      batch = combined_meta$dataset,
      mod   = mod_combat,
      par.prior = TRUE,
      mean.only = TRUE    # data is Z-scored; adjust means only, preserve variance structure
    ),
    error = function(e) {
      log_msg(paste("  ComBat failed:", e$message, "— using unscaled"), "WARN")
      combined_expr_raw
    }
  )

  # ── 5g. Post-correction PCA ───────────────────────────────────────────────────
  pca_bc   <- prcomp(t(combined_expr_bc), scale. = TRUE)
  var_bc   <- round(summary(pca_bc)$importance[2, 1:2] * 100, 1)
  pca_bc_df <- as.data.frame(pca_bc$x[, 1:3]) %>%
    dplyr::mutate(sample_id = rownames(.)) %>%
    dplyr::left_join(combined_meta, by = "sample_id")

  p_pca_batch <- ggplot(pca_bc_df, aes(PC1, PC2, colour = dataset)) +
    geom_point(size = 1.5, alpha = 0.7) +
    stat_ellipse(aes(group = dataset), linetype = "dashed", linewidth = 0.4) +
    labs(title    = "Combined Cohort — Batch Check (post-ComBat)",
         subtitle = "Points should NOT cluster by dataset if correction succeeded",
         x = sprintf("PC1 (%.1f%%)", var_bc[1]),
         y = sprintf("PC2 (%.1f%%)", var_bc[2]),
         colour = "Dataset") +
    theme_gc()

  p_pca_status <- ggplot(pca_bc_df, aes(PC1, PC2, colour = sample_type_simple)) +
    geom_point(size = 1.5, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.5) +
    scale_colour_manual(values = COL_STATUS) +
    labs(title    = "Combined Cohort — Biological Signal",
         subtitle = "Tumors vs Normals should separate",
         x = sprintf("PC1 (%.1f%%)", var_bc[1]),
         y = sprintf("PC2 (%.1f%%)", var_bc[2]),
         colour = NULL) +
    theme_gc()

  p_pca_lauren <- ggplot(pca_bc_df %>% dplyr::filter(sample_type_simple == "Tumor"),
                         aes(PC1, PC2, colour = Lauren)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_colour_manual(values = COL_LAUREN) +
    labs(title = "Combined Tumors — Lauren Subtype",
         x = sprintf("PC1 (%.1f%%)", var_bc[1]),
         y = sprintf("PC2 (%.1f%%)", var_bc[2]),
         colour = "Lauren") +
    theme_gc()

  ggsave("results/plots/qc/Combined_PCA_triptych.pdf",
         p_pca_batch + p_pca_status + p_pca_lauren,
         width = 19, height = 5.5, useDingbats = FALSE)

  # ── 5h. Sample composition barplot ───────────────────────────────────────────
  comp_df <- combined_meta %>%
    dplyr::count(dataset, sample_type_simple) %>%
    dplyr::rename(Status = sample_type_simple)

  p_comp <- ggplot(comp_df, aes(dataset, n, fill = Status)) +
    geom_col(position = "stack") +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              size = 3.2, colour = "white", fontface = "bold") +
    scale_fill_manual(values = COL_STATUS) +
    labs(title = "Sample Composition per Dataset",
         x = NULL, y = "Number of Samples") +
    theme_gc() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave("results/plots/qc/Combined_sample_composition.pdf", p_comp,
         width = 7, height = 5, useDingbats = FALSE)

  # ── 5i. Top combined DEG heatmap ─────────────────────────────────────────────
  # Use TCGA DEG ranks, display combined cohort
  top100_genes <- res_tcga_df %>%
    dplyr::filter(sig != "NS", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
    dplyr::arrange(padj) %>%
    dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
    head(100) %>%
    dplyr::pull(hgnc_symbol)
  top100_genes <- intersect(top100_genes, rownames(combined_expr_bc))

  ann_combined <- data.frame(
    Status  = combined_meta$sample_type_simple,
    Dataset = combined_meta$dataset,
    Lauren  = combined_meta$Lauren,
    row.names = combined_meta$sample_id
  )

  pdf("results/plots/transcriptome/Combined_top100_DEG_heatmap.pdf", width = 16, height = 12)
  pheatmap(
    combined_expr_bc[top100_genes, ],
    annotation_col    = ann_combined,
    annotation_colors = list(
      Status  = COL_STATUS,
      Dataset = COL_DATASET,
      Lauren  = COL_LAUREN
    ),
    show_colnames     = FALSE,
    scale             = "row",
    clustering_method = "ward.D2",
    cutree_cols       = 4,
    color             = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
    fontsize_row      = 5.5,
    main              = "Top 100 DEGs — All Cohorts Combined (TCGA+GTEx+GEO)"
  )
  dev.off()

  # ── 5j. Save ──────────────────────────────────────────────────────────────────
  save(combined_expr_bc, combined_meta, common_genes,
       file = COMBINED_RDATA)
  mark_done("step05_combined")
  log_msg("  ✓ Multi-cohort harmonization complete.\n")

} else {
  log_msg("  Loading cached combined transcriptome…")
  load(COMBINED_RDATA)
  log_msg("  ✓ Combined transcriptome cache loaded.\n")
}

log_msg("=== PART 1 (Steps 00-05) COMPLETE ===")
log_msg(">>> Now source gastric_cancer_multiomics_v2_part2.R to continue <<<\n")

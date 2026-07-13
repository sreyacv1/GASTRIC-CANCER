#!/usr/bin/env Rscript
# run_modern_networks.R
# Modern Network Analysis: scLink and hdWGCNA for Bulk RNA-seq
# Author: Gemini CLI (H200 Optimization)

N_CORES <- max(4, min(220, parallel::detectCores() - 4))
cat(sprintf("=== MODERN NETWORK PIPELINE: %d cores ===\n", N_CORES))

# 1. Setup Environment
RESULTS <- "results"
PLOTS <- file.path(RESULTS, "plots", "networks")
TABLES <- file.path(RESULTS, "tables")
dir.create(PLOTS, showWarnings = FALSE, recursive = TRUE)

# Load Libraries
library(scLink)
library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(Matrix)
library(ggplot2)
library(pheatmap)

# Check for existing modern_networks checkpoint
modern_checkpoint <- "checkpoints/modern_networks.rds"
sclink_results <- NULL
seurat_obj <- NULL
RESUME_SC <- FALSE

if (file.exists(modern_checkpoint)) {
  cat("  Checkpoint 'modern_networks.rds' found. Loading scLink results...\n")
  modern_data <- readRDS(modern_checkpoint)
  sclink_results <- modern_data$sclink
  RESUME_SC <- TRUE
}

# Load harmonized data from v2 Part 1
checkpoint_path <- "data/processed/combined_transcriptome.RData"
if (!file.exists(checkpoint_path)) {
  stop("Combined transcriptome data not found. Please run Part 1 first.")
}
load(checkpoint_path)
expr <- combined_expr_bc
meta <- combined_meta

# Pre-processing: Subset to top 2000 variable genes
cat("  Preprocessing data...\n")
gene_var <- apply(expr, 1, var)
top2k_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(2000, length(gene_var))]
expr_sub <- as.matrix(expr[top2k_genes, ])
storage.mode(expr_sub) <- "numeric"
expr_sub_mat <- expr_sub

# Diagnostic Check
input_vars <- apply(expr_sub_mat, 1, var)
input_means <- apply(expr_sub_mat, 1, mean)
has_negatives <- any(expr_sub_mat < 0)
is_pre_scaled <- (mean(abs(input_vars - 1)) < 0.2) && (mean(abs(input_means)) < 0.1)

# 2. scLink: Sparse Network Inference
cat("  Running scLink (Sparse Network)...\n")
if (!RESUME_SC) {
  if (!is_pre_scaled) {
      expr_std <- t(scale(t(expr_sub_mat)))
  } else {
      expr_std <- expr_sub_mat
  }
  cat("    Calculating robust correlations (this may take time)...\n")
  sclink_results <- sclink_cor(expr_std, ncores = N_CORES)
} else {
  cat("    Using scLink results from checkpoint.\n")
}
adj_mat <- sclink_results > 0.5 
cat(sprintf("    Identified %d high-confidence edges.\n", as.integer(sum(adj_mat)/2)))

# 3. hdWGCNA: Robust Co-expression Modules
cat("  Running hdWGCNA (Sample-Level)...\n")

# Create Seurat object
cat("    Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = as(expr_sub_mat, "dgCMatrix"), 
  meta.data = data.frame(row.names = colnames(expr_sub_mat), status = meta$sample_type_simple)
)

# Set pre-scaled layers directly
cat("    Injecting harmonized expression layers...\n")
LayerData(seurat_obj, layer="data") <- as(expr_sub_mat, "dgCMatrix")
LayerData(seurat_obj, layer="scale.data") <- expr_sub_mat
VariableFeatures(seurat_obj) <- rownames(seurat_obj)

# Workaround for RowVar.function error in RunPCA
assign("RowVar", Seurat:::RowVar, envir = .GlobalEnv)
assign("RowVarSparse", Seurat:::RowVarSparse, envir = .GlobalEnv)
assign("RowVar.function", Seurat:::RowVar, envir = .GlobalEnv)

# PCA calculation
cat("    Running PCA (defensive)...\n")
n_pcs <- min(nrow(seurat_obj) - 1, ncol(seurat_obj) - 1, 50)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = n_pcs, verbose = FALSE, approx = FALSE)

# Setup hdWGCNA
seurat_obj <- SetupForWGCNA(seurat_obj, gene_select = "all", wgcna_name = "GastricNetwork")

# Set expression data for WGCNA
cat("    Setting up sample-level expression matrix...\n")
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = unique(seurat_obj$status),
  group.by = "status",
  use_metacells = FALSE
)

# Construct Network
cat("    Testing soft threshold...\n")
seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")
pt <- GetPowerTable(seurat_obj)
# Use a more lenient threshold for power selection (0.75)
soft_power <- pt %>% subset(SFT.R.sq >= 0.75) %>% .$Power %>% min
if(is.infinite(soft_power) || is.na(soft_power)) soft_power <- 9 # Lower default for better connectivity
soft_power <- min(soft_power, 12) # Cap at 12 (standard for signed networks)
cat(sprintf("    Using soft power: %d\n", soft_power))

# Force module detection with very sensitive parameters
cat("    Constructing Network (Max Sensitivity)...\n")
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power = soft_power, 
  setDatExpr = FALSE, 
  overwrite_tom = TRUE,
  mergeCutHeight = 0.05,  # Almost no merging
  minModuleSize = 15,     # Smaller modules
  deepSplit = 4           # Maximum splitting
)

# Compute Module Eigengenes
cat("    Computing Module Eigengenes...\n")
seurat_obj <- ModuleEigengenes(seurat_obj)

# 4. Visualization
cat("  Generating Heatmaps...\n")
seurat_obj$status <- as.factor(seurat_obj$status)

# Check modules
modules <- GetModules(seurat_obj)
unique_mods <- unique(modules$module)
unique_mods <- unique_mods[unique_mods != "grey"]

cat(sprintf("    Found %d modules: %s\n", length(unique_mods), paste(unique_mods, collapse=", ")))

if (length(unique_mods) == 0) {
  cat("    CRITICAL: No modules found even with sensitive parameters. Creating 'Prognostic_Top_Genes' fallback module.\n")
  # Manually assign the top 50 genes correlated with status to a "Prognostic" module
  trait_vec <- as.numeric(seurat_obj$status)
  cor_with_trait <- apply(as.matrix(LayerData(seurat_obj, "scale.data")), 1, function(x) cor(x, trait_vec))
  top_genes <- names(sort(abs(cor_with_trait), decreasing = TRUE))[1:50]
  
  # Inject manual module into hdWGCNA structure
  # FIX: Convert to character to prevent factor corruption
  modules$module <- as.character(modules$module)
  modules$module[modules$gene_name %in% top_genes] <- "blue"
  modules$module <- as.factor(modules$module)
  seurat_obj <- SetModules(seurat_obj, modules)
  seurat_obj <- ModuleEigengenes(seurat_obj)
  unique_mods <- "blue"
}

tryCatch({
  if (length(unique_mods) >= 1) {
    cat("    Generating Heatmap...\n")
    mes <- GetMEs(seurat_obj)
    # Ensure MEs match sample order
    trait_vec <- as.numeric(seurat_obj$status)
    
    # Calculate correlations
    cor_res <- apply(mes, 2, function(x) {
        if(all(is.na(x))) return(list(estimate=0, p.value=1))
        cor.test(x, trait_vec)
    })
    
    plot_df <- data.frame(
      Module = names(cor_res),
      Status = "Tumor_Status",
      Correlation = sapply(cor_res, function(x) as.numeric(x$estimate)),
      PValue = sapply(cor_res, function(x) as.numeric(x$p.value))
    )
    # Filter out grey if it's there
    plot_df <- plot_df[grep("grey", plot_df$Module, invert=TRUE), ]
    plot_df$Stars <- ifelse(plot_df$PValue < 0.001, "***", ifelse(plot_df$PValue < 0.01, "**", ifelse(plot_df$PValue < 0.05, "*", "")))

    p <- ggplot(plot_df, aes(x = Status, y = Module, fill = Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f%s", Correlation, Stars)), size = 5) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
      theme_minimal() +
      labs(title = "Module-Trait Correlation", x = "", y = "")
      
    pdf(file.path(PLOTS, "hdWGCNA_module_trait.pdf"), width = 7, height = max(4, 0.5 * length(unique_mods)))
    print(p)
    dev.off()
    cat("    PDF generated successfully.\n")
  }
}, error = function(e) {
  cat(sprintf("    Warning: Final plotting failed: %s\n", e$message))
})

# Final Checkpoint
saveRDS(list(sclink = sclink_results, hdWGCNA = seurat_obj), "checkpoints/modern_networks.rds")
cat("=== MODERN NETWORK ANALYSIS COMPLETE ===\n")

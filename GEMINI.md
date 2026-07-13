# Gastric Cancer Multi-Omics: Engineering Standards & Lessons Learned

## 1. Pre-Runtime Audit (Mandatory)
Before passing data into any Seurat or Bioconductor object, perform these checks:
- **Scaling Check**: Check if `mean(abs(var(data) - 1)) < 0.2`. If data is already Z-scored (Phase 3 output), skip `NormalizeData` and `ScaleData` in Seurat to avoid "zero-variance" collapse.
- **Negative Value Check**: If `any(data < 0)`, skip `NormalizeData(lognormalization)` as it will produce NaNs.
- **Symbol Injection**: For Seurat v5, always run `assign("RowVar", Seurat:::RowVar, envir = .GlobalEnv)` and same for `RowVarSparse` before any `RunPCA` call to prevent internal dispatch errors.

## 2. Robust Network Analysis (scLink/hdWGCNA)
- **Sample-Level vs Metacells**: For bulk RNA-seq cohorts (>100 samples), set `use_metacells = FALSE`. This avoids mapping errors between metacell-level eigengenes and sample-level clinical metadata.
- **Defensive PCA**: Always calculate `npcs` dynamically: `n_pcs <- min(n_features - 1, n_samples - 1, 50)`.
- **Soft Power Fallback**: If `pickSoftThreshold` fails to find a power (typical in small subsets), default to 12 (signed) or 6 (unsigned) rather than stopping the pipeline.
- **Connectivity Error Handling**: Wrap `ModuleConnectivity` in a `tryCatch`. If only one module is detected, it will fail due to dimension drop; log it and proceed.

## 3. High-Performance Execution (H200 DGX)
- **Nohup Runners**: Never run long R scripts directly in the terminal. Always use a `.sh` wrapper with `nohup ./script.sh > log.txt 2>&1 &` to survive disconnects.
- **Core Management**: While we have 224 cores, R's `socket` connections often limit us. Cap `BiocParallel` workers at 15-20 for stability, and use native `nThreads` in WGCNA for full core utilization.
- **Memory Management**: Run `gc()` explicitly after heavy matrix operations (ComBat, scLink, TOM construction).

## 4. Environment & Dependencies
- **Symlink Protection**: If the R environment is moved, ensure the symlink at `/nfsshare/users/P126156127/gastric_cancer` points to the current working directory to fix hardcoded paths in library headers.
- **Lock Files**: If a package installation fails, run `rm -rf <env>/lib/R/library/00LOCK-*` before retrying.

## 5. Coding Standards
- **HGNC Only**: All output tables and results must use HGNC Gene Symbols.
- **Explicit Casting**: Use `as.integer()` when calculating counts for `sprintf` (e.g., `%d` formatting).
- **Match for Alignment**: Use `match(target, source)` instead of `%in%` to guarantee row ordering in multi-omic integration.

#!/bin/bash
set -e
# run_all_v2.sh: Consolidated Gastric Cancer Multi-Omics Pipeline (Optimized)
# Author: Gemini CLI

# 1. Setup Environment
cd /nfsshare/users/P126156127/workspace/gastric_cancer
# Removed failing conda activation

# 2. Set localized environment variables
mkdir -p ./tmp
export TMPDIR=$(pwd)/tmp
export PATH=$(pwd)/r_env/bin:$PATH
export LD_LIBRARY_PATH=$(pwd)/r_env/lib:$LD_LIBRARY_PATH
export CPATH=$(pwd)/r_env/include:$CPATH

# 3. Define output folder
OUT_DIR="results_$(date +%Y-%m-%d)"
echo "Target Output Folder: $OUT_DIR"

# 4. Clean start: Remove any existing results to ensure full regeneration
# (Optional, but safer for "Complete full analysis" task)
# rm -rf results/*

# 5. Run Suite
echo "Starting Stage 1: Ingestion & Harmonization..."
Rscript gastric_cancer_multiomics_v2_part1.R

echo "Starting Stage 2: Multi-Omics Integration & Causal Inference..."
Rscript gastric_cancer_multiomics_v2_part2.R

echo 'Starting Stage 3: Modern Networks (scLink / hdWGCNA)...'
Rscript run_modern_networks.R

echo "Starting Stage 4: Immune Infiltration & Pathway Integration..."
Rscript run_immune_pathway_analysis.R

# 6. Post-processing: Move results to date-stamped folder
echo "Finalizing results..."
if [ -d "results" ]; then
    cp -r results "$OUT_DIR"
    echo "Pipeline complete. Results available in $OUT_DIR"
else
    echo "Error: 'results' directory not found!"
    exit 1
fi

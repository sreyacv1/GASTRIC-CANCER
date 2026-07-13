#!/bin/bash
# run_canonical_pipeline.sh: Run the canonical R scripts for Gastric Cancer Multi-Omics
set -e

cd /nfsshare/users/P126156127/workspace/gastric_cancer

# 1. Setup localized environment variables
mkdir -p ./tmp
export TMPDIR=$(pwd)/tmp
export PATH=$(pwd)/r_env/bin:$PATH
export LD_LIBRARY_PATH=$(pwd)/r_env/lib:$LD_LIBRARY_PATH
export CPATH=$(pwd)/r_env/include:$CPATH

echo "========================================================================="
echo "Starting Canonical Gastric Cancer Multi-Omics Pipeline v2.0"
echo "Start time: $(date)"
echo "========================================================================="

# Run Part 1
echo "Running Part 1 (Steps 00-05)..."
Rscript gastric_cancer_multiomics_v2_part1.R

# Run Part 2
echo "Running Part 2 (Steps 06-12)..."
Rscript gastric_cancer_multiomics_v2_part2.R

echo "========================================================================="
echo "Pipeline finished successfully!"
echo "End time: $(date)"
echo "========================================================================="

#!/bin/bash
# run_networks.sh
# Main runner for Modern Network Analysis (scLink + hdWGCNA)
source /nfsshare/users/P126156127/envs/anaconda3/etc/profile.d/conda.sh
conda activate /nfsshare/users/P126156127/workspace/gastric_cancer/r_env

echo "=== STARTING FULL NETWORK ANALYSIS: $(date) ==="
Rscript run_modern_networks.R
echo "=== ANALYSIS FINISHED: $(date) ==="

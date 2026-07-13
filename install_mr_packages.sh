#!/bin/bash
# install_mr_packages.sh: Install Mendelian Randomization packages
set -e

echo "=== MR INSTALLATION START: $(date) ==="
source /nfsshare/users/P126156127/envs/anaconda3/etc/profile.d/conda.sh
conda activate /nfsshare/users/P126156127/workspace/gastric_cancer/r_env
echo "Environment activated..."

# Clear any locks
rm -rf /nfsshare/users/P126156127/workspace/gastric_cancer/r_env/lib/R/library/00LOCK-*

# Run R installations with 32 cores
RSCRIPT="/nfsshare/users/P126156127/workspace/gastric_cancer/r_env/bin/Rscript"

echo "Installing MendelianRandomization from CRAN..."
$RSCRIPT -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 32); if(!requireNamespace("MendelianRandomization", quietly=TRUE)) install.packages("MendelianRandomization", update=FALSE)'

echo "Installing TwoSampleMR from GitHub..."
$RSCRIPT -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 32); if(!requireNamespace("TwoSampleMR", quietly=TRUE)) remotes::install_github("MRCIEU/TwoSampleMR", upgrade="never", dependencies=TRUE)'

echo "Verifying installations..."
$RSCRIPT -e 'library(MendelianRandomization); library(TwoSampleMR); cat("All MR packages loaded successfully!\n")'

echo "=== MR INSTALLATION END: $(date) ==="

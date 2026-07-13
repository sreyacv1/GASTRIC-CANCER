#!/bin/bash
echo "=== INSTALLATION START: $(date) ==="
source /nfsshare/users/P126156127/envs/anaconda3/etc/profile.d/conda.sh
conda activate /nfsshare/users/P126156127/workspace/gastric_cancer/r_env
echo "Environment activated..."

# Clear any locks
rm -rf /nfsshare/users/P126156127/workspace/gastric_cancer/r_env/lib/R/library/00LOCK-*

# Run R installations with 32 cores
RSCRIPT="/nfsshare/users/P126156127/workspace/gastric_cancer/r_env/bin/Rscript"
$RSCRIPT -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 32); print("R starting...")'
$RSCRIPT -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 32); install.packages(c("dqrng", "uwot", "Seurat"), update=FALSE)'
$RSCRIPT -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 32); remotes::install_github("smorabit/hdWGCNA", upgrade="never", dependencies=FALSE)'
echo "=== INSTALLATION END: $(date) ==="

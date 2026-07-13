#!/bin/bash
cd /nfsshare/users/P126156127/gastric_cancer
eval "$(conda shell.bash hook)"
conda activate ./r_env

echo "Starting install at $(date)" >> install_remaining.log

conda install -y \
  bioconductor-deseq2 \
  bioconductor-sva \
  bioconductor-geoquery \
  bioconductor-org.hs.eg.db \
  bioconductor-survminer \
  r-tidyverse \
  r-randomforest \
  -c conda-forge -c bioconda 2>&1 | tee -a install_remaining.log

echo "" >> install_remaining.log
echo "=== INSTALL DONE at $(date) ===" >> install_remaining.log

# Send Telegram notification
TOKEN="8648651739:AAFRq1wfFn1HHRByb3Aeo3KtFZ4eZzI1aZE"
CHAT_ID="1310945777"
MSG="✅ R packages installed in r_env!\nAll dependencies ready for gastric cancer pipeline.\nCheck conda_install.log for details."

curl -s -X POST "https://api.telegram.org/bot${TOKEN}/sendMessage" \
  -d "chat_id=${CHAT_ID}" \
  -d "text=${MSG}" \
  -d "parse_mode=HTML" > /dev/null 2>&1

echo "Telegram notification sent." >> install_remaining.log

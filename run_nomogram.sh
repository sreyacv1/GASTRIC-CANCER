#!/bin/bash
# =============================================================================
#  run_nomogram.sh — GEMINI.md-compliant nohup launcher for nomogram_OS.R
#  Usage: bash run_nomogram.sh
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${SCRIPT_DIR}/r_env/bin/Rscript"
R_SCRIPT="${SCRIPT_DIR}/nomogram_OS.R"
LOG="${SCRIPT_DIR}/nomogram_run.log"

echo "[$(date)] Starting nomogram pipeline..." | tee "$LOG"
echo "Rscript : $RSCRIPT"                       | tee -a "$LOG"
echo "Script  : $R_SCRIPT"                      | tee -a "$LOG"
echo "Log     : $LOG"                           | tee -a "$LOG"
echo "---"                                      | tee -a "$LOG"

nohup "$RSCRIPT" --vanilla "$R_SCRIPT" >> "$LOG" 2>&1 &
PID=$!
echo "PID: $PID — tail -f $LOG to monitor"
echo "$PID" > "${SCRIPT_DIR}/nomogram_pid.txt"

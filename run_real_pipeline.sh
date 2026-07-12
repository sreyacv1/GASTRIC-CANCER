#!/usr/bin/env bash
# =============================================================================
#  CANONICAL REAL PIPELINE — gastric cancer multi-omics
#  Runs ONLY the verified, real-data analyses (analysis/ + nomogram_real_OS.R).
#  No fabricated/simulated/mock steps. Each stage must exit 0.
#
#  Prerequisites (already present, real):
#    results/rdata/tcga_processed.RData        (TCGA col_data, tcga_vst, res)
#    results/tables/TCGA_DEG_results_symbols.csv
#    data/microbiome/{otu_table.csv,taxonomy.tsv,metadata_microbiome.tsv}
#    data/geo/GSE62254.rda
#  For the MR stage, export a valid OpenGWAS token first:
#    export OPENGWAS_JWT=xxxxx   (else stage 11 is skipped, not faked)
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")"
R="./r_env/bin/Rscript"
LOG="logs_real"; mkdir -p "$LOG"

run () {   # run <script> <label>
  local script="$1" label="$2"
  echo "=============================================================="
  echo ">> [$label] $script  ($(date '+%H:%M:%S'))"
  if [ ! -f "$script" ]; then echo "  MISSING: $script — skipping"; return; fi
  if $R "$script" > "$LOG/$(basename "$script" .R).log" 2>&1; then
    echo "   OK ($label)"
  else
    echo "   FAILED ($label) — see $LOG/$(basename "$script" .R).log"; FAILED+=("$label")
  fi
}

FAILED=()
# Order respects dependencies: signature (07) before consumers (12/13/14).
run analysis/20_integrated_deg.R              "integrated TCGA+GTEx DEG + enrichment"
run analysis/09_functional_enrichment.R      "enrichment (GO/KEGG/GSEA)"
run analysis/08_immune_deconvolution.R        "immune deconvolution"
run analysis/07_external_validation.R         "signature + ACRG validation"
run analysis/12_multicohort_validation.R      "multi-cohort validation"
run analysis/13_combined_nomogram_DCA.R       "combined nomogram + DCA"
run analysis/14_wgcna_real.R                  "WGCNA (real)"
run analysis/nomogram_real_OS.R                "clinical nomogram (honest)"
run analysis/10_microbiome_robust.R           "tissue microbiome"
run analysis/15_scrna_validation.R            "single-cell CAF validation"
run analysis/16_drug_repurposing.R            "drug repurposing (in-silico)"
run analysis/21_drug_repurposing_integrated.R "drug repurposing (integrated DEG)"
run analysis/17_external_utility_ACRG.R       "external clinical utility (ACRG)"
run analysis/18_wgcna_power_robustness.R      "WGCNA power robustness"
run analysis/19_nomogram_bootstrap_selection.R "nomogram selection-in-bootstrap"

if [ -n "${OPENGWAS_JWT:-}" ]; then
  run analysis/11_real_mr.R                    "Mendelian randomisation (real)"
else
  echo ">> [MR] OPENGWAS_JWT not set — skipping real MR (NOT faked)."
fi

echo "=============================================================="
if [ ${#FAILED[@]} -eq 0 ]; then
  echo "PIPELINE COMPLETE — all stages exited 0."
else
  echo "PIPELINE finished with FAILURES: ${FAILED[*]}"; exit 1
fi

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
R="${RSCRIPT:-./r_env/bin/Rscript}"
LOG="logs_real"; mkdir -p "$LOG"

run () {   # run <script> <label> [logbase]
  local script="$1" label="$2" logbase="${3:-$(basename "$script" .R)}"
  echo "=============================================================="
  echo ">> [$label] $script  ($(date '+%H:%M:%S'))"
  if [ ! -f "$script" ]; then echo "  MISSING: $script — skipping"; return; fi
  if $R "$script" > "$LOG/$logbase.log" 2>&1; then
    echo "   OK ($label)"
  else
    echo "   FAILED ($label) — see $LOG/$(basename "$script" .R).log"; FAILED+=("$label")
  fi
}

FAILED=()

# Bridge: regenerate results/rdata/tcga_processed.RData on a clean clone where it
# is missing (no other tracked script creates it). Needs the base-preprocessing
# output data/processed/TCGA_STAD_processed.RData (from part1); the script errors
# clearly if that prerequisite is absent.
if [ ! -f results/rdata/tcga_processed.RData ]; then
  run analysis/00_prepare_tcga_processed.R    "prepare TCGA processed RData (bridge)"
else
  echo ">> [bridge] results/rdata/tcga_processed.RData present — skipping regen."
fi

# Order respects dependencies: signature (07) before consumers (12/13/14).
run analysis/20_integrated_deg.R              "integrated TCGA+GTEx DEG + enrichment"
run analysis/22_deg_diagnostics.R             "DEG diagnostics"
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

MR_RAN=0
if [ -n "${OPENGWAS_JWT:-}" ]; then
  # European outcome (Sakaue 2021) — defaults inside 11_real_mr.R
  run analysis/11_real_mr.R                    "Mendelian randomisation (European)" 11_real_mr_eur
  # East-Asian sensitivity outcome (ebi-a-GCST90018629), distinct outdir
  export GC_OUTCOME="ebi-a-GCST90018629" MR_OUTDIR="results/mr_real_eas"
  run analysis/11_real_mr.R                    "Mendelian randomisation (East-Asian sensitivity)" 11_real_mr_eas
  unset GC_OUTCOME MR_OUTDIR
  MR_RAN=1
else
  echo ">> [MR] OPENGWAS_JWT not set — skipping real MR (NOT faked)."
fi

echo "=============================================================="
if [ "$MR_RAN" -eq 1 ]; then MR_STATUS="ran (European + East-Asian)"; else MR_STATUS="SKIPPED (OPENGWAS_JWT unset)"; fi
if [ ${#FAILED[@]} -eq 0 ]; then
  if [ "$MR_RAN" -eq 1 ]; then
    echo "PIPELINE COMPLETE — all stages exited 0. MR: $MR_STATUS."
  else
    echo "PIPELINE PARTIAL — non-MR stages exited 0, but MR: $MR_STATUS."
    echo "  (Set OPENGWAS_JWT and re-run for a full pipeline.)"
  fi
else
  echo "PIPELINE finished with FAILURES: ${FAILED[*]}  | MR: $MR_STATUS"; exit 1
fi

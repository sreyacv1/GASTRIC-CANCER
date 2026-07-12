# Analysis pipeline

Every result in `PAPER.md` maps to a script under `analysis/` and an output under
`results/`. Run the whole pipeline with `bash run_real_pipeline.sh` (R 4.3.x; see
`package_versions.csv`). The Mendelian-randomisation stage needs an OpenGWAS token
(`export OPENGWAS_JWT=<token>`).

## Stages

| # | Script | Produces | Key inputs |
|---|--------|----------|-----------|
| base | `analysis/gastric_cancer_multiomics_v2_part1.R` | TCGA/GTEx/GEO processing + harmonisation | TCGA-STAD, GTEx, GEO |
| 20 | `analysis/20_integrated_deg.R` | integrated TCGA+GTEx tumour-vs-normal DEG + Hallmark/GO/KEGG enrichment; GEO cohorts as validation-concordance | harmonised TCGA+GTEx matrix |
| 09 | `analysis/09_functional_enrichment.R` | GO/KEGG/GSEA tables + plots | DEG tables |
| 08 | `analysis/08_immune_deconvolution.R` | MCP-counter/xCell scores, validated vs measured leukocyte % | TCGA expression |
| 07 | `analysis/07_external_validation.R` | 25-gene LASSO-Cox signature + ACRG validation | TCGA + GSE62254 |
| 12 | `analysis/12_multicohort_validation.R` | validation in GSE15459, GSE84437 | signature + GEO |
| 13 | `analysis/13_combined_nomogram_DCA.R` | combined model, DCA, time-AUC, NRI/IDI (in-sample) | signature + TCGA |
| 17 | `analysis/17_external_utility_ACRG.R` | external DCA/IDI/NRI/C-index in ACRG | signature + GSE62254 |
| 14 | `analysis/14_wgcna_real.R` | modules, module–trait, hub genes, survival module | TCGA expression |
| 18 | `analysis/18_wgcna_power_robustness.R` | module robustness across soft-powers 3–12 | TCGA expression |
| nomo | `analysis/nomogram_real_OS.R` | clinical survival nomogram (complete-case) | TCGA clinical |
| 19 | `analysis/19_nomogram_bootstrap_selection.R` | selection-inside-bootstrap C-index + EPV | TCGA clinical |
| 10 | `analysis/10_microbiome_robust.R` | diversity, PERMANOVA, differential abundance, co-abundance | tissue 16S (PRJDB20660) |
| 15 | `analysis/15_scrna_validation.R` | single-cell CAF localisation | GSE134520 |
| 16 | `analysis/16_drug_repurposing.R` | in-silico drug-repurposing (TCGA DEG) | DEGs + Enrichr/LINCS |
| 21 | `analysis/21_drug_repurposing_integrated.R` | drug-repurposing on the integrated TCGA+GTEx DEG | integrated DEG + Enrichr/LINCS |
| 11 | `analysis/11_real_mr.R` | two-sample Mendelian randomisation | IEU OpenGWAS |

## Mendelian-randomisation GWAS accessions (IEU OpenGWAS)

Exposures (European): anti–*H. pylori* IgG seropositivity `ebi-a-GCST90006910`
(Butler-Laporte 2020); genus *Streptococcus* `ebi-a-GCST90017070`, *Prevotella*
`ebi-a-GCST90017045`, *Veillonella* `ebi-a-GCST90017088`, *Lactobacillus*
`ebi-a-GCST90017030` (MiBioGen 2021); *Fusobacterium* `ebi-a-GCST90032406` (Qin 2022).
Outcome: gastric cancer, European `ebi-a-GCST90018849` (Sakaue 2021); East-Asian
sensitivity `ebi-a-GCST90018629`.

## Reproducibility

- R 4.3.3; package versions pinned in `package_versions.csv`.
- The MR token is read from `OPENGWAS_JWT` and is never stored in the repo.
- Large primary data are public and not committed; see `data/README.md`.

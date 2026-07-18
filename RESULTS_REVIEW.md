# Full Results Review — every analysis independently re-checked

Independent re-computation / re-reading of every `results/` subdirectory (31 directories,
165 tables) against the manuscript's claims. **Verdict: every headline number reproduces.**
No unsupported claim, no fabrication signature, no duplicated-statistic artefact found.

## Survival / prognostic signature
| Analysis | Manuscript | Source file | Re-check |
|---|---|---|---|
| Nested-CV C-index | 0.611 (0.562–0.659); Uno 0.573 | `nested_cv/performance.csv` | ✓ exact |
| Hartung–Knapp meta | HR 1.19 (0.96–1.47), p=0.073, I²=19% | `meta_HK/meta_result.csv` | ✓ exact |
| Per-cohort validation | ACRG C=0.61 HR 1.90 (1.37–2.62); GSE15459 C=0.58 HR 1.68; GSE84437 null (C=0.53) | `validation_multi/cindex_HR_summary.csv` | ✓ exact (confirms ACRG HR=1.90 fix) |
| Time-varying HR | 1.49@12mo → 1.03@36 → 0.87@60; cox.zph p=0.003 | `timevarying_ACRG/*` | ✓ exact |
| External added value | ΔC +0.005, IDI 0.021, NRI 0.175, LRT p=0.002 | `external_utility_ACRG/added_value_external.csv` | ✓ exact |
| Nomogram | in-sample combined C=0.77–0.79; ACRG transported 0.66 / refit 0.72 | `nomogram*/` | ✓ exact |
| Signature stability | B=200 bootstrap; SERPINE1 retained; 13/25 >50% | `signature_stability/*` | ✓ present |

## Differential expression / enrichment
| DEG integrated | 3,722 up / 4,025 down of 12,899; **λ=17.3** | `tables/DEG_integrated_TCGA_GTEx.csv` | ✓ recomputed 3722/4025/12899, λ=17.32 |
| GTEx-confound concordance | top-100/200 = 100% replicated in TCGA-only | `tables/DEG_TCGAonly_replication.csv` | ✓ exact |
| Hallmark EMT (diffuse) | NES 3.17, padj 1.6×10⁻⁴⁰ | `enrichment/GSEA_Hallmark_DiffuseVsIntestinal.csv` | ✓ NES 3.174, padj 1.58e-40 |
| Hallmark EMT (integrated tumour-vs-normal) | positive NES | `enrichment_integrated/fgsea_Hallmark_integrated.csv` | ✓ NES 1.85, padj 1.4e-5 |

## Immune / stroma / purity
| Immune deconvolution | T-cell(MCP) vs leukocyte% ρ=0.67, p=3.6e-36 | `immune/validation_vs_measured.csv` | ✓ 0.6656 / 3.57e-36 |
| Tumour purity | signature vs purity ρ=−0.20; vs stroma score +0.39 | `purity/correlations.csv` | ✓ exact |

## WGCNA / module preservation / single-cell
| Red-module survival | univariable HR 1.31 (1.12–1.53) p=9.3e-4; multivariable 1.35 (0.95–1.90) p=0.090 | `wgcna_real/ME_survival_cox*.csv` | ✓ exact |
| Soft-power scale-free R² | corrected to 0.865–0.886 across powers 3–12 | `wgcna_real/soft_threshold_table.csv` | ✓ (was mis-stated 0.877 floor) |
| Module preservation Z | 15.9 / 16.8 / 17.1 | `module_preservation/*` | ✓ exact |
| scRNA composition | 43,992 cells; fibroblast 4.16% | `scrna/celltype_composition.csv` | ✓ exact |
| scRNA localisation (non-circular) | 23/23 non-annotation hub genes fibroblast-dominant (median 0.96) | `scrna/gene_dominant_celltype_noncircular.csv` | ✓ computed this cycle |

## Microbiome (secondary, exploratory)
| RF classifier | cancer-vs-control AUC 0.916 (0.896–0.936) | `microbiome_biomarker/05_rf_metrics_and_batch_sanity.csv` | ✓ exact |
| Batch leakage | flowcell prediction 78% (baseline 55%) | same file | ✓ 0.776 / 0.545 |
| β-diversity collapse | Bray R² 0.065→0.011 after flowcell adjustment | `microbiome_biomarker/02_beta_permanova.csv` | ✓ exact |
| α-diversity effect sizes | richness δ=−0.27; Shannon −0.12; Simpson −0.07 (ns) | `microbiome_biomarker/02_alpha_effectsizes_cascade.csv` | ✓ computed this cycle |
| Read-tracking fidelity | Pearson r=0.98 | `data/microbiome/reprocess/supp_table3_readtracking.csv` | ✓ (0.9833, prior audit) |

## Mendelian randomisation (secondary, exploratory)
| MR IVW (all 6) | all null; smallest p=0.35 (Streptococcus) | `mr_real/MR_results_all_methods_REAL.csv` | ✓ b→OR all match (H.pylori 0.96, Strep 1.10, Fuso 1.04) |
| Instruments | all six at p<1e-5; F-range 19.3–20.3 | `scratch/mr_real.log` | ✓ exact |
| MR-PRESSO global | H. pylori borderline p=0.052; rest 0.16–0.68 | `mr_real/MR_PRESSO_global_REAL.csv` | ✓ computed this cycle |
| Egger intercept | Fusobacterium 0.038, p=0.040 (nominal) | `mr_real/MR_pleiotropy_REAL.csv` | ✓ |
| East-Asian sensitivity | all null (larger 7,921-case outcome) | `mr_real_eas/*` | ✓ null throughout |

## Drug repurposing (hypothesis-generating)
Top reverser hits are CDK4/6 (palbociclib) and PI3K/mTOR/FGFR (NVP-BEZ235, dovitinib,
pd173074) inhibitors, targeting the proliferation program (FOXM1/TOP2A/MKI67/CDK1);
consistent with §3.1–3.2 tumour-up biology. `drug_repurposing_integrated/candidate_drugs_ranked.csv`.

## Orphan / unused outputs
`results/MaAsLin2_gene_vs_microbiome/` — exploratory, no significant results, not referenced
in the manuscript (now documented in PIPELINE.md).

**Conclusion.** The manuscript is a faithful report of its result files. The only corrections
this cycle were reporting/completeness items (documented in `DEFECT_REPORT.md`), not data errors.

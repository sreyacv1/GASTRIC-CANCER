## Appendix: Data-Acquisition & Reproducibility Checklist

*Gastric-cancer multi-omics study. Every dataset used, its accession, access route,
version/build, the raw→processed entry point, and the pipeline script that consumes it.
Datasets are listed in the order the pipeline uses them. This lets a reader reacquire
every input from its primary repository and re-enter the pipeline at the correct stage.*

## A. Host transcriptome (differential expression, WGCNA, signature, immune, drug-repurposing)

- [ ] **TCGA-STAD** — primary tumour + solid-tissue normal RNA-seq.
  - Accession/portal: GDC / UCSC Xena (project TCGA-STAD).
  - Build: STAR gene counts, GENCODE v36; 448 samples (412 tumour, 36 normal).
  - Clinical: Lauren class, AJCC stage, grade, TCGA molecular subtype, driver mutations,
    TMB, measured leukocyte/lymphocyte %, overall survival.
  - Raw→processed entry point: count matrix → `results/rdata/tcga_processed.RData`.
  - Consuming scripts: DEG (`analysis/…DEG…`), WGCNA (`analysis/…wgcna…`),
    immune deconvolution (`08_immune_deconvolution.R`), enrichment,
    microbe-response GSEA (`31_microbe_response_enrichment.R`).
- [ ] **GTEx v10** — normal stomach reference (DE denominator).
  - Portal: gtexportal.org. Combined with TCGA normals for the integrated tumour-vs-normal contrast.

## B. External survival / DE validation cohorts (transcriptome)

- [ ] **GSE62254 (ACRG)** — 300 patients; Affymetrix. Primary external validation of the signature.
  - Access: GEO. Local: `data/geo/GSE62254.rda`. Consumer: `12_multicohort_validation.R`, `32_nested_cv_signature.R`.
- [ ] **GSE15459** — 191 patients. Access: GEO. Local: `data/geo/GSE15459_es.rds` + `GSE15459_outcome.xls`.
- [ ] **GSE84437** — 431 patients. Access: GEO. Local: `data/geo/GSE84437_es.rds`. (Reported negative — retained transparently.)
- [ ] **GSE27342, GSE63089** — independent DE-concordance cohorts. Access: GEO.
  Consumer: DEG cross-cohort concordance (`results/tables/DEG_GEO_*.csv`).

## C. Single-cell transcriptome (module localisation)

- [ ] **GSE134520** — premalignant-to-early gastric-cancer scRNA-seq.
  - Access: GEO. Processing: Seurat v4 (QC 200–6,000 genes, MT<20% → 43,992 cells;
    LogNormalize; 2,000 HVGs; PCA 30; SNN res 0.5; UMAP). 8 cell types by canonical markers.
  - Consumer: scRNA localisation → `results/scrna/`.

## D. Tissue microbiome (16S)

- [ ] **PRJDB20660 (Japan, discovery)** — 944 libraries, 16S V3–V4.
  - Access: DDBJ, **raw FASTQ**. Processing: DADA2 from raw → ASV table
    (897 samples × 314 genera after host/off-target removal).
  - Read-tracking validation: per-sample counts vs published Supplementary Table 3, Pearson r=0.98
    (`results/microbiome_biomarker/00_readtracking_concordance.csv`).
  - Consumer: `23_dada2_16S.R` → `24_microbiome_real.R`; reference DB `silva_nr99_v138.1`.
- [ ] **PRJNA641258 (Italy, validation)** — paired tumour/non-tumour, V3–V4. Access: NCBI/ENA. (Ravegnini 2020, PMC7766162.)
- [ ] **PRJNA413125 (Portugal, validation)** — gastritis vs carcinoma, V5–V6. Access: NCBI/ENA. (Ferreira 2018.)
  - **Batch caveat:** the discovery tumour/normal split is confounded with sequencing flowcell (§3.8);
    validation cohorts are the batch-clean test.

## E. GWAS summary statistics (Mendelian randomisation)

*All accessions independently verified against GWAS-Catalog records; raw `gwasinfo()`/lookup output
archived at `results/mr_real/gwas_catalog_verification.json`.*

- [ ] **Exposures (IEU OpenGWAS):**
  - `ebi-a-GCST90006910` — anti-*H. pylori* IgG seropositivity (Butler-Laporte 2020; European, PMID 33204752).
  - `ebi-a-GCST90017070` — genus *Streptococcus* (id.1853; MiBioGen/Kurilshikov 2021; European, PMID 33462485).
  - `ebi-a-GCST90017045` — genus *Prevotella 9* (id.11183; MiBioGen; European, PMID 33462485).
  - `ebi-a-GCST90017088` — genus *Veillonella* (id.2198; MiBioGen; European, PMID 33462485).
  - `ebi-a-GCST90017030` — genus *Lactobacillus* (id.1837; MiBioGen; European, PMID 33462485).
  - `ebi-a-GCST90032406` — *Fusobacterium* A abundance in stool (Qin 2022; European/Finland, PMID 35115689).
- [ ] **Outcomes:**
  - `ebi-a-GCST90018849` — gastric cancer, European, 1,029 cases (Sakaue 2021, PMID 34594039). **Primary.**
  - `ebi-a-GCST90018629` — gastric cancer, East-Asian, 7,921 cases (Sakaue 2021). Cross-ancestry power sensitivity only.
  - Access token: OpenGWAS JWT via `OPENGWAS_JWT` env var. Consumer: `11_real_mr.R`.
  - Harmonised per-SNP instrument tables archived (`results/mr_real/`) so exact SNP sets survive upstream reprocessing.

## F. Environment & reproducibility

- [ ] R 4.3.3; pinned package versions in `package_versions.csv`; full `sessionInfo()` with the code.
- [ ] Complete-case analysis throughout; no imputation.
- [ ] Bundled `r_env/` has a build-time hardcoded path — after relocation use the `R_HOME`-override
      invocation or a system R (documented in `PIPELINE.md`).
- [ ] Every `results/` output maps to its generating script in `PIPELINE.md`;
      figure→source mapping in `FIGURE_SOURCES.md`.

## G. Analyses NOT performed (stated to prevent over-reading)

- [ ] **No genuinely paired per-patient tissue-microbiome + host-transcriptome cohort exists** in this
      study (16S and transcriptome are from different, unpaired cohorts). A per-genus host-gene MaAsLin2
      association was run in a superseded pipeline on a bridged, non-verifiable sample set; its outputs are
      not reproducible from retained files, so it is excluded and no host-gene↔genus claim is made. The
      microbe-response result (§3.9) is a pathway-level GSEA on the TCGA tumour-vs-normal ranking only —
      hypothesis-generating context, not a host–microbiome axis.

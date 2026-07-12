# Gastric Cancer Multi-Omics

An integrative multi-omics analysis of gastric cancer (GC) combining host
transcriptomics, tumour-tissue microbiome, and human genetics (Mendelian
randomisation), with an externally-validated prognostic signature and a
single-cell-corroborated stromal/cancer-associated-fibroblast (CAF) programme.

All results derive from **real public data only** — no simulated or imputed
values. Every figure/number in the manuscript traces to a file under `results/`.

## Repository layout

```
├── README.md                 # this file
├── PAPER.md                  # the manuscript (IMRaD)
├── PIPELINE.md               # stage-to-script-to-output map
├── run_real_pipeline.sh      # one-command pipeline runner
├── package_versions.csv      # pinned R package versions (R 4.3.3)
├── analysis/                 # all analysis scripts + notebook
├── data/                     # data sources (README) + curated microbiome tables
└── results/                  # result tables and figures
```

See `PIPELINE.md` for the full stage-by-stage description and the GWAS accessions
used for Mendelian randomisation.

## Data

Primary data are public and are **not** committed (large / downloadable):
TCGA-STAD, GTEx stomach, GEO GSE27342/63089/62254/15459/84437, single-cell
GSE134520, tissue 16S PRJDB20660, and IEU OpenGWAS summary statistics.
The curated microbiome genus tables used by the pipeline are under
`data/microbiome/`.

## Reproduce

```bash
# R 4.3.x environment expected (see package_versions.csv)
export OPENGWAS_JWT=<your_token>   # for the Mendelian randomisation stage
bash run_real_pipeline.sh
```

## Analysis stages (real)

Differential expression (TCGA tumour vs normal) → functional enrichment (GO/KEGG/
GSEA) → immune deconvolution (validated against measured leukocyte fraction) →
LASSO-Cox prognostic signature + multi-cohort external validation → WGCNA
(power-robust CAF module) → combined clinical+genomic model with external
decision-curve analysis → single-cell CAF localisation → tissue-microbiome
compositional analysis → two-sample Mendelian randomisation (honest null) →
in-silico drug repurposing.

See `CANONICAL_PIPELINE.md` for the full stage-to-script-to-output map.

# Data sources

All primary data are public. Large files are not committed to the repository;
download them from the sources below into the paths indicated.

| Dataset | Modality | Source / accession | Path |
|---------|----------|--------------------|------|
| TCGA-STAD | RNA-seq + clinical | GDC / Xena | `data/host/` |
| GTEx v10 stomach | RNA-seq (normal) | gtexportal.org | `data/host/`, `data/gtex/` |
| GSE27342, GSE63089 | microarray | GEO | `data/host/`, `data/geo/` |
| GSE62254 (ACRG) | microarray + survival | GEO | `data/geo/` |
| GSE15459, GSE84437 | microarray + survival | GEO | downloaded at runtime |
| GSE134520 | single-cell RNA-seq | GEO | `data/scrna/` |
| PRJDB20660 | tissue 16S (genus table) | study supplementary | `data/microbiome/` |
| MiBioGen / Qin 2022 / Butler-Laporte / Sakaue | GWAS summary stats | IEU OpenGWAS | fetched via API |

The curated microbiome genus abundance tables actually consumed by the pipeline
(`otu_table.csv`, `taxonomy.tsv`, `metadata_microbiome.tsv`) are included under
`data/microbiome/`. Everything else is downloaded via `analysis/` scripts or the
sources above.

`download_all_data.py` (in `analysis/`) helps fetch the host transcriptomic data.

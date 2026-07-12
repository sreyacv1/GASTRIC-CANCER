#!/usr/bin/env python3
"""
Download all transcriptomics datasets for gastric cancer multi-omics pipeline.
- TCGA-STAD RNA-seq (STAR counts) via GDC API
- GTEx Stomach TPM from Google Storage
- GEO GSE27342 and GSE63089 via GEOparse
- TCGA-STAD clinical metadata via GDC API
"""

import os, sys, json, gzip, time, requests
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path("data/host")
BASE.mkdir(parents=True, exist_ok=True)

def log(msg):
    print(f"[DOWNLOAD] {msg}")

# ============================================================================
# 1. TCGA-STAD RNA-seq via GDC API
# ============================================================================
def download_tcga_stad():
    log("TCGA-STAD: Querying GDC API for STAR counts...")
    
    manifest_path = BASE / "gdc_manifest_tcga_stad.txt"
    data_dir = BASE / "tcga_raw"
    data_dir.mkdir(exist_ok=True)
    
    # Query TCGA-STAD STAR counts
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-STAD"]}},
            {"op": "in", "content": {"field": "data.data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "data.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "experimental_strategy", "value": ["RNA-Seq"]}}
        ]
    }
    
    # Get file UUIDs
    response = requests.post(
        "https://api.gdc.cancer.gov/files",
        json={"filters": filters, "fields": "file_id,file_name,cases.case_id", "size": 500},
        headers={"Content-Type": "application/json"}
    )
    
    if response.status_code != 200:
        log(f"GDC API query failed: {response.status_code}")
        return None
    
    data = response.json()
    files = data.get("data", {}).get("hits", [])
    log(f"Found {len(files)} TCGA-STAD STAR count files")
    
    if not files:
        return None
    
    # Build manifest for bulk download
    manifest_lines = ["id\tfile_name\tuuid\ttype"]
    case_ids = []
    for f in files:
        fid = f["file_id"]
        fname = f["file_name"]
        case_id = f["cases"][0]["case_id"] if f.get("cases") else "unknown"
        manifest_lines.append(f"{fid}\t{fname}\t{fid}\tgene_expression")
        case_ids.append(case_id)
    
    with open(manifest_path, "w") as f:
        f.write("\n".join(manifest_lines))
    
    log(f"Manifest saved: {manifest_path}")
    
    # Download files individually (GDC client not available)
    downloaded = {}
    for i, f in enumerate(files):
        fid = f["file_id"]
        fname = f["file_name"]
        case_id = f["cases"][0]["case_id"] if f.get("cases") else "unknown"
        
        out_path = data_dir / fname
        if out_path.exists():
            # Read the file
            pass
        else:
            url = f"https://api.gdc.cancer.gov/data/{fid}"
            try:
                resp = requests.get(url, stream=True, timeout=120)
                if resp.status_code == 200:
                    with open(out_path, "wb") as fh:
                        for chunk in resp.iter_content(chunk_size=8192):
                            fh.write(chunk)
                    log(f"  [{i+1}/{len(files)}] Downloaded {fname}")
                else:
                    log(f"  [{i+1}/{len(files)}] Failed: {resp.status_code}")
                    continue
            except Exception as e:
                log(f"  [{i+1}/{len(files)}] Error: {e}")
                continue
        
        # Parse the STAR counts file (TSV format)
        try:
            if out_path.exists() and out_path.stat().st_size > 1000:
                df = pd.read_csv(out_path, sep="\t", index_col=0)
                # Columns: gene_id, raw_count, normalized_count, length, etc.
                # We need raw counts
                if "unstranded" in df.columns:
                    counts = df["unstranded"]
                elif "raw_count" in df.columns:
                    counts = df["raw_count"]
                else:
                    counts = df.iloc[:, 0]
                downloaded[case_id] = counts
        except Exception as e:
            log(f"  Error parsing {fname}: {e}")
        
        if (i + 1) % 50 == 0:
            log(f"  Progress: {i+1}/{len(files)} files processed")
    
    if downloaded:
        # Build count matrix
        log(f"Building count matrix from {len(downloaded)} samples...")
        tcga_matrix = pd.DataFrame(downloaded)
        tcga_matrix.index.name = "gene_id"
        tcga_matrix.to_csv(BASE / "tcga_stad_star_counts.csv.gz", compression="gzip")
        log(f"TCGA-STAD matrix saved: {tcga_matrix.shape[0]} genes x {tcga_matrix.shape[1]} samples")
        return tcga_matrix
    else:
        log("No TCGA files downloaded")
        return None

# ============================================================================
# 2. GTEx Stomach TPM
# ============================================================================
def download_gtex_stomach():
    log("GTEx Stomach: Downloading TPM matrix...")
    
    tpm_url = "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"
    tpm_path = BASE / "GTEx_v10_tpm.gct.gz"
    
    if not tpm_path.exists():
        log(f"  Downloading GTEx v10 TPM (~1.5 GB)...")
        try:
            resp = requests.get(tpm_url, stream=True, timeout=600)
            if resp.status_code == 200:
                total = 0
                with open(tpm_path, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192*1024):
                        f.write(chunk)
                        total += len(chunk)
                        if total % (100*1024*1024) < 8192*1024:
                            log(f"    Downloaded: {total/1024/1024:.0f} MB")
                log(f"  GTEx TPM downloaded: {tpm_path.stat().st_size/1024/1024:.0f} MB")
            else:
                log(f"  GTEx download failed: {resp.status_code}")
                return None
        except Exception as e:
            log(f"  GTEx download error: {e}")
            return None
    else:
        log(f"  GTEx TPM already exists: {tpm_path.stat().st_size/1024/1024:.0f} MB")
    
    # Download sample attributes
    attr_url = "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_SampleAttributesDS.txt"
    attr_path = BASE / "GTEx_v10_attrs.txt"
    if not attr_path.exists():
        log("  Downloading sample attributes...")
        try:
            resp = requests.get(attr_url, timeout=120)
            if resp.status_code == 200:
                attr_path.write_bytes(resp.content)
                log("  Attributes downloaded")
        except Exception as e:
            log(f"  Attributes download error: {e}")
    
    # Parse TPM matrix - extract stomach samples
    log("  Parsing GTEx TPM matrix...")
    try:
        attrs = pd.read_csv(attr_path, sep="\t", low_memory=False)
        stomach_samples = attrs[attrs["SMTSD"] == "Stomach"]["SAMGID"].tolist()
        log(f"  Found {len(stomach_samples)} stomach samples")
        
        # Read GCT file
        with gzip.open(tpm_path, "rt") as f:
            # Skip first 2 lines (version + dimensions)
            f.readline()
            f.readline()
            # Read header
            header = f.readline().strip().split("\t")
            sample_cols = header[2:]  # Skip Name and Description
            
            # Find stomach sample indices
            stomach_idx = [i for i, s in enumerate(sample_cols) if s in stomach_samples]
            log(f"  {len(stomach_idx)} stomach samples in TPM matrix")
            
            if not stomach_idx:
                log("  WARNING: No stomach samples found in TPM matrix columns!")
                return None
            
            # Read data
            rows = []
            gene_ids = []
            for line in f:
                parts = line.strip().split("\t")
                gene_ids.append(parts[0])
                # Only extract stomach columns
                row_data = [float(parts[i+2]) for i in stomach_idx]
                rows.append(row_data)
        
        gtex_stomach = pd.DataFrame(
            rows,
            index=gene_ids,
            columns=[sample_cols[i] for i in stomach_idx]
        )
        gtex_stomach.index.name = "gene_id"
        
        # Log2 transform
        gtex_stomach = np.log2(gtex_stomach + 1)
        
        out_path = BASE / "gtex_stomach_tpm_log2.csv.gz"
        gtex_stomach.to_csv(out_path, compression="gzip")
        log(f"  GTEx Stomach saved: {gtex_stomach.shape[0]} genes x {gtex_stomach.shape[1]} samples")
        
        # Save sample metadata
        stomach_meta = attrs[attrs["SAMGID"].isin(stomach_samples)]
        stomach_meta.to_csv(BASE / "gtex_stomach_metadata.csv", index=False)
        log(f"  Stomach metadata saved: {stomach_meta.shape[0]} samples")
        
        return gtex_stomach
    except Exception as e:
        log(f"  GTEx parsing error: {e}")
        import traceback
        traceback.print_exc()
        return None

# ============================================================================
# 3. GEO Datasets via GEOparse
# ============================================================================
def download_geo_datasets():
    log("GEO: Downloading GSE27342 and GSE63089...")
    results = {}
    
    for gse_id in ["GSE27342", "GSE63089"]:
        log(f"  Downloading {gse_id}...")
        try:
            import GEOparse
            gse = GEOparse.get_GEO(geo=gse_id, destdir=str(BASE), silent=False)
            
            # Get expression matrix (first GPL)
            gpl = list(gse.gpls.keys())[0]
            matrix = gse.gpls[gpl].table.pivot_table(
                index="ID_REF", columns="GSM", values="VALUE"
            )
            matrix = matrix.apply(pd.to_numeric, errors="coerce")
            
            # Get phenotype data
            phenotype = pd.DataFrame({
                gsm: gse.gsms[gsm].metadata
                for gsm in gse.gsms.keys()
            }).T
            
            # Save
            matrix.to_csv(BASE / f"{gse_id}_matrix.csv.gz", compression="gzip")
            phenotype.to_csv(BASE / f"{gse_id}_phenotype.csv")
            
            results[gse_id] = {"matrix": matrix, "phenotype": phenotype}
            log(f"  {gse_id}: {matrix.shape[0]} genes x {matrix.shape[1]} samples")
        except Exception as e:
            log(f"  {gse_id} download error: {e}")
            import traceback
            traceback.print_exc()
    
    return results if results else None

# ============================================================================
# 4. TCGA Clinical Metadata
# ============================================================================
def download_tcga_clinical():
    log("TCGA-STAD: Downloading clinical metadata...")
    
    response = requests.post(
        "https://api.gdc.cancer.gov/cases",
        json={
            "filters": {
                "op": "in",
                "content": {"field": "cases.project.project_id", "value": ["TCGA-STAD"]}
            },
            "fields": "case_id,demographic,diagnoses.tumor_stage,diagnoses.primary_diagnosis,"
                      "diagnoses.age_at_diagnosis,family_histories.relative_with_cancer_history,"
                      "exposures.cigarettes_per_day,exposures.years_smoked,"
                      "diagnoses.classification_of_tumor,diagnoses.days_to_last_follow_up,"
                      "diagnoses.days_to_death,diagnoses.vital_status,"
                      "summary.data_categories.data_category,summary.data_categories.case_count,"
                      "submitter_id",
            "size": 500
        },
        headers={"Content-Type": "application/json"}
    )
    
    if response.status_code != 200:
        log(f"Clinical data query failed: {response.status_code}")
        return None
    
    data = response.json()
    cases = data.get("data", {}).get("hits", [])
    log(f"Found {len(cases)} TCGA-STAD cases with clinical data")
    
    # Flatten clinical data
    clinical_rows = []
    for c in cases:
        row = {"case_id": c.get("case_id"), "submitter_id": c.get("submitter_id")}
        
        demo = c.get("demographic", {}) or {}
        row.update({
            "gender": demo.get("gender"),
            "race": demo.get("race"),
            "ethnicity": demo.get("ethnicity"),
            "year_of_birth": demo.get("year_of_birth"),
        })
        
        for dx in c.get("diagnoses", []) or []:
            row.update({
                "tumor_stage": dx.get("tumor_stage"),
                "primary_diagnosis": dx.get("primary_diagnosis"),
                "age_at_diagnosis": dx.get("age_at_diagnosis"),
                "classification": dx.get("classification_of_tumor"),
                "days_to_last_follow_up": dx.get("days_to_last_follow_up"),
                "days_to_death": dx.get("days_to_death"),
                "vital_status": dx.get("vital_status"),
            })
            break
        
        clinical_rows.append(row)
    
    clinical_df = pd.DataFrame(clinical_rows)
    clinical_df.to_csv(BASE / "tcga_stad_clinical.csv", index=False)
    log(f"Clinical metadata saved: {clinical_df.shape[0]} cases")
    return clinical_df

# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    log("=" * 60)
    log("Starting full data download pipeline")
    log("=" * 60)
    
    # 1. TCGA-STAD
    tcga = download_tcga_stad()
    
    # 2. GTEx Stomach
    gtex = download_gtex_stomach()
    
    # 3. GEO
    geo = download_geo_datasets()
    
    # 4. Clinical
    clinical = download_tcga_clinical()
    
    log("=" * 60)
    log("Download summary:")
    log(f"  TCGA-STAD: {'SUCCESS' if tcga is not None else 'FAILED'}")
    log(f"  GTEx Stomach: {'SUCCESS' if gtex is not None else 'FAILED'}")
    log(f"  GEO: {'SUCCESS' if geo is not None else 'FAILED'}")
    log(f"  Clinical: {'SUCCESS' if clinical is not None else 'FAILED'}")
    log("=" * 60)

# The gastric tumour-tissue microbiome is confounded by sequencing batch: a re-analysis of the PRJDB20660 cohort with independent validation

*(Focused microbiome study — discovery in the Japanese gastric-tissue 16S cohort, validated in one independent batch-clean cohort. No other data used.)*

---

## Abstract

**Background.** Tumour-tissue microbiome studies increasingly nominate dysbiotic "biomarkers" of gastric cancer, but gastric tissue is a low-biomass niche in which technical artefacts can masquerade as biology. We re-analysed a large gastric-tissue 16S cohort from raw sequencing and asked whether its tumour-associated microbial signal is genuine and generalisable.

**Methods.** We reprocessed all 944 libraries of the Japanese gastric-tissue cohort **PRJDB20660** from raw FASTQ through DADA2 (V3–V4), recovering a genuine amplicon-sequence-variant (ASV) table (validated against the study's published read-tracking, Pearson r=0.98). The cohort spans the gastric-cancer sequence: non-ulcer (n=299) and ulcer (n=103) controls, cancer-adjacent normal (GCN, n=219) and cancer tumour (GCT, n=323). We assessed α/β-diversity, a cancer-vs-control classifier, and — critically — the relationship between phenotype and **sequencing flowcell**. We then tested whether the tumour signal replicates in one **independent, batch-clean** cohort (**PRJNA641258**, Italy; 20 gastric-cancer vs 20 matched controls, V3–V4, tumour and control sequenced together).

**Results.** Tumour (GCT) libraries were sequenced on **separate flowcells** from all non-tumour samples (of 216 tumour/adjacent-normal pairs, only 1 shared a flowcell): phenotype is nested within batch. A cancer-vs-control classifier appeared excellent (cross-validated AUC **0.916**) but was **batch-driven** — its top features were environmental contaminants, and the same features predicted the flowcell among biologically-similar samples at 78% accuracy (vs 54% baseline). β-diversity separation by phenotype was small and shrank further after modelling flowcell (Bray R² 0.065→0.011). In the **independent batch-clean cohort**, the tumour signal **did not replicate**: no diversity difference (Shannon p=0.25), no compositional difference (Bray PERMANOVA R²=0.018, **p=0.80**), oral taxa trending the opposite way, and **0/344 genera** significant.

**Conclusions.** The tumour-associated microbiome "signature" in this cohort is largely a **sequencing-batch artefact** and does not generalise. Only the less-confounded gastritis-to-cancer diversity gradient (which shares flowcells) is defensible. We present this as a cautionary, honest negative: **cross-cohort validation and explicit batch auditing are prerequisites** for any gastric-tissue microbiome biomarker claim.

---

## 1. Introduction

The gastric mucosa harbours a low-biomass microbial community, and its dysbiosis has been repeatedly proposed as a biomarker of gastric carcinogenesis — typically as an enrichment of oral-origin taxa (*Streptococcus, Fusobacterium, Peptostreptococcus, Veillonella*) accompanying reduced diversity. Because microbial reads in low-biomass tissue are vulnerable to reagent/kit contamination and sequencing-run batch effects, a differential-abundance signature between tumour and normal can arise from *how* samples were processed rather than *what* they contained. The decisive test of whether such a signature is biological is **replication in an independent cohort processed separately**: batch artefacts are cohort-specific and cannot transfer. We apply this test to a large gastric-tissue 16S cohort.

## 2. Methods

**Discovery cohort (Japan, PRJDB20660).** All 944 paired-end libraries (Illumina MiSeq, V3–V4; primers 341F/805R) were downloaded from DDBJ and processed through DADA2 (primer trimming, `truncLen` 260/220, error learning, ASV inference, chimera removal, SILVA v138.1 taxonomy). Host mitochondrial (~34% of reads), chloroplast and non-bacterial ASVs were removed; features were prevalence-filtered and agglomerated to genus (final 897 samples × 314 genera). Denoising fidelity was confirmed against the published per-sample read-tracking (Pearson r=0.98). Sample phenotype (Non-ul/Ul/GCN/GCT) and the **sequencing flowcell** of each library (from FASTQ headers) were recorded.

**Statistics.** α-diversity (Observed/Shannon/Simpson) was tested across the ordered cascade (Jonckheere–Terpstra, Kruskal–Wallis, Spearman). β-diversity (Bray–Curtis and CLR/Aitchison) was tested by PERMANOVA **with and without adjustment for flowcell**. A random-forest classifier (cancer vs control, CLR genus features, 5-fold CV) quantified apparent discrimination; a **batch sanity check** trained the same features to predict flowcell among biologically-similar samples.

**Validation cohort (Italy, PRJNA641258).** 40 gastric-biopsy libraries (V3–V4, MiSeq) — 20 gastric-cancer vs 20 matched controls, with tumour and control **sequenced on shared runs** so phenotype is *not* nested within batch — were processed through the identical DADA2/SILVA pipeline and analysed for the same tumour-vs-control contrasts.

## 3. Results

### 3.1 Phenotype is nested within sequencing batch (discovery cohort)

Cross-tabulating phenotype against flowcell exposed a near-total confound: tumour (GCT) libraries reside almost entirely on two flowcells (L3RVN, L848P) that contain essentially **no** control or adjacent-normal samples, while controls and adjacent-normal share three other flowcells.

| Flowcell | Non-ul | Ul | GCN | GCT |
|---|---|---|---|---|
| L848P | 0 | 0 | 1 | 200 |
| L3RVN | 0 | 0 | 0 | 95 |
| LJDKG | 144 | 74 | 103 | 0 |
| L7Y62 | 67 | 28 | 78 | 0 |
| L3RVR | 72 | 1 | 22 | 0 |

Of 216 tumour/adjacent-normal patient pairs, **only 1** has both tissues on the same flowcell. Tumour status therefore cannot be separated from batch by any statistical adjustment.

### 3.2 The apparent tumour signal is batch-driven

A cancer-vs-control classifier achieved a striking cross-validated **AUC = 0.916** (95% CI 0.896–0.936). However, its most important features were environmental/skin genera (*Dietzia, Serinicoccus, Methylobacterium, Sphingomonas*), and a sanity check confirmed the signal is technical: the same genus features predicted the **sequencing flowcell** among biologically-similar control/adjacent-normal samples at **77.6% accuracy** (majority baseline 54.5%; multiclass AUC 0.73). β-diversity separation by phenotype was modest and **halved after flowcell adjustment** (Bray R² 0.065→0.011; flowcell R² 0.030 exceeded phenotype R²). The tumour "oralization" pattern is thus largely an artefact of the batch structure in §3.1.

### 3.3 What is defensible: the less-confounded diversity gradient

Restricting to the flowcell-*sharing* groups (Non-ul → Ul → GCN), microbial richness **declined** monotonically along the gastritis-to-cancer-adjacent sequence (Observed Jonckheere–Terpstra Z=−4.94, **p=7.7×10⁻⁷**; Shannon Z=−2.66, p=0.008), and *Helicobacter* rose toward cancer-adjacent tissue — both consistent with prior gastric literature. These are the only signals not structurally confounded with batch.

### 3.4 The tumour signal does not replicate (independent batch-clean cohort)

In the Italian cohort, where tumour and control were sequenced together, the tumour-vs-control contrast was **null on every axis**:

| Test | Result |
|---|---|
| α-diversity (Observed / Shannon / Simpson) | p = 0.17 / 0.25 / 0.26 (no difference) |
| β-diversity (Bray PERMANOVA) | R² = 0.018, **p = 0.80** |
| β-diversity (Aitchison PERMANOVA) | R² = 0.025, p = 0.51 |
| Oral-taxa enrichment | trended **opposite** (lower in tumour) |
| Differential abundance | **0 / 344** genera significant (BH) |

Because per-run batch effects do not transfer between independent laboratories and countries, this clean non-replication shows that the discovery cohort's tumour-associated microbiome signature is **not a generalisable biological biomarker**.

## 4. Discussion

Re-analysing a widely-usable gastric-tissue 16S cohort from raw data, we find that its headline tumour-vs-normal dysbiosis is **inseparable from a sequencing-batch confound** and, tested directly, **fails to replicate** in an independent batch-clean cohort. The apparent AUC of 0.916 is a cautionary example: high classifier performance in a low-biomass tissue-microbiome study can reflect contamination and batch structure rather than biology, and only cross-cohort validation reveals which.

Our honest positive finding is narrow but real — a **declining-diversity gradient along the gastritis-to-cancer sequence** on the less-confounded axis, with *Helicobacter* enrichment toward cancer-adjacent tissue. We do **not** claim a tumour-specific oral-taxa biomarker on this evidence.

**Limitations.** The discovery cohort's batch structure is irreparable post hoc; only new sequencing with balanced designs could resolve the tumour contrast. The validation cohort is modest (n=40) and V3–V4 like discovery but from a different population; a larger confirmatory cohort would strengthen the negative. Gastric tissue is low-biomass and lacked negative controls, precluding formal `decontam` correction.

**Conclusion.** For gastric-tissue microbiome biomarkers, **batch auditing and independent cross-cohort validation are not optional** — without them, a compelling but spurious dysbiosis signature is easy to produce and hard to distinguish from biology. We report this analysis as an honest negative in that spirit.

---

*Data: discovery PRJDB20660 (DDBJ); validation PRJNA641258 (SRA). Reprocessing and analysis code: `analysis/23_dada2_16S.R`, `analysis/28_validation_IT.R`, `results/microbiome_biomarker/`.*

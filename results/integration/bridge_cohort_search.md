# Bridge-cohort search: gastric cohorts with same-patient microbiome + host RNA-seq

Goal: find a PUBLIC gastric-cancer cohort with BOTH microbial profiling
(16S rRNA or shotgun metagenomics) AND host RNA-seq measured on the SAME
patients, to correlate microbe abundance with host gene expression per patient.
All candidates below are independently-assayed dual-omics (NOT TCGA/TCMA-derived).

Verification method: for each BioProject, library strategies were checked
directly against the ENA filereport API
(`result=read_run&fields=library_strategy,library_source,sample_alias`).
Data-availability statements read from the primary papers.

Date: 2026-07-13.

---

## RANKED SHORTLIST (genuine same-patient dual-omics, verified)

### #1 — Frontiers Cell Infect Microbiol 2026 (Northwest China, advanced GC)
"Integrated multi-omics profiling ... paired tumor and peritumoral tissues
of advanced gastric cancer patients from Northwest China"
(10.3389/fcimb.2026.1763765)

- Accessions (both verified in ENA):
  - **PRJNA1067082** = shotgun metagenomics — ENA confirms all runs
    `library_strategy=WGS, library_source=METAGENOMIC` (human gut metagenome),
    ~162 runs.
  - **PRJNA1013242** = host transcriptome — ENA confirms all runs
    `library_strategy=RNA-Seq, library_source=TRANSCRIPTOMIC` (Homo sapiens),
    ~146 runs.
- n: 88 advanced GC patients. Paired tumoral + peritumoral mucosa per patient.
  Integrated subset with BOTH metagenome AND transcriptome:
  **53 peritumoral + 60 tumoral samples** (paper's own paired-omics count).
  73 tumoral + 73 peritumoral had transcriptome (58 patients fully paired T/N).
- Tumour vs normal: YES, explicit — tumoral vs peritumoral (>=5 cm margin).
- Microbe type: shotgun metagenomics (not 16S) — species/strain + function.
- Downloadable: YES, both projects public on SRA/ENA (fastq via ENA).
- Per-patient linkage: metagenome aliases `gastric _tumor_mucosa_N` /
  `gastric mucosa_normal_N`; RNA aliases `gastric _tumor_patient_N`.
  The numbering scheme differs between the two omics, so the microbe<->RNA
  per-patient map is NOT 1:1 by alias string — it must come from the paper's
  supplementary metadata (which exists, since the authors ran the integration).
- Why #1: gastric-cancer tumour tissue, clean paired T/N design, BOTH omics
  independently assayed and fully public + verified, largest confirmed
  paired-omics n. Best bridge cohort overall.

### #2 — iScience 2022 (Korea, Hanyang) — TRUE 16S + RNA-seq
"Multi-omics reveals microbiome, host gene expression, and immune landscape
in gastric carcinogenesis" (PMC8898972; 10.1016/j.isci.2022.103956)

- Accessions (both verified in ENA):
  - **PRJNA703470** = microbiome — ENA confirms `library_strategy=AMPLICON`
    (16S rRNA), 110 runs.
  - **PRJNA703469** = host transcriptome — ENA confirms
    `library_strategy=RNA-Seq, TRANSCRIPTOMIC` (Homo sapiens), 110 runs.
- n: 70 participants across a biopsy cohort (30) + surgery cohort (40).
  Both bacterial DNA (16S) and human mRNA taken from the SAME gastric antrum /
  surgical specimen per participant.
- Tumour vs normal: disease spectrum — healthy stomach, gastritis, gastric
  cancer; GC distinguishable from adjacent severe gastritis. Tumour tissue in
  the surgery cohort. (Not a strict paired tumour/adjacent-normal design.)
- Microbe type: **genuine 16S rRNA amplicon** (the exact "16S + RNA-seq" match).
- Downloadable: YES, both projects public on SRA/ENA.
- Per-patient linkage: 16S aliases `Hanyang__SRA_2020_<group>_<n>` (e.g. C_03,
  D_02); RNA aliases `Hanyang__SRA_rna_2020_<group>_<n>` (e.g. GC_30, A_10).
  Group-letter/number encodes disease group + subject; the 16S<->RNA map is
  recoverable but requires the paper's sample table, not a raw string match.
- Why #2: the only fully-public cohort here with real 16S (not shotgun) AND
  RNA-seq on the same subjects. Ranked below #1 only because it is a
  carcinogenesis spectrum (antrum biopsies) rather than a clean paired
  tumour/normal GC design.

### #3 — Cell Discovery 2024 (S. anginosus, 609 samples) — 16S public, RNA-seq NOT verifiable
"Characterization of the landscape of the intratumoral microbiota reveals that
Streptococcus anginosus increases the risk of gastric cancer initiation and
progression" (PMC11589709; 10.1038/s41421-024-00746-0)

- Accession: **PRJNA1061213** — ENA confirms 695 runs, ALL
  `library_strategy=AMPLICON / library_source=METAGENOMIC` (16S). Metabolome:
  MetaboLights **MTBLS9211**.
- n: 609 samples; 16S on 290 tumour + 319 adjacent-normal (274 paired patients).
  Transcriptome reported on 108 paired tumour/normal; metabolome on 90 paired.
- Tumour vs normal: YES, excellent — 16S sample aliases are `CT###` (cancer
  tumour) vs `NT###` (adjacent normal), so tumour/normal and patient index are
  explicit in the alias itself (cleanest linkage of the three for 16S).
- CAVEAT (blocks it from #1): the paper states the transcriptome is under the
  same PRJNA1061213, but ENA/SRA return ONLY the 695 amplicon runs there — NO
  human RNA-Seq is present under that accession, and no sibling human
  transcriptome BioProject was locatable (the nearest hit, PRJNA1103097, is an
  unrelated mouse Fusobacterium study). The 274-paired 16S is a first-class
  public resource; the 108-paired host RNA-seq is NOT verifiably downloadable
  ("all supporting data available upon request"). Usable as a bridge cohort
  only if the RNA-seq is obtained by request from the authors.

---

## Excluded (checked, failed a hard criterion)

- **PMC11456469** (Mucosal microbiota, 170 tissues, 85 paired GC patients,
  PRJNA1032279): host gene expression is **qRT-PCR only**, not RNA-seq /
  transcriptome. Fails criterion 2.
- **TCGA / TCMA**: microbes are extracted from TCGA RNA-seq, i.e. NOT an
  independently-assayed 16S. Explicitly out of scope (handled separately).
  None of the three cohorts above are TCGA-derived.

## What was searched
WebSearch (multiple phrasings: matched/paired/same-patient multi-omics,
16S/metagenomic + RNA-seq/transcriptome, integrative microbiome+transcriptome
gastric, GEO SuperSeries, dbGaP/SRA dual-strategy BioProjects); primary-paper
data-availability sections; ENA filereport API library-strategy verification
for every candidate BioProject; NCBI eutils for sibling/umbrella projects.

## Bottom line
Same-patient microbiome+RNA-seq gastric cohorts are rare but they DO exist
publicly. Two are fully verified end-to-end and immediately usable:
**#1 Frontiers 2026 (PRJNA1067082 shotgun + PRJNA1013242 RNA-seq, 88 AGC,
paired T/N)** and **#2 iScience 2022 (PRJNA703470 16S + PRJNA703469 RNA-seq,
70 subjects)**. For both, per-patient microbe<->gene linkage needs the paper's
supplementary sample table (the raw SRA aliases use different numbering per
omics). #3 Cell Discovery 2024 has the largest, cleanest public 16S
(PRJNA1061213, CT/NT-coded, 274 paired) but its host RNA-seq is not
verifiably public.

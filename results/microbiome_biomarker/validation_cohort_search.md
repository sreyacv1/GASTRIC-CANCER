# Independent gastric-tissue 16S validation-cohort search

Goal: find an INDEPENDENT public gastric-mucosa 16S cohort (cancer vs
non-cancer) to validate the oral-taxa dysbiosis biomarker (genus-level
enrichment of Streptococcus, Fusobacterium, Peptostreptococcus,
Prevotella, Veillonella, Parvimonas) discovered in DDBJ PRJDB20660
(Japan, V3-V4) — whose tumour-vs-normal contrast is fatally
batch-confounded (groups on separate flowcells).

Date: 2026-07-13. Verified live against ENA portal API + ENA/NCBI
BioSample XML. All FASTQ confirmed downloadable from ENA.

---

## TL;DR RECOMMENDATION

Validate on TWO complementary independent, non-Asian cohorts:

1. **PRJNA641258 (Italy, Bologna/Rimini — IJMS 2020, PMC7766162)** —
   PRIMARY / batch-safest. Gastric biopsies, **V3-V4 MiSeq 2x300 PE**
   (matches our discovery region), **PAIRED** tumour vs matched
   non-tumour (10 ADC + 10 SRCC tumour vs 20 matched controls). Labels
   in BioSample `host_disease` (CTRL / ADC / SRCC). 0.40 GB.
   The paired design is the strongest available structural guard
   against the exact batch-confound that killed PRJDB20660: a patient's
   tumour and matched control are essentially always sequenced in the
   same run. Downside: small n (40).

2. **PRJNA413125 (Portugal — Ferreira 2018, Gut)** — HIGHER-POWERED
   co-validation. 135 gastric samples: **81 Chronic_gastritis vs 54
   Gastric_carcinoma**. Labels in BioSample `isolation_source`
   (cleanest possible; counts match the paper exactly). V5-V6 Ion
   Torrent PGM, 7.14 GB. This is the canonical gastric-dysbiosis paper
   (it defined the Microbial Dysbiosis Index) — directly on-point.
   Downside: region V5-V6 (not V3-V4; genus-level still fine); the
   non-cancer arm is chronic gastritis, not adjacent-normal.

Use Italy as the batch-safe, region-matched, exact-contrast replicate;
use Ferreira for statistical power and the reference dysbiosis
comparison. Do NOT rely on FASTQ-header flowcell probing — SRA has
anonymized read headers for every SRA/ENA-hosted candidate (see
"Batch caveat" below).

---

## RANKED SHORTLIST (pass criteria 1-4)

### #1  PRJNA641258 — Italy (Bologna/Rimini)
- Study: Palmieri et al., IJMS 2020, "Gastric Adenocarcinomas and
  Signet-Ring Cell Carcinoma..." (PMC7766162).
- Tissue: gastric biopsy — `env_medium = Gastric Biopsy [NCIT:C51685]`,
  `env_broad_scale = stomach [UBERON:0000945]`. CONFIRMED gastric.
- Design: 40 samples = 10 ADC + 10 SRCC (tumour) + 10 CTRL-ADC +
  10 CTRL-SRCC (matched non-tumour). Paired tumour/normal.
- 16S region: **V3-V4**, Illumina MiSeq 2x300 PE.
- Downloadable: YES, 0.40 GB (fastq_ftp with bytes present).
- Label recovery: BioSample `host_disease` = {ADC, SRCC, CTRL ADC,
  CTRL SRCC}. Verified across all 40 samples (10/10/10/10).
- Batch: flowcell IDs NOT recoverable (SRA-anonymized headers), BUT
  paired case-control design + interleaved run accessions
  (SRR12072987 SRCC next to SRR12072988 CTRL_SRCC) strongly imply
  co-sequencing. Best structural guarantee against confounding.
- Independence: Italy — fully independent of Asian discovery cohort.
- Score: 9/10.

### #2  PRJNA413125 — Portugal (Ferreira 2018, Gut)
- Study: Ferreira et al. 2018, Gut, "Gastric microbial community
  profiling reveals a dysbiotic cancer-associated microbiota"
  (PMC5868293). Defined the Microbial Dysbiosis Index (MDI).
- Tissue: `organism = stomach metagenome`, geo Portugal:Porto.
  `source_material_id` e.g. "Gastritis63". CONFIRMED gastric.
- Design: 135 samples = **81 Chronic_gastritis + 54 Gastric_carcinoma**
  (verified by sweeping all 135 BioSamples; matches paper exactly).
- 16S region: V5-V6, Ion Torrent PGM (single-end).
- Downloadable: YES, 7.14 GB.
- Label recovery: BioSample `isolation_source` =
  {Chronic_gastritis, Gastric_carcinoma}. Cleanest binary label of any
  candidate; directly maps run->group.
- Batch: Ion Torrent (no Illumina flowcell structure); labels
  interleaved across run-accession order (no evidence of blocked
  batching); headers anonymized so not provable. Likely acceptable.
- Independence: Portugal — fully independent.
- Score: 9/10 (power + label cleanliness; -1 region mismatch).

### #3  PRJNA375772 — China (Xi'an), 454 GS
- Tissue: `env_material = human gastric tissue`. CONFIRMED gastric.
- Design: 389 runs; sample_title encodes group, e.g. "YJ_Ca6"
  (cancer), "YJ_N1" (normal), "HK12Dou243". Labels recoverable by
  PARSING sample_title (Ca vs N) — cryptic, needs decoding/paper.
- Region: 454 GS (pyrosequencing). Downloadable: YES, 2.93 GB.
- Batch: unassessed. Asian cohort (less independent than #1/#2).
- Score: 6/10 (label recovery requires title parsing; old 454 tech).

### #4  PRJEB26931 — China, gastric mucosa, MiSeq
- Tissue: `environment (material) = gastric mucosa`. CONFIRMED gastric.
- Design: 311 runs. BUT BioSample attributes are generic (no disease/
  tumour tag) — label recovery FAILS at the metadata level; would need
  the paper's supplementary table mapping runs->groups. Not yet
  located. Downloadable: YES, 1.05 GB.
- Score: 4/10 (criterion 4 unmet from metadata; supplement required).

---

## REJECTED

| Accession | Reason |
|-----------|--------|
| PRJNA678413 | `isolation_source = stool` — FECAL, not gastric (crit 1). |
| PRJNA239281 (Eun 2014) | labeled "human gut metagenome"; tiny n=31; 454. |
| PRJNA428883 | all runs map to one BioSample id; "gut metagenome" — unusable. |
| PRJEB11763 | 454 "gastric mucosal"; disease labels not verified in metadata. |
| PRJNA310127 (Yu 2017) | gastric but all 134 = GC, no non-cancer arm (crit 2). |
| PRJEB21497 (Castano 2017) | all 12 = GC, no controls (crit 2). |
| Coker 2018 (Gut) | no public accession — "available from author on request". |
| Sung 2019 | no accession; GC=0 (no cancer samples). |
| PRJNA532731 | already rejected (V4-V5, labels not in metadata). |
| PRJDB20660 | our confounded discovery cohort. |

---

## Batch caveat (applies to ALL candidates)

Criterion-5 probing via FASTQ header `INSTRUMENT:RUN:FLOWCELL:LANE`
is NOT possible for any SRA/ENA-hosted candidate: SRA re-headers reads
to `@<ACCESSION>.<n> <n>/1`, stripping the original instrument/flowcell
string. Verified empirically on PRJNA641258 (`@SRR12072987.1 1/1`) and
PRJNA413125 (`@SRR6151146.1 1/1`). Batch confounding must therefore be
argued structurally, not from headers:
- PRJNA641258: PAIRED tumour/normal per patient => same run => tumour/
  normal axis cannot be flowcell-confounded. Strongest guard.
- PRJNA413125: interleaved run accessions across groups, single-center
  Ion Torrent => no evidence of blocked batching.

---

## Source datasets

Candidate pool assembled from two gastric-microbiome meta-analyses:
- PMC9308235 (10 datasets): PRJEB11763, PRJEB21497, PRJEB26931,
  PRJNA239281, PRJNA310127, PRJNA375772, PRJNA428883, PRJNA532731,
  PRJNA641258, PRJNA678413.
- PMC9270228 (6 gastric-tissue studies): Coker_2018, Sung_2019,
  Ferreira_2018 (PRJNA413125), Yu_2017 (PRJNA310127),
  Eun_2014 (PRJNA239281), Castano-Rodriguez_2017 (PRJEB21497).

# NW China matched-cohort project — full shotgun-metagenomics + RNA-seq analysis blueprint

**Data:** PRJNA1067082 (shotgun metagenomics, ~162 runs) + PRJNA1013242 (host RNA-seq, ~146 runs); 88 advanced-gastric-cancer patients, paired **tumoral + peritumoral** mucosa; integrated paired-omics subset 53 peritumoral + 60 tumoral. Same patients → **true integration.** *(Separate project; this is the exhaustive menu — we will prioritise, not necessarily run all.)*

---

## PHASE 0 — Setup & QC

- **Reads QC/trim:** fastp/FastQC → adapter + quality trimming, dedup.
- **Host decontamination (metagenome):** map to human GRCh38 (Bowtie2/bbmap), keep unmapped = microbial reads. Report human fraction (gastric tissue is high-host).
- **RNA-seq QC:** fastp; rRNA/host handled at quantification.
- Reference DBs to stage: Kraken2 standard+fungi+viral (~100 GB), HUMAnN3 (ChocoPhlAn+UniRef), CARD, VFDB, GTDB, human genome+transcriptome index.

## PHASE 1 — Microbiome: taxonomy (all domains)

1. **Taxonomic profiling** — Kraken2+Bracken *and* MetaPhlAn4 (consensus); **bacteria, archaea, fungi, protozoa, viruses** at every rank (phylum→species→**strain**).
2. **Diversity:** α (Shannon, Simpson, richness, Faith's PD), β (Bray–Curtis, Aitchison/CLR, weighted UniFrac) — tumour vs peritumoral, PERMANOVA/PERMDISP.
3. **Differential abundance:** ANCOM-BC2 (+ subject random effect for pairing), MaAsLin2, ALDEx2, LEfSe — concordance across methods, BH-FDR.
4. **Community typing:** Dirichlet-multinomial mixture "microbial subtypes" (à la the base paper) + dysbiosis index.

## PHASE 2 — Virome (high value for GC)

5. **Human oncoviruses:** detect + quantify load of **EBV, HPV, HBV, HCV, CMV, HHV-8** (Kraken2 viral, bwa to viral refs). *EBV defines a TCGA GC molecular subtype* → correlate EBV⁺ with host program.
6. **Bacteriophages:** viral-contig calling (geNomad/VirSorter2/Cenote-Taker3), phage diversity, **phage–bacterial-host linkage**, tumour vs normal.
7. **EBV strain/typing** from reads where depth allows.

## PHASE 3 — Mycobiome (fungi)

8. Fungal community (Kraken2 fungi / dedicated) — *Candida, Malassezia, Aspergillus* etc.; tumour-vs-normal; fungal–bacterial correlation (trans-kingdom).

## PHASE 4 — Strain-level & microbial SNPs

9. **Strain tracking:** StrainPhlAn — strain identity, tumour↔peritumoral strain sharing within patient.
10. **SNV / microdiversity:** inStrain (per-species nucleotide diversity, SNVs, fixation), MIDAS2 — which strains/alleles associate with tumour.
11. **Pathogen genotyping from reads:** *H. pylori* **cagA / vacA** genotype, other virulence alleles.
12. **Microbial "GWAS":** associate species SNVs/strains with tumour status + host program.

## PHASE 5 — Functional & metabolic potential

13. **HUMAnN3:** gene families (UniRef90), **MetaCyc/KEGG pathways** (abundance + coverage), tumour-vs-normal DA.
14. **Targeted carcinogenic metabolism:** nitrosation, LPS biosynthesis, bile-acid metabolism, SCFA (butyrate/acetate) genes, lactate, polyamines.
15. **Genotoxins / oncogenic loci:** **pks island (colibactin)**, cytolethal distending toxin (**CDT**), *B. fragilis* toxin (**bft**) — direct carcinogenesis links.

## PHASE 6 — Resistome, virulome, mobilome

16. **Antibiotic-resistance genes:** CARD/RGI, AMRFinderPlus, ResFinder.
17. **Virulence factors:** VFDB.
18. **Mobile elements/plasmids:** plasmid detection, integrons, HGT signals.

## PHASE 7 — Assembly & novel genomes (MAGs)

19. Assembly (MEGAHIT/metaSPAdes) → binning (MetaBAT2+MaxBin2+CONCOCT, DAS_Tool) → QC (CheckM2) → dereplication (dRep) → taxonomy (GTDB-Tk) → annotation (Prokka/Bakta).
20. **Pangenome** of key species (e.g. *H. pylori*, streptococci); MAG-level functional/virulence annotation; novel-taxa discovery.

## PHASE 8 — Microbial ecology / networks

21. Co-abundance networks (SPIEC-EASI/SparCC), keystone taxa, modules; trans-kingdom (bacteria–fungi–phage) networks.
22. ML classifier (tumour vs peritumoral) on taxa/functions with nested-CV + feature importance + **batch/contamination audit** (the rigor we built here).

## PHASE 9 — Host transcriptome (matched RNA-seq)

23. Quantify (salmon/STAR) → DEG (tumour vs peritumoral), GSEA/pathways, **immune deconvolution** (CIBERSORTx/xCell/MCP), TCGA molecular-subtype assignment (EBV/MSI/GS/CIN).
24. **Host viral transcripts** (confirm EBV/virus *activity*, not just DNA presence).
25. Prognostic modelling if survival available; CAF/stromal program scoring.

## PHASE 10 — THE INTEGRATION (same-patient — the novelty)

26. **Per-patient microbe↔gene correlation:** Spearman/HAllA across all microbial features × host genes (FDR-controlled); sparse-CCA.
27. **Multi-omics factor models:** MOFA+ and mixOmics **DIABLO** (blocks: taxa, functions, virome, host expression) → co-varying multi-omics modules defining tumour state.
28. **Microbe→host-pathway associations:** MaAsLin2 with host pathway scores as outcome, microbial features as predictors.
29. **Mediation:** microbe → host inflammatory/CAF pathway → tumour phenotype (formal causal-mediation).
30. **Virome→host:** EBV load → host interferon/immune program; phage → bacterial-host → gene.
31. **Function-coupling:** microbial metabolic pathway ↔ host metabolic transcriptome (e.g. microbial SCFA genes ↔ host butyrate-response).
32. **Bipartite microbe–gene networks**; keystone-microbe → hub-gene mapping.
33. **Combined predictive model:** microbiome + host features → diagnostic/prognostic classifier, nested-CV, vs each omics alone (does integration beat single-omics?).

## PHASE 11 — Clinical & biological synthesis

34. Associate every layer with clinical variables (stage, Lauren, location, survival, Correa stage).
35. Trans-kingdom + host synthesis: a mechanistic model (dysbiosis/virome → inflammation → CAF/EMT → prognosis) — now testable **per patient**.
36. Reproducibility: full env lock, per-tool versions, reference-DB hashes, clean-clone runner.

---

## Rigor carried over from this project (non-negotiable)
- **Batch/contamination audit first** (sequencing-run structure, negative-control awareness, `decontam`/host-fraction).
- **Cross-validation & nested resampling** for any classifier/signature.
- **Independent validation** where possible (e.g. the Korea 16S+RNA-seq cohort, PRJNA703470/703469, as a second matched set).
- **Honest scoping** — report nulls; distinguish DNA presence from activity; don't over-read low-biomass signal.

## Realistic prioritisation (when we start)
**Tier 1 (core paper):** Phases 0,1,2,5,9,10 (taxonomy + virome + function + host + integration).
**Tier 2 (depth):** Phases 4,6,8 (strains/SNPs, resistome, networks).
**Tier 3 (if time/compute):** Phase 7 (MAGs), Phase 3 mycobiome deep-dive.

Feasibility check (data size, disk quota, tool/DB install) is **step 1** before any download.

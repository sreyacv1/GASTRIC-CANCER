# The Project End to End — Design, Biology, Measurement, and Full Pipeline

Companion to `LEARN_MY_PROJECT.md` (which teaches the concepts from DNA upward).
This document explains **why the study is designed the way it is**, **what is physically
being measured at each step**, **the biology being interrogated**, and **the complete
pipeline script by script**. Every number is from your own result files.

---

# PART A — THE STUDY DESIGN, AND WHY

## A.1 The two questions the project asks

The project is deliberately **two-armed**, and the arms have different evidential status.

**Arm 1 (primary, transcriptomic):** *Which gene-expression programmes in gastric tumour
tissue predict how long a patient survives, and which cells produce them?*

**Arm 2 (secondary, microbial/causal):** *Is the tumour-associated bacterial community
genuinely altered in gastric cancer, and if so, is it a cause or a consequence?*

Arm 1 produced a robust, externally-replicated finding. Arm 2 produced two
well-characterised negatives. The paper is honest that these are different in strength —
that asymmetry is the design working as intended, not a failure.

## A.2 Why comparative designs are needed at all

You cannot learn anything from tumour tissue alone. If a tumour expresses COL1A1 at
level 11.4, that number is meaningless in isolation. Meaning only comes from a
**contrast**: 11.4 in tumour *versus* 8.6 in normal stomach. Every finding in this
project is a contrast, and the quality of a finding depends almost entirely on **whether
the two groups being contrasted differ only in the thing you care about.**

That principle is what drives every design decision below.

---

# PART B — WHY PAIRED SAMPLES

This is the design question with the most subtlety, and the answer is **different in your
two arms**. It is worth understanding properly because it is exactly what a reviewer or
examiner will probe.

## B.1 What "paired" means

A **paired** (or matched) design compares two measurements **from the same individual** —
here, tumour tissue and adjacent normal tissue from the same patient's stomach. An
**unpaired** (independent-groups) design compares tumours from one set of people with
normal tissue from a different set of people.

## B.2 Why pairing is statistically powerful

Between-person variation in gene expression is **large**. Two healthy people can differ
substantially in the expression of the same gene because of genetics, age, diet,
inflammation, medication, and *H. pylori* status.

When you compare tumours from group A against normals from group B, that between-person
variability lands in your error term. It is noise you must overcome to detect the
tumour effect.

When you compare tumour vs normal **within the same person**, that person's genetic
background, age, diet, and infection history are **identical in both samples** — so they
cancel out. Each patient acts as their own control.

Concretely, an unpaired test asks "is mean(tumour) different from mean(normal)?" A paired
test computes the **within-patient difference** for each patient and asks "is the mean
difference different from zero?" The second question has far less noise in it.

The consequences:
- **Higher statistical power** — you detect real effects with fewer samples.
- **Automatic control of every stable patient-level confounder** — including ones you
  never measured and could not adjust for.
- **Cleaner biological interpretation** — a difference is attributable to the tissue's
  malignant state, not to who the donor was.

## B.3 What pairing does NOT fix

Pairing controls **patient-level** confounders. It does **not** control:
- **Technical/batch confounders** — if a patient's tumour and normal samples were
  processed on different sequencing runs, pairing does not save you (this is exactly the
  trap your microbiome arm fell into; §B.6).
- **Field effects** — "adjacent normal" tissue from a cancer patient is not truly
  healthy tissue. It sits in an organ with chronic inflammation and pre-malignant
  change (the Correa cascade). It may already carry early molecular alterations. This
  makes the tumour-vs-adjacent-normal contrast **conservative**: it can *understate* how
  different a tumour is from genuinely healthy stomach.

That second limitation is a real reason to *also* have unpaired healthy controls, which
brings us to your actual design.

## B.4 What your transcriptomic arm actually has

I checked your data directly rather than assuming. From
`data/processed/TCGA_STAD_processed.RData` (18,419 genes × 448 samples):

| | Count |
|---|---|
| Tumour samples | **412** |
| Normal samples | **36** |
| Unique tumour patients | 412 |
| Unique normal patients | 36 |
| **Patients with BOTH tumour and normal (true pairs)** | **33** |
| Normals with no matched tumour | 3 |

So TCGA-STAD contains **33 genuine tumour/adjacent-normal pairs**, but 412 tumours.

**This is the central design tension, and here is how your project resolves it.**

If you restricted the analysis to paired samples only, you would have a beautifully
controlled comparison — and you would throw away 379 tumours, i.e. most of your survival
information. You cannot build or validate a prognostic model on 33 patients.

So the design uses **the right structure for each question**:

1. **For differential expression** (tumour vs normal), the limiting factor is that 36
   normals is a small reference set. You therefore **augmented the normals with GTEx
   healthy stomach tissue** — genuinely healthy donors, not adjacent-normal. This
   increases power and escapes the field-effect problem in B.3.
2. **For the prognostic model**, you used **all 412 tumours** with survival data. This
   analysis does not need normals at all — it correlates expression *among tumours* with
   outcome.
3. **For the microbiome arm**, you used an explicitly **paired** differential-abundance
   analysis (§B.6).

**The trade-off you accepted, stated honestly:** adding GTEx means the normal group now
differs from the tumour group both biologically *and* by study of origin — an unpaired,
cross-study contrast. That inflates test statistics. You measured the inflation
(**λ = 17.3**, where 1.0 means none), ran a **permutation null** to see what λ arises by
chance (mean 1.06, 95th percentile 1.58), and confirmed the gene ranking against
**TCGA-only** analysis and two independent GEO cohorts (concordance r = 0.62–0.81).

That is the defensible way to make this trade: choose power, quantify the cost, disclose
it, and triangulate around it.

## B.5 So why not just use the 33 pairs as the main analysis?

Because **n = 33 answers a different, smaller question.** With 33 pairs you can establish
"these genes change in tumours" with excellent internal validity, but you cannot:
- build a 25-gene LASSO model (you would overfit catastrophically),
- estimate survival hazard ratios with usable precision,
- run WGCNA (co-expression networks need hundreds of samples to estimate correlations
  stably),
- or stratify by Lauren subtype or stage.

The 33 pairs are best used as a **confirmatory check** on the direction of key genes,
which is effectively what the TCGA-only concordance analysis provides.

## B.6 Pairing in the microbiome arm — and why it was not enough

Your microbiome design *does* use pairing. The differential-abundance analysis
`04b_DA_GCN_vs_GCT_paired.csv` is a **paired** comparison of gastric-cancer tumour tissue
(GCT) against that same patient's adjacent non-tumour tissue (GCN) — 61 genera tested,
18 significant.

This is the correct design choice: bacterial communities vary enormously between people
(far more than gene expression does), so pairing is even more valuable here.

**And yet the paired design was defeated by a technical confounder.** Here is the
crosstab of sample group against sequencing flowcell (`01_confound_crosstab_final.csv`):

| Group | KVF37 | L2HBK | L3RVN | L3RVR | L7Y62 | L848P | LJDKG |
|---|---|---|---|---|---|---|---|
| Non-ulcer control | 0 | 11 | 0 | 72 | 67 | 0 | 144 |
| Ulcer | 0 | 0 | 0 | 1 | 28 | 0 | 74 |
| Cancer adjacent-normal (GCN) | 0 | 0 | 0 | 22 | 78 | 1 | 103 |
| **Cancer tumour (GCT)** | **1** | **0** | **95** | **0** | **0** | **200** | **0** |

Look at the bottom row. **Every tumour sample sits on flowcells L3RVN or L848P, and
almost nothing else does.** The tumour samples and the control samples were sequenced on
essentially non-overlapping machinery runs.

This means **"tumour vs control" and "flowcell X vs flowcell Y" are the same comparison**
in this dataset. They are **perfectly confounded** — statistically inseparable, no matter
which test you use. Of 216 tumour/adjacent-normal pairs, **only 1 shares a flowcell**, so
even the paired design cannot rescue it: pairing removes the *patient* effect but leaves
the *run* effect intact.

Your evidence that this matters:
- A classifier separates cancer from control at **AUC 0.9161** — apparently excellent.
- The same data predicts **which flowcell** a sample came from at **77.6% accuracy**
  versus a **54.5%** majority baseline.
- Adjusting for flowcell **collapses** the community difference: Bray–Curtis
  **R² 0.0647 → 0.0109**.
- The top discriminating genera are **Dietzia, Serinicoccus, Methylobacterium,
  Microbacterium, Sphingomonas, Serratia, Cutibacterium** — environmental and skin
  contaminants, not gastric flora. Real gastric biology would be led by gut organisms.

**The lesson, and it is the most valuable methodological content in your paper:** a paired
design controls the confounder it was built to control, and nothing else. Always audit
batch structure explicitly. Most published tissue-microbiome studies do not run the
flowcell-predictability check you ran, which is why this literature is littered with
"dysbiosis" findings that may be reagent contamination.

## B.7 Design summary table

| Analysis | Design | Why |
|---|---|---|
| Integrated DEG | Unpaired, cross-study (412 TCGA tumours vs 36 TCGA + GTEx normals) | Needs many normals for power; escapes field effect. Cost: inflation λ=17.3, quantified and disclosed |
| TCGA-only DEG | Unpaired within one study (412 vs 36) | Concordance control for the above |
| Prognostic signature | Tumours only, n=412 → 383 with survival | Survival modelling needs no normals |
| External validation | Independent cohorts, n=922 | The only real test of generalisation |
| WGCNA | Tumours only, 5,000 most variable genes | Co-expression needs many samples |
| Microbiome DA | **Paired** (GCT vs same patient's GCN) | Bacterial communities vary hugely between people |
| Microbiome cascade | Unpaired across Correa stages | Different patients occupy different stages |
| MR | Genetic instruments, two-sample | Genotype is randomised at conception, immune to reverse causation |

---

# PART C — WHAT IS PHYSICALLY BEING MEASURED

A recurring source of confusion: each assay measures a **different molecule** and
therefore answers a different question. Here is what is actually in the tube.

## C.1 Bulk RNA-seq (TCGA) — measuring mRNA abundance

**Physical process:**
1. A tissue biopsy is homogenised — all cells destroyed together.
2. Total RNA is extracted; mRNA is captured (usually via its poly-A tail).
3. mRNA is reverse-transcribed to **cDNA** (more stable) and fragmented.
4. Adapters are ligated; fragments are amplified by PCR and loaded onto a flowcell.
5. The sequencer reads ~100-base stretches by detecting fluorescent nucleotides added
   one base at a time — millions of **reads** in parallel.
6. Reads are **aligned** to the human reference genome to determine origin.
7. Reads overlapping each gene are **counted**.

**What the number means:** read count for a gene ≈ how many mRNA molecules of that gene
were in the tissue. More count = more transcription.

**Normalisation is mandatory.** A sample with 30 million reads shows roughly double the
counts of a 15-million-read sample for *every* gene. Library-size correction (TPM/FPKM,
or DESeq2/limma size factors) removes this, and values are **log₂-transformed** because
expression spans several orders of magnitude and the log scale makes fold-changes additive
and variance more uniform.

**Critical limitation — this is the key to your primary finding:** bulk RNA-seq measures
the **average across every cell in the biopsy**. A biopsy is 40–80% cancer cells with the
rest fibroblasts, immune cells, endothelium. So a "high collagen tumour" could mean *more
fibroblasts* or *the same fibroblasts working harder* — bulk data cannot distinguish
these. This ambiguity is precisely why Parts F and G (deconvolution and single-cell)
exist in your project.

## C.2 Microarrays (your GEO cohorts) — measuring hybridisation intensity

A glass chip carries millions of short DNA **probes** of known sequence. Labelled sample
cDNA is washed over it; complementary sequences **hybridise** (A–T, G–C base pairing);
unbound material is washed away; a laser excites the label and a scanner records
**fluorescence intensity** per probe. More target present → brighter spot.

**Differences from RNA-seq that matter to you:**
- Arrays measure only genes with probes on the chip (**closed** platform); RNA-seq is open.
- Arrays measure **relative intensity**, not counts, with a compressed dynamic range and
  a fluorescence saturation ceiling.
- Absolute values are **not comparable** between platforms.

Your cohorts: TCGA = RNA-seq; ACRG/GSE62254 and GSE15459 = Affymetrix; GSE84437 =
Illumina array. **This is why you z-score within each cohort** before applying the
signature — converting each gene to "how many SDs above this cohort's mean," which is
comparable across platforms. It also guarantees no information leaks between cohorts.

## C.3 Single-cell RNA-seq (GSE134520) — measuring per-cell transcriptomes

**Physical process:** tissue is dissociated into a **suspension of individual cells**.
Each cell is captured in a droplet with a bead carrying a unique **cell barcode** and
**UMIs** (unique molecular identifiers). Inside the droplet the cell lyses and its mRNA is
reverse-transcribed with that barcode attached. All droplets are pooled and sequenced;
the barcode assigns each read back to its cell of origin.

**What you get:** a genes × **cells** matrix instead of genes × samples. Yours: **43,992
cells**, resolving **68.6% epithelial, 4.2% fibroblast** and 6 other types.

**Limitations:** very **sparse** (most genes show zero counts in any given cell — capture
efficiency is low); dissociation stresses cells and destroys spatial context; fragile
cell types are lost preferentially.

## C.4 16S rRNA amplicon sequencing — measuring bacterial gene copies

**Why 16S:** every bacterium carries the 16S ribosomal RNA gene. It contains regions
**conserved** across all bacteria (so one primer pair amplifies everything) interspersed
with **variable regions** (V3–V4, V5–V6) whose sequence differs between taxa. So a single
PCR gives you a taxonomic census.

**Physical process:** extract total DNA from tissue → PCR-amplify a variable region with
universal primers → sequence the amplicons → cluster/denoise into exact sequence variants
(**ASVs**) → match against a reference database (**SILVA v138.1**) to assign taxonomy.

**What the number means:** count of 16S gene copies assigned to a taxon. **Not** cell
count — bacteria carry 1–15 copies of the 16S gene depending on species, so abundance is
biased by copy number.

**Limitations that matter to your interpretation:**
- Usually only **genus**-level resolution — you cannot resolve species or strains.
- Measures **who is there**, not what they are doing (no functional information).
- Cannot distinguish live from dead bacteria.
- **Extremely vulnerable to contamination.** Tumour tissue has very low bacterial
  biomass, so reagent, kit, and skin contaminants can dominate. This is exactly what your
  random-forest importance ranking exposed.

## C.5 GWAS summary statistics — measuring genotype–phenotype association

Not a tissue assay. A **genome-wide association study** genotypes hundreds of thousands
of people at millions of **SNPs** (single-nucleotide polymorphisms — positions where the
population carries different bases) and, for each SNP, regresses the trait on the
genotype. The published **summary statistics** give per-SNP effect size (β), standard
error, and p-value.

You did not run a GWAS; you **reused published summary statistics** from IEU OpenGWAS as
inputs to Mendelian randomisation. Your outcome dataset: **1,029 gastric-cancer cases,
475,087 controls**.

## C.6 CRISPR knockout screens (DepMap) — measuring gene essentiality

**Physical process:** a pooled library of guide RNAs (targeting every gene) is delivered
to cancer cells expressing Cas9. Each cell receives one guide, which directs Cas9 to cut
that gene, disabling it. Cells grow for ~3 weeks. Guides are then sequenced. If a gene
was **essential**, cells carrying its guide died, so that guide is **depleted**.

**What the number means:** the Chronos/CERES **gene-effect score** — how much killing
that gene harms the cell, scaled so ~0 = no effect and ~−1 = a median common-essential
gene. Below about **−0.5** is called dependent.

This is **experimental** evidence, not computational prediction, which is why your DepMap
section materially strengthens the drug arm.

---

# PART D — THE BIOLOGY BEING INTERROGATED

## D.1 The organ

The stomach lining (**gastric mucosa**) is a specialised epithelium: mucus-secreting cells
protect against acid; **parietal cells** secrete hydrochloric acid; **chief cells**
secrete pepsinogen for protein digestion. Beneath the epithelium lies the **lamina
propria** — connective tissue with fibroblasts, immune cells, and vessels.

Normal gastric tissue therefore has a characteristic expression profile: high acid- and
digestion-related genes, high specialised metabolic genes. **Loss of this differentiated
programme is one of the two axes of your findings** — your down-regulated fatty-acid
metabolism and oxidative phosphorylation signals, and the "turquoise" WGCNA module
(r = −0.34 with tumour status), are this loss.

## D.2 How gastric cancer develops — the Correa cascade

Gastric adenocarcinoma is the endpoint of a decades-long inflammatory sequence:

**Normal mucosa → chronic gastritis → atrophic gastritis → intestinal metaplasia →
dysplasia → carcinoma**

- **Chronic gastritis** — persistent inflammation, usually from *H. pylori*, which
  survives gastric acid by producing urease and burrowing into the mucus layer.
- **Atrophic gastritis** — sustained inflammation destroys glands, including
  acid-producing parietal cells. Acid output falls (**hypochlorhydria**).
- **Intestinal metaplasia** — gastric epithelium is replaced by intestine-like
  epithelium: an adaptive but pre-malignant change.
- **Dysplasia** — cells become architecturally and cytologically abnormal.
- **Carcinoma** — invasion through the basement membrane.

**Why this matters for your microbiome arm — and it is important biology, not a
footnote.** Stomach acid is the principal barrier to bacterial colonisation. When
atrophic gastritis destroys acid production, the pH rises and **the stomach becomes
colonisable by bacteria that could not previously survive there.**

So a genuine reduction in *H. pylori* with an increase in oral- and gut-type organisms is
the *expected* consequence of the cascade. **Which means the direction of causality is
genuinely ambiguous on observational data**: altered flora could contribute to
carcinogenesis, or simply reflect the acid-free environment the disease created. This
biological ambiguity is the entire motivation for your Mendelian-randomisation arm.

Your diversity findings are consistent with the cascade: richness falls from a median of
**47 → 30** taxa (Cliff's δ = −0.27, p = 3.3×10⁻⁷), Shannon declines modestly
(2.176 → 1.923, p = 0.021), and **Simpson does not change significantly** (p = 0.172).
The honest biological reading: **rare taxa are lost while the dominant community
structure is largely preserved** — consistent with a niche becoming more restrictive,
not with a wholesale community replacement.

## D.3 The Lauren classification — two different diseases

- **Intestinal type** — cells retain adhesion and form gland-like structures; typically
  arises via the Correa cascade in older patients; better prognosis.
- **Diffuse type** — cells lose adhesion and infiltrate the wall individually; younger
  patients; worse prognosis; often driven by loss of **CDH1** (E-cadherin), the protein
  that physically glues epithelial cells together.

**Your data rediscovers this distinction from expression alone.** EMT
(epithelial–mesenchymal transition) is the top Hallmark set in diffuse tumours,
**NES +3.17, padj 1.6×10⁻⁴⁰**.

**The biology of EMT:** epithelial cells are normally polarised, adherent, and stationary.
EMT is a developmental programme (essential in embryogenesis and wound healing) in which
cells down-regulate adhesion molecules like E-cadherin, up-regulate mesenchymal markers
like vimentin and N-cadherin, dissolve their attachment to the basement membrane, and
become **migratory and invasive**. Cancer cells hijack it to metastasise.

So "EMT is enriched in diffuse tumours" is not an abstract statistic — it is the
molecular description of exactly what the pathologist sees down the microscope:
non-cohesive cells scattering through the gastric wall. **When a molecular result
independently recovers a century-old histological classification, that is strong evidence
your pipeline is measuring real biology.**

## D.4 The proliferation programme

Your integrated GSEA (`fgsea_Hallmark_integrated.csv`):

| Hallmark set | NES | Biological meaning |
|---|---|---|
| E2F targets | **+3.75** | Cell-cycle entry machinery up |
| G2M checkpoint | **+3.64** | Mitotic machinery up |
| MYC targets | **+2.67** | Master growth programme up |
| Oxidative phosphorylation | **−2.54** | Mitochondrial ATP production down |
| Fatty-acid metabolism | **−2.58** | Differentiated metabolic function down |

**The cell cycle, briefly.** A cell divides through G1 → S (DNA replication) → G2 → M
(mitosis). Progression is gated by **checkpoints**. The **RB–E2F axis** controls the
G1→S transition: RB protein restrains **E2F** transcription factors; when
**CDK4/CDK6–cyclin D** phosphorylates RB, E2F is released and switches on the genes
needed to replicate DNA. **This is why E2F-target enrichment is the transcriptional
fingerprint of proliferation — and why CDK4/6 emerge as drug targets in your DepMap
analysis (CDK4 −0.825, CDK6 −0.548, both dependent).** The pathway analysis and the
dependency screen are pointing at the same biology from two directions.

**MYC** is a master transcription factor driving ribosome biogenesis, metabolism, and
growth; it is among the most frequently activated oncogenes in human cancer.

**The Warburg effect.** Down-regulated oxidative phosphorylation alongside up-regulated
proliferation looks paradoxical — dividing cells need energy. The resolution: tumours
shift toward **aerobic glycolysis**, which yields far less ATP per glucose but produces
carbon **intermediates** for building nucleotides, amino acids, and lipids. A dividing
cell needs *biomass* more than it needs maximal ATP efficiency. Your data shows the
transcriptional signature of that trade.

## D.5 The tumour microenvironment and CAFs — your primary finding

A tumour is an **ecosystem**, not a clone. Its non-cancer components:

- **Cancer-associated fibroblasts (CAFs)** — the focus of your finding.
- **Immune cells** — T cells, macrophages, often re-educated to support the tumour.
- **Endothelial cells** — new blood vessels (**angiogenesis**).
- **Extracellular matrix (ECM)** — the structural scaffold of collagens, fibronectin,
  proteoglycans.

**What fibroblasts normally do:** they are the connective-tissue maintenance cells. They
synthesise and remodel ECM, and during wound healing they activate, proliferate, secrete
collagen, and contract to close the wound — then normally stand down.

**What CAFs are:** fibroblasts recruited and permanently activated by tumour-derived
signals (notably **TGF-β**, PDGF, FGF). They behave like wound-healing fibroblasts that
never switch off — hence the description of a tumour as "a wound that does not heal."

**What CAFs do that worsens prognosis — the mechanism behind your numbers:**
1. **Build a dense collagen scaffold (desmoplasia)** that stiffens the tissue. Mechanical
   stiffness itself promotes malignant behaviour through mechanotransduction.
2. **Secrete matrix-degrading enzymes** (MMPs) that cut paths through basement membrane —
   creating **invasion tracks**.
3. **Physically impede drug delivery.** A stiff, high-interstitial-pressure stroma
   compresses vessels and blocks drug penetration — a direct route to treatment resistance.
4. **Exclude immune cells.** Dense collagen forms a barrier that keeps cytotoxic T cells
   from reaching tumour cells — a major mechanism of immunotherapy resistance.
5. **Secrete pro-tumour signalling factors** (HGF, IL-6, CXCL12) that drive proliferation
   and survival.

**So a stroma-rich tumour is a more invasive, more drug-resistant, more
immune-evasive tumour.** That is *why* the stromal programme predicts survival.

**Your WGCNA red module is precisely this programme.** Its hub genes:

| Gene | Protein | Function |
|---|---|---|
| **COL1A1, COL1A2, COL3A1** | Collagens I, III | The main fibrous ECM scaffold |
| **POSTN** | Periostin | Matrix protein; promotes invasion and metastatic niche |
| **FAP** | Fibroblast activation protein | The canonical CAF marker; a protease |
| **FN1** | Fibronectin | Adhesion scaffold guiding migration |
| **LUM, DCN, BGN, VCAN** | Proteoglycans | Organise collagen fibril assembly |
| **SPARC** | Osteonectin | Matrix remodelling, collagen binding |
| **THBS2** | Thrombospondin-2 | Matrix-cell signalling |
| **CDH11** | OB-cadherin | Mesenchymal adhesion |
| **MATN3** | Matrilin-3 | ECM structural protein |
| **SERPINE1** | PAI-1 | Regulates matrix degradation and invasion |

This is a coherent, single-biology module: **activated fibroblasts building and remodelling
extracellular matrix.**

**And note the crucial honest detail in your own results:** in multivariable Cox adjusting
for stage, age, leukocyte fraction and Lauren subtype (n = 183, 57 events), the red module
attenuated to **HR 1.35 per SD (0.95–1.90), p = 0.090 — no longer significant**, while
stage remained dominant (HR 1.82). The biological interpretation is that **desmoplastic
content increases with invasion depth**, so the stromal programme substantially *tracks*
stage rather than adding information independent of it. Your paper says this. It is the
correct reading, and it explains why the clinical added value is small (ΔC +0.005).

## D.6 Tumour purity — why this is not an artefact

If the prognostic signal comes from stroma, a sceptic asks: are you just measuring **how
much tumour was in the biopsy**? A sample with few cancer cells has proportionally more
stroma, so a "stromal score" might be a purity proxy.

You tested this. The signature **survived adjustment for ABSOLUTE-estimated purity
(HR 2.97, 95% CI 2.28–3.86, p = 7.6×10⁻¹⁶), and purity itself was non-significant
(p = 0.35)**. Correlation with purity was modest and negative (ρ = −0.20), and with xCell
stromal score positive (ρ = 0.39). CAF score differed strongly across Lauren subtypes
(Kruskal p = 1.9×10⁻⁸).

**Conclusion:** the signal is stromal *biology*, not a sampling artefact. Anticipating and
closing this objection before a reviewer raises it is exactly what a strong paper does.

## D.7 The immune microenvironment

Immune infiltration cuts both ways: cytotoxic **CD8+ T cells** can kill tumour cells, but
**tumour-associated macrophages** are frequently polarised to an M2-like,
tissue-remodelling, immunosuppressive state that supports tumour growth.

Your deconvolution found a **macrophage-weighted** microenvironment, consistent with the
CAF/desmoplasia biology — a stroma-rich, immune-excluded, macrophage-rich tumour.

You validated the estimates against an independent non-computational ground truth:
predicted T-cell content vs **pathologist-measured leukocyte percentage**,
**Spearman ρ = 0.6656, p = 3.57×10⁻³⁶ (n = 272)**. And you disclosed the limitation:
CD8 estimates track leukocyte percentage (ρ = 0.47) but **not** lymphocyte-infiltration
percentage (ρ = 0.02, p = 0.74) — so the deconvolution captures **overall immune burden
rather than lymphocyte subset composition.**

---

# PART E — THE FULL PIPELINE, SCRIPT BY SCRIPT

All scripts live in `analysis/`. R 4.3.3; `sessionInfo.txt` and `package_versions.csv`
are committed. Order below is logical dependency order.

## Stage 0 — Data acquisition and harmonisation

**`00_prepare_tcga_processed.R`**, **`gastric_cancer_multiomics_v2_part1.R`**

- Download TCGA-STAD expression + clinical + survival; GTEx stomach; GEO series.
- Map probe/gene identifiers to common symbols (`gene_map`).
- Filter low-expression genes; VST/log-transform; assemble the harmonised matrix.
- **Output:** `data/processed/TCGA_STAD_processed.RData` — **18,419 genes × 448 samples**
  (412 tumour, 36 normal, 33 true pairs).

**Why first:** every downstream stage reads this object, so identifier harmonisation and
normalisation must be settled once, centrally.

## Stage 1 — Differential expression

**`20_integrated_deg.R`** (primary), **`22_deg_diagnostics.R`** (diagnostics)

- limma linear model, tumour vs normal, on the integrated TCGA+GTEx matrix.
- Empirical-Bayes variance shrinkage; **Benjamini–Hochberg FDR** across all genes.
- **Output:** **3,722 up / 4,025 down** of 12,899 tested. TCGA-only contrast:
  2,134 up / 2,362 down of 21,446.
- **Diagnostics:** inflation **λ = 17.3**; **permutation null** (100 label shuffles) →
  mean 1.06, 95th pct 1.58; cross-cohort concordance r = 0.62–0.81.

**Why the diagnostics exist:** a cross-study contrast inflates statistics. Rather than
ignore it, you quantified it, benchmarked it against a null, and confirmed the ranking
independently.

## Stage 2 — Functional interpretation

**`09_functional_enrichment.R`**

- **fgsea** on the full ranked gene list (Hallmark, GO, KEGG, C2CP) — no cutoff needed.
- ORA as a secondary view.
- **Output:** E2F +3.75, G2M +3.64, MYC +2.67, OXPHOS −2.54, FAO −2.58; EMT-diffuse
  +3.17 (padj 1.6×10⁻⁴⁰).
- **Disclosed caveat:** prominent olfactory/sensory-perception ORA terms reflect the large
  olfactory-receptor gene family, not gastric biology.

## Stage 3 — Prognostic signature

**`07_external_validation.R`** → **`12_multicohort_validation.R`** →
**`32_nested_cv_signature.R`** → **`26_signature_stability.R`** → **`34_meta_HK.R`** →
**`33_timevarying_ACRG.R`**

1. **Fit:** LASSO-Cox on TCGA tumours; λ by cross-validation → **25 genes**
   (16 risk-increasing, 9 protective).
2. **Honest internal estimate (`32_`):** **nested** CV — gene selection and λ tuning
   inside the inner loop only. **Harrell C = 0.6112 (0.5620–0.6589)**,
   Uno C = 0.5727, IBS 0.1789. Apparent C was **0.72**; the **~0.11 gap is leakage
   optimism**.
3. **External validation (`12_`):** z-score within each cohort, then apply the frozen model.

   | Cohort | n (events) | C | HR high-vs-low, median split | HR per SD (age/stage-adj) |
   |---|---|---|---|---|
   | ACRG/GSE62254 | 300 (152) | 0.608 | 1.90 (1.37–2.62) | 1.30 (1.10–1.54) |
   | GSE15459 | 191 (95) | 0.575 | 1.68 (1.11–2.54) | 1.20 (0.97–1.48) |
   | GSE84437 | 431 (207) | 0.530 | 1.11 (0.84–1.46) | 1.11 (0.97–1.27) |

   **Always state which HR you mean** — the median-split and per-SD numbers are different
   quantities.
4. **Stability (`26_`):** 200 bootstraps. **13/25 genes selected in >50%, none >80%**
   (median 0.505). So the *programme* is robust but the *specific gene list* is not — a
   limitation you state.
5. **Meta-analysis (`34_`):** REML + **Hartung–Knapp** on the per-SD adjusted estimates →
   **pooled HR 1.19 (0.96–1.47), p = 0.073, I² = 19.2% — not significant.**
6. **PH assumption (`33_`):** `cox.zph` violated (p = 0.0018 signature alone; p = 0.003
   age/stage-adjusted). Time-varying model: **HR 1.49 at 12 months → 1.03 at 36 → 0.87 at
   60**. So it is an **early**-risk marker, not a durable gradient.

## Stage 4 — Clinical utility

**`13_combined_nomogram_DCA.R`** (in-sample), **`17_external_utility_ACRG.R`**
(external), **`nomogram_real_OS.R`**, **`19_nomogram_bootstrap_selection.R`**

- Clinical model (age + stage) vs clinical + signature: ΔC, likelihood-ratio test,
  IDI/NRI, **decision-curve analysis**.
- **External result:** ΔC **+0.0045**, LRT p = 0.0020, IDI@3y 0.021, cNRI@3y 0.175.
- **Interpretation:** statistically detectable, **clinically negligible**. DCA shows no
  *incremental* net benefit over staging.

## Stage 5 — Network analysis (the primary finding)

**`14_wgcna_real.R`** → **`18_wgcna_power_robustness.R`** → **`25_module_preservation.R`**

1. **Build:** 5,000 most variable genes, signed-hybrid, **bicor**, soft power **β = 3**
   (scale-free R² = 0.88; near-flat 0.865–0.886 across powers 3–12, so chosen on the
   plateau/connectivity criterion). Topological overlap → clustering → 10 modules →
   eigengenes.
2. **Module–trait:** turquoise most associated with tumour status (r = −0.34,
   p = 1.8×10⁻¹³); **red** most prognostic (univariable **HR 1.31/SD, 1.12–1.53,
   p = 9.3×10⁻⁴**); red **attenuates to p = 0.090** after adjustment (§D.5).
3. **Power robustness (`18_`):** across β = 3, 6, 9, 12 the 12 CAF hub genes stay
   co-clustered (11–12/12) in a module that stays prognostic (HR/SD 1.28–1.31, all
   p < 0.0025). Only the colour label changes.
4. **Preservation (`25_`) — the strongest result in the project:**

   | Cohort | Zsummary | Eigengene HR/SD | p |
   |---|---|---|---|
   | ACRG/GSE62254 | **15.853** | 1.274 (1.091–1.487) | 0.0022 |
   | GSE15459 | **16.797** | 1.548 (1.248–1.919) | 6.8×10⁻⁵ |
   | GSE84437 | **17.083** | 1.237 (1.080–1.416) | 0.0021 |

   All Z > 10 = **strong preservation**. Note it replicates in **GSE84437, where the
   25-gene signature failed** — the biology is more robust than the gene score. This is
   why the paper is framed on the programme, not the signature.

## Stage 6 — Cellular attribution

**`08_immune_deconvolution.R`** → **`15_scrna_validation.R`**

1. **Deconvolution:** MCP-counter, xCell, ESTIMATE on bulk TCGA. Validated against
   **pathologist-measured leukocyte %**: ρ = 0.666, p = 3.6×10⁻³⁶ (n = 272), with the
   CD8/lymphocyte-subset limitation disclosed.
2. **Single cell (GSE134520):** Seurat QC (min 3 cells/gene, 200–6,000 features/cell,
   <20% mitochondrial), 2,000 HVGs, 30 PCs, resolution 0.5 → **43,992 cells**, 8 types
   (68.6% epithelial, 4.2% fibroblast).
3. **Non-circular localisation — the step that makes this evidence:** using COL1A1 to
   *label* a cluster "fibroblast" and then declaring COL1A1 fibroblast-specific proves
   nothing. So all annotation genes were **removed** and only the remainder tested.
   **Result: 23/23 remaining hub genes fibroblast-dominant, median fraction 0.960.**

Three independent methods — bulk co-expression, deconvolution, single-cell — converge on
fibroblasts. That convergence is what makes the finding credible.

## Stage 7 — Microbiome

**`23_dada2_16S.R`** → **`24_microbiome_real.R`** → **`28_validation_IT.R`** /
**`29_validation_PT.R`** → **`31_microbe_response_enrichment.R`**

1. **DADA2 (`23_`):** raw FASTQ (944 libraries, **13,487,331 reads**) → primer trim →
   truncLen 260/220 → error model → ASV inference → chimaera removal → **SILVA v138.1** →
   drop mitochondria (786)/chloroplast (82)/non-bacterial (111) → **897 samples,
   314 genera**. **Reproducibility check: Pearson r = 0.983** against the original
   study's published read-tracking table.
2. **Diversity and composition (`24_`):** alpha (Observed/Shannon/Simpson) with **Cliff's
   delta**; beta via Bray–Curtis **and Aitchison** (compositionally aware); PERMANOVA;
   **CLR** transform for differential abundance. Paired DA (GCT vs same patient's GCN):
   61 genera tested, 18 significant. Unpaired control-vs-GCN: 61 tested, 44 significant.
3. **Batch audit (`24_`) — the key result:** cancer-vs-control AUC **0.9161**, but
   flowcell predictable at **77.6%** vs **54.5%** baseline; flowcell adjustment collapses
   Bray R² **0.0647 → 0.0109**; top discriminators are environmental/skin contaminants.
   Only 1 of 216 pairs shares a flowcell → batch and biology **inseparable**.
4. **Independent replication:** Italy (`28_`) **null** (Bray R² 0.018, p = 0.80);
   Portugal (`29_`) **replicates reduced diversity** (Shannon p = 0.0044; Bray R² 0.145,
   p = 0.001). So the diversity finding has external support; the composition finding
   does not.

## Stage 8 — Causal inference

**`11_real_mr.R`**

- Two-sample MR, **TwoSampleMR**. Instruments: p < 5×10⁻⁸ with a 1×10⁻⁵ cascade,
  clumped r² = 0.001 / 10,000 kb; harmonise action 2.
- **Instrument strength:** F = 19.3–22.8 (all > 10) ✓
- Six exposures → gastric cancer (EUR **1,029 cases / 475,087 controls**; East-Asian
  sensitivity cohort).
- **Full sensitivity suite:** IVW, MR-Egger, weighted median, weighted mode, **MR-PRESSO**,
  leave-one-out, Cochran's Q.
- **Result: all six IVW null; smallest p = 0.348** (*Streptococcus*). MR-PRESSO global
  non-significant for **all six** (p = 0.052–0.676; *H. pylori* is the closest at 0.052,
  still above the 0.05 threshold).
- **Correct interpretation:** at this instrument and outcome power, **no evidence of a
  causal effect** — not evidence of no effect. Disclosed tension: *H. pylori* is an
  established cause epidemiologically, but the instrument proxies **antibody response**,
  not lifetime infection burden.

## Stage 9 — Therapeutic hypotheses

**`16_drug_repurposing.R`**, **`21_drug_repurposing_integrated.R`**,
**`26_depmap_dependency.R`**, **`30_basepaper_5gene.R`**

1. **Signature reversal** against LINCS/Enrichr → PI3K/mTOR inhibitors (NVP-BEZ235,
   PI-103, GDC-0941, torin-2), CDK4/6 (palbociclib), FGFR (PD-173074, dovitinib).
   **Caveat stated:** the up-signature is proliferation-dominated, so hits are broadly
   anti-proliferative rather than gastric-specific, and unvalidated experimentally.
2. **DepMap (`26_`)** — 35 gastric lines, CRISPR gene-effect:
   **PIK3CA −0.742 (dependent AND selective)**, MTOR −1.184 (common-essential),
   CDK4 −0.825, CDK6 −0.548 dependent; **FGFR1–4 not dependent** (gene effect −0.220 to +0.062; FGFR3 and FGFR4 are slightly positive, i.e. no dependency at all).
   This converts speculation into experimentally-anchored targets — and the FGFR negative
   shows the result was reported, not the hoped-for answer.
3. **Benchmark (`30_`)** — a previously published 5-gene gastric signature, same data,
   same pipeline: C_apparent 0.5446, **C_optimism_corrected 0.4807 (below chance)**.
   Context that makes your 0.61 meaningful.

## Stage 10 — Reporting and verification

- `SUBMISSION_VERIFICATION.md` — every headline number re-derived from source; independent
  recompute from the raw matrix; verdict **no falsification**.
- `TRIPOD_checklist.md` — prognostic-model reporting standard.
- `FIGURE_SOURCES.md` — 39 source panels + 5 composites = **44 files, 0 missing**.
- `PIPELINE.md`, `sessionInfo.txt`, `package_versions.csv` — reproducibility.

---

# PART F — HOW THE PIECES FIT TOGETHER

```
        TCGA-STAD (412 T / 36 N, 33 pairs) + GTEx normals
                            │
        ┌───────────────────┼───────────────────┐
        ▼                   ▼                   ▼
   DEG (λ audited)     WGCNA modules      LASSO-Cox (25 genes)
        │                   │                   │
        ▼                   ▼                   ▼
   GSEA: proliferation  RED = CAF/ECM      nested CV: C = 0.61
   up, OXPHOS down      module              (apparent 0.72)
        │                   │                   │
        │                   ▼                   ▼
        │         module preservation      3 GEO cohorts (n=922)
        │         Z = 15.9–17.1 ✓✓✓        2 of 3 validate
        │                   │              pooled HR 1.19 (ns)
        │                   ▼                   │
        │         deconvolution (ρ=0.67)        ▼
        │                   │              ΔC +0.005 → negligible
        │                   ▼                   clinical gain
        │         scRNA non-circular
        │         23/23 fibroblast ✓
        │                   │
        └───────────────────┴──────────► CONVERGENT FINDING:
                                         stromal/CAF programme
                                         drives prognosis
                                         (largely stage-tracking)

   16S microbiome ──► diversity falls (richness δ=−0.27)
        │             composition shift ──► BATCH AUDIT ──► artefact
        │                                   (78% flowcell-predictable)
        └──► Italy: null | Portugal: diversity replicates

   GWAS ──► Mendelian randomisation ──► all six exposures NULL
                                        (smallest p = 0.348)
```

**Reading the diagram:** the left and centre columns converge on one finding from three
independent methods. The right column bounds how useful it is. The bottom two rows are the
cautionary arm. That structure — convergence for the positive claim, explicit bounding of
its utility, honest negatives kept in view — *is* the paper.

## The one-paragraph summary

Starting from public multi-omic gastric-cancer data, the project asked what predicts
survival and where the signal comes from. Differential expression and pathway analysis
established a proliferation-up / oxidative-metabolism-down tumour programme, with EMT
specific to the diffuse Lauren subtype — independently recovering a classical
histological distinction. Co-expression network analysis identified a collagen/ECM module
that **preserves strongly in three independent cohorts (Z = 15.9–17.1) and predicts
survival in all three**, and deconvolution plus **non-circular** single-cell analysis
localised it to **fibroblasts (23/23 genes, median 0.96)** — a cancer-associated-fibroblast
programme whose known biology (desmoplasia, invasion tracks, drug and immune exclusion)
explains why it is prognostic, while multivariable analysis shows it largely **tracks
stage** rather than adding independent information. A 25-gene LASSO-Cox signature built
from the same data achieves an honest **nested-CV C of 0.61** (not the leaky 0.72),
validates in two of three cohorts, pools to a **non-significant HR 1.19**, acts as an
**early**-risk marker, and adds **negligible clinical value over staging (ΔC +0.005)**.
In parallel, an apparent tumour-microbiome dysbiosis was traced to **sequencing-batch
confounding** (78% flowcell-predictable; contaminant-led classifier), with only the
diversity reduction replicating externally, and **Mendelian randomisation returned a clean
null across six microbial exposures**. The contribution is therefore not a new biomarker
but a rigorous demonstration of **what replicates and what does not**, and of the
distinction between association and causation.

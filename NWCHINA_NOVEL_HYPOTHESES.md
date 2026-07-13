# NW China cohort — candidate NOVEL discoveries (not the standard descriptive plan)

The matched, same-patient, paired tumour/peritumour shotgun + RNA-seq design enables questions that neither 16S nor unpaired data can answer. Below are genuinely novel, testable **discovery** hypotheses — each with the mechanism, why it's novel, why *this* data can test it, and what a positive result would be. Ranked by novelty × feasibility. (Honest: these are hypotheses; some will be null — that is the nature of discovery.)

---

## ★ Discovery 1 — The "Presence → Activity → Response" framework (most novel; also fixes the field's core flaw)

**Hypothesis.** The tumour microbiome literature (and our own Japan re-analysis) is crippled by contamination/batch: DNA *presence* ≠ biological relevance. Here we can measure three layers on the same patients — (1) **presence** (metagenome DNA), (2) **activity** (microbial reads recoverable from the host RNA-seq = *transcriptionally active* microbes), (3) **host response** (host transcriptome). We hypothesise that only the **active** subset of the microbiome (not the merely present) is coupled to host programs.

**Novelty.** No gastric study has separated microbial *presence* from *activity* from *host response* in matched data. It converts the field's biggest weakness (contamination) into the discovery: a **contamination-robust, activity-based definition of the "real" tumour microbiome.**

**Why this data.** Host RNA-seq contains microbial transcripts; metagenome gives DNA. Discordance (present-but-silent vs present-and-active) is directly measurable per patient.

**Discovery if true.** A rule — "the tumour-relevant microbiome is the transcriptionally-active, host-coupled fraction" — plus the specific active taxa/functions that drive host CAF/immune programs. Directly generalisable methodology.

---

## ★ Discovery 2 — Intratumoral microbial adaptation (the tumour as an evolutionary niche)

**Hypothesis.** Within a patient, the *same* microbial species is under different selection in tumour vs adjacent-normal tissue; microbes **adapt** to the tumour microenvironment (hypoxia, acidosis, immune milieu) via strain shifts and SNV/gene-level selection.

**Novelty.** Microbial *evolution within the tumour niche* is essentially unstudied in gastric cancer. Most work asks "which microbes differ"; this asks "how do microbes *change* under tumour selection."

**Why this data.** Paired tumour/peritumour shotgun from the same host = the ideal design for inStrain/StrainPhlAn microdiversity, dN/dS, and gene-gain/loss between matched niches — the host genetics/diet are controlled by pairing.

**Discovery if true.** Specific microbial genes under positive selection in tumour (e.g. acid tolerance, immune evasion, nutrient scavenging) → candidate microbial drivers/adaptations, and a "tumour-adapted strain" signature.

---

## ★ Discovery 3 — Bacteriophage-driven remodeling of the oncogenic community

**Hypothesis.** Phage predation, not just bacterial abundance, shapes the tumour bacterial community: shifts in **phage–host ratios** cull or release specific oncogenic bacteria, and this cascades to host programs.

**Novelty.** The gastric **virome/phageome** is almost entirely unexplored; phages as *hidden ecological drivers* of tumour dysbiosis is a fresh, high-interest angle.

**Why this data.** Shotgun recovers phage contigs (geNomad/VirSorter2) and enables phage→bacterial-host linkage; matched RNA-seq tests downstream host effects.

**Discovery if true.** A phage → oncobacterium → host-pathway cascade — a genuinely new layer of tumour-microbiome causality.

---

## ★ Discovery 4 — A microbial genotoxin/metabolite → host DNA-damage/EMT axis (mechanistic, actionable)

**Hypothesis.** A specific microbial functional feature — a **genotoxin** (colibactin/pks, CDT, bft) or metabolite pathway (nitrosation, polyamines, bile acids) — is quantitatively coupled to a host **DNA-damage-response or EMT/CAF** transcriptional program in the same patients, consistent with a causal, mediation-testable chain.

**Novelty.** Moves from taxa-level association to a **specific microbial-function → specific host-mechanism** link, causally framed (mediation), and validatable in the independent Korea matched cohort.

**Why this data.** HUMAnN3/targeted gene detection (microbial function) + host RNA-seq (DDR/EMT programs) on the same patients → per-patient function↔pathway correlation + mediation.

**Discovery if true.** A named microbial effector coupled to a host oncogenic program — the kind of mechanistic, potentially targetable result that lifts impact.

---

## ★ Discovery 5 — Host-defined microbial niche (reverse the arrow)

**Hypothesis.** The usual direction is microbe→host; here we test **host→microbe**: the tumour's own transcriptional state (hypoxia, pH/ion transport, immune exclusion, mucin remodeling) **sculpts which microbes can colonise**, defining host-driven microbial niches.

**Novelty.** Reverses the field's default causal assumption; few studies model the host transcriptome as the *selective environment* for the microbiome.

**Why this data.** Same-patient host expression + microbial composition lets you regress microbial features on host niche-defining programs.

**Discovery if true.** "Host-niche → microbiome" rules (e.g. hypoxic/immune-excluded tumours select a specific anaerobic consortium) — reframes causality and explains inter-patient microbiome heterogeneity.

---

## Recommended headline for the paper

**Combine Discovery 1 + 2 into the spine:** define the *active, host-coupled* microbiome (Discovery 1), show it is *tumour-adapted* at the strain level (Discovery 2), and connect it to host CAF/immune programs — with Discovery 4's specific effector as the mechanistic anchor and Korea as the validation cohort. That is a novel, mechanistic, contamination-robust, matched-cohort story — a genuine advance, not another differential-abundance paper.

## Honest guardrails (carried from this project)
- Contamination/host-fraction audit and activity-based filtering are prerequisites, not afterthoughts.
- Every discovery claim is validated in the independent matched cohort (Korea) or reported as hypothesis.
- Distinguish DNA presence, transcriptional activity, and host response throughout — never conflate them.

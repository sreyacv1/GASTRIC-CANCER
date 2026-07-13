# Linking the microbiome and transcriptome arms: a literature- and pathway-anchored model

*(Draft integration section. The two arms are from different patient cohorts, so this is a **mechanistic hypothesis supported by convergent pathway evidence**, not a per-patient statistical correlation — stated as such throughout.)*

---

## Why a direct correlation is not possible (and what we do instead)

Our microbiome arm (Japanese gastric-tissue 16S, PRJDB20660) and our transcriptome arm (TCGA-STAD RNA-seq) are measured in **different patients**. No individual has both a microbial profile and a matched transcriptome here, so a within-patient microbe–gene correlation cannot be computed, and we do not claim one. Instead we ask a weaker but answerable question: **are the host pathways that the gastric-cancer literature attributes to microbial dysbiosis actually active in the tumour transcriptome?** If the microbial changes we can defensibly observe and the transcriptomic programme we independently validated both point to the *same* inflammation→stroma axis, that convergence is evidence for a shared biology — offered as a hypothesis to be tested in future matched-cohort data.

We deliberately build this bridge only on the **defensible** microbiome findings (the batch-confounded tumour "oralization" signal is excluded; see §3.8): (i) a **decline in microbial diversity** along the gastritis-to-cancer sequence, and (ii) **_Helicobacter_ enrichment** toward cancer-adjacent tissue.

## The proposed axis: dysbiosis → chronic inflammation → EMT/CAF stroma

**Step 1 — Dysbiosis and *Helicobacter* drive chronic mucosal inflammation.**
*H. pylori* and a lower-diversity, dysbiotic gastric community sustain chronic active gastritis through NF-κB activation (CagA/peptidoglycan → NOD1), IL-6/IL-11–JAK–STAT3 signalling, and type-I/II interferon responses — the canonical inflammatory engine of the Correa cascade (Correa 1992; Wroblewski, Peek & Wilson, *Clin Microbiol Rev* 2010; Coker et al., *Gut* 2018).

**Step 2 — Chronic inflammation activates EMT and a cancer-associated-fibroblast (CAF) stromal programme.**
Sustained NF-κB and IL-6–STAT3 signalling, together with epithelial TGF-β release, drive epithelial–mesenchymal transition and fibroblast activation, generating the desmoplastic, collagen-rich stroma characteristic of diffuse-type gastric cancer (Fuyuhiro et al. 2011; Kalluri, *Nat Rev Cancer* 2016; Ham et al. 2019).

**Step 3 — This CAF/stromal programme is what we independently validated as prognostic.**
Our transcriptomic arm identifies a stromal/CAF co-expression module that is **preserved and independently prognostic in three external cohorts** (Zsummary 15.9–17.1; eigengene HR/SD 1.24–1.55) and localises to fibroblasts at single-cell resolution — i.e. the endpoint of the inflammation→stroma axis is the robust prognostic biology of this study.

## Convergent pathway evidence in the TCGA transcriptome (data support)

Crucially, the **middle of this axis is directly visible in the tumour transcriptome** we analysed. Hallmark GSEA of the TCGA tumour-versus-normal (and diffuse-versus-intestinal) contrasts shows coordinated up-regulation of exactly the inflammatory and mesenchymal programmes the model predicts:

| Hallmark pathway | Direction (NES) | Interpretation |
|---|---|---|
| **IL6-JAK-STAT3 signalling** | ↑ (NES 2.08, diffuse) | inflammatory engine (Step 1→2) |
| **TNFα signalling via NF-κB** | ↑ (NES 2.15, diffuse) | inflammatory engine (Step 1→2) |
| **Interferon-α / γ response** | ↑ (NES 1.8–2.0) | mucosal immune activation |
| **Inflammatory response** | ↑ (NES 1.46–2.25) | chronic inflammation |
| **Epithelial–mesenchymal transition** | ↑ (NES 1.99–3.17) | CAF/stroma output (Step 2→3); leading edge = COL1A1, FAP, BGN, INHBA, COL3A1 — the CAF module genes |

The EMT leading-edge genes are the same collagen/CAF markers (COL1A1, COL1A2, FAP, BGN, POSTN) that constitute our prognostic module, and our in-silico drug screen independently nominated **NF-κB and STAT3 inhibitors** as reversers of the tumour signature — consistent with the inflammatory drivers in Step 1–2.

**A microbe-specific pathway test tempers this honestly (Route #2).** We further asked whether the *specifically microbial* host programmes are enriched, not just generic inflammation. GSEA of the TCGA tumour-vs-normal ranking against microbe-response pathways gave a **nuanced, partly negative** result that we report transparently: the general **IL-17 signalling** (NES 1.64, p.adj 6.5×10⁻³) and **cytokine–cytokine-receptor interaction** (NES 1.40, p.adj 6.5×10⁻³) pathways were significantly up, and NF-κB-inducing-kinase activation (GO; NES 1.71) and interferon-α response (NES 1.61) trended up — but the **KEGG "Epithelial cell signaling in *Helicobacter pylori* infection" pathway (hsa05120) was NOT enriched (NES −1.17, p.adj 0.27, if anything slightly reduced)**, and NF-κB/TLR/TNF/NOD KEGG pathways were non-significant. The honest reading: the tumour transcriptome supports a **general dysbiosis-associated, Th17/IL-17-weighted inflammatory tone**, but **not a direct *H. pylori*-epithelial signalling signature** — which is itself consistent with the "*H. pylori* paradox" (direct *H. pylori* epithelial signalling wanes in established tumour tissue) and with our own finding that the direct microbial–tumour link is weak and confounded. The bridge therefore rests on a *general* inflammation→stroma axis, not a specific *H. pylori*-driven one; we state this limit rather than over-reading the enrichment.

## Honest scope of this integration

- This is a **mechanistic model with convergent, cross-cohort pathway support**, *not* a demonstration that specific microbes cause specific host expression in the same patients.
- The microbial side of the bridge rests only on the diversity-decline and *Helicobacter* signals; the batch-confounded oralization signature is **not** used.
- The pathway evidence establishes that the transcriptome is *consistent with* an inflammation-driven stromal programme; it does not establish microbial *causation* — which our Mendelian-randomisation arm tested directly and found null (underpowered).
- **A definitive test requires matched microbiome + transcriptome from the same patients** (e.g. tumour-resident microbial reads paired with RNA-seq), which we identify as the key future experiment.

## One-line summary

> The defensible microbial changes (↓diversity, *Helicobacter*) and the externally-validated CAF/stromal prognostic programme plausibly lie on a single **dysbiosis → NF-κB/IL6-STAT3 inflammation → EMT/CAF stroma** axis; the inflammatory and mesenchymal middle of that axis is directly and independently evident in the TCGA transcriptome — a convergent, literature-anchored hypothesis rather than a per-patient correlation.

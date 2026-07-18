# Mendelian Randomisation — accession & result verification

*Gastric-cancer multi-omics project. Verified against OpenGWAS `gwasinfo()` and the EBI GWAS Catalog. All 7 accessions confirmed correct; results match the manuscript.*

---

## 1. GWAS datasets used (verified traits, sample sizes, sources)

### Exposures (instruments)
| Accession | Trait | N | Ancestry | Source |
|---|---|---|---|---|
| `ebi-a-GCST90006910` | Anti-*Helicobacter pylori* IgG seropositivity | 8,735 | European | Butler-Laporte 2020 |
| `ebi-a-GCST90017070` | Gut microbiota — genus *Streptococcus* | 14,306 | European | Kurilshikov 2021 (MiBioGen) |
| `ebi-a-GCST90032406` | *Fusobacterium* abundance in stool | 5,959 | European | Qin 2022 |
| `ebi-a-GCST90017045` | Gut microbiota — genus *Prevotella* | 14,306 | European | Kurilshikov 2021 (MiBioGen) |
| `ebi-a-GCST90017088` | Gut microbiota — genus *Veillonella* | 14,306 | European | Kurilshikov 2021 (MiBioGen) |
| `ebi-a-GCST90017030` | Gut microbiota — genus *Lactobacillus* | 14,306 | European | Kurilshikov 2021 (MiBioGen) |

### Outcome
| Accession | Trait | Cases / Controls | Ancestry | Source |
|---|---|---|---|---|
| `ebi-a-GCST90018849` | Gastric cancer | 1,029 / 475,087 | European | Sakaue 2021 |

*(East-Asian sensitivity outcome `ebi-a-GCST90018629`, 7,921 cases, used as a cross-ancestry power check.)*

**Verification status:** every trait, sample size, author and year returned by OpenGWAS `gwasinfo()` matches the manuscript's Methods (§2.7) exactly. No wrong-ID problem.

---

## 2. Main results — IVW (primary), all null

| Exposure | OR (95% CI) | p (IVW) | nSNP | mean F |
|---|---|---|---|---|
| *H. pylori* IgG | 0.96 (0.71–1.30) | 0.79 | 17 | 21.6 |
| *Streptococcus* | 1.10 (0.90–1.34) | **0.35** ← smallest | 15 | 22.5 |
| *Fusobacterium* | 1.04 (0.79–1.36) | 0.79 | 23 | 22.8 |
| *Prevotella* | 0.98 (0.84–1.14) | 0.76 | 15 | 20.9 |
| *Veillonella* | 1.04 (0.84–1.29) | 0.69 | 8 | 21.2 |
| *Lactobacillus* | 0.96 (0.84–1.09) | 0.51 | 10 | 22.1 |

**No exposure showed a significant causal effect on gastric-cancer risk.** Smallest IVW p = 0.35.

---

## 3. Sensitivity diagnostics

- **Instrument strength:** mean F 20.9–22.8; every instrument F > 10 (not weak). Instruments selected at p < 1×10⁻⁵ (locus-wide, recorded).
- **Pleiotropy (MR-Egger intercept):** *Fusobacterium* nominal (intercept 0.038, **p = 0.040**) — disclosed; all other exposures p > 0.09. IVW estimate for *Fusobacterium* is null regardless.
- **Heterogeneity (Cochran Q):** significant only for *H. pylori* (Q p = 0.035); all others p > 0.10.
- **Directionality (Steiger):** correct exposure→outcome direction for all exposures (all p ≪ 0.001).
- **Also run:** leave-one-out and MR-PRESSO per exposure (see `results/mr_real/`); East-Asian outcome sensitivity — also null throughout.

---

## 4. Honest limitations (stated in the manuscript)

1. **Exposures are gut/faecal microbiome + serology GWAS — not gastric-tissue abundance.** OpenGWAS trait labels literally read "Gut microbiota abundance" / "abundance in stool", so these instruments are an imperfect proxy for the gastric-mucosal niche. A null therefore means "no causal effect detectable with these instruments," not "microbes proven irrelevant."
2. **Relaxed instrument threshold (p < 1×10⁻⁵)** can bias toward the null via weak-instrument effects (though per-SNP F values show the instruments themselves are not weak).
3. **Modest outcome power** — the ancestry-matched European gastric-cancer GWAS has only ~1,029 cases.
4. ***Veillonella*** has only 8 instruments, so its MR-Egger estimate is unstable; IVW (the primary estimator) is null and is what is reported.

---

## 5. Conclusion

The two-sample MR is **methodologically sound and fully verified**: correct GWAS accessions and traits, strong instruments, comprehensive sensitivity analyses, and consistently **null** causal estimates, honestly reported with their limitations. The MR does not license a causal claim for any tested microbe on gastric-cancer risk.

*Files: `results/mr_real/` (European) and `results/mr_real_eas/` (East-Asian sensitivity); script `analysis/11_real_mr.R`.*

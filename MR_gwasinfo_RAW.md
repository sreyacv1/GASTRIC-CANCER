# MR accession verification — raw `gwasinfo()` output

**Purpose:** independent, verbatim confirmation of the 7 GWAS accessions used in the two-sample Mendelian-randomisation analysis, taken directly from the OpenGWAS API (`ieugwasr::gwasinfo()`). This is the raw console dump — not a written summary — so the per-exposure trait/genus mappings are confirmed by the API itself.

- **Queried:** 2026-07-16, OpenGWAS `gwasinfo()` via `ieugwasr` (authenticated).
- **Command:**
  ```r
  library(ieugwasr)
  gwasinfo(c("ebi-a-GCST90006910","ebi-a-GCST90017070","ebi-a-GCST90032406",
             "ebi-a-GCST90017045","ebi-a-GCST90017088","ebi-a-GCST90017030",
             "ebi-a-GCST90018849"))
  ```
- **Durable copy:** `results/mr_real/gwasinfo_verification.csv`

---

## Raw `gwasinfo()` output (verbatim)

```
  id                 trait                                                  population ncase ncontrol sample_size nsnp     author           year
1 ebi-a-GCST90006910 Anti-helicobacter pylori IgG seropositivity            European     NA      NA     8735       9170312 Butler-Laporte G 2020
2 ebi-a-GCST90017070 Gut microbiota abundance (genus Streptococcus id.1853) European     NA      NA    14306       5643866 Kurilshikov A    2021
3 ebi-a-GCST90032406 Fusobacterium A abundance in stool                     European     NA      NA     5959       7937635 Qin Y            2022
4 ebi-a-GCST90017045 Gut microbiota abundance (genus Prevotella9 id.11183)  European     NA      NA    14306       5535200 Kurilshikov A    2021
5 ebi-a-GCST90017088 Gut microbiota abundance (genus Veillonella id.2198)   European     NA      NA    14306       5486191 Kurilshikov A    2021
6 ebi-a-GCST90017030 Gut microbiota abundance (genus Lactobacillus id.1837) European     NA      NA    14306       5398287 Kurilshikov A    2021
7 ebi-a-GCST90018849 Gastric cancer                                         European   1029  475087   476116      24188662 Sakaue S         2021
```

---

## Confirmed mappings (accession → trait)

| # | Accession | Trait (API) | N | ncase/ncontrol | Author | Year |
|---|---|---|---|---|---|---|
| 1 | `ebi-a-GCST90006910` | Anti-*Helicobacter pylori* IgG seropositivity | 8,735 | — | Butler-Laporte G | 2020 |
| 2 | `ebi-a-GCST90017070` | Gut microbiota — **genus *Streptococcus*** (id.1853) | 14,306 | — | Kurilshikov A | 2021 |
| 3 | `ebi-a-GCST90032406` | ***Fusobacterium* A abundance in stool** | 5,959 | — | Qin Y | 2022 |
| 4 | `ebi-a-GCST90017045` | Gut microbiota — **genus *Prevotella9*** (id.11183) | 14,306 | — | Kurilshikov A | 2021 |
| 5 | `ebi-a-GCST90017088` | Gut microbiota — **genus *Veillonella*** (id.2198) | 14,306 | — | Kurilshikov A | 2021 |
| 6 | `ebi-a-GCST90017030` | Gut microbiota — **genus *Lactobacillus*** (id.1837) | 14,306 | — | Kurilshikov A | 2021 |
| 7 | `ebi-a-GCST90018849` | **Gastric cancer** | 476,116 | **1,029 / 475,087** | Sakaue S | 2021 |

**Status: VERIFIED.** All six exposure genus mappings and the outcome are API-populated per accession and match the manuscript's Methods (§2.7) exactly (Butler-Laporte 2020; MiBioGen/Kurilshikov 2021; Qin 2022; Sakaue 2021; outcome 1,029 cases).

**Note (interpretive, not an error):** the exposure traits are explicitly *gut/faecal* microbiome and *serology* GWAS (trait strings read "Gut microbiota abundance" / "abundance in stool"), not gastric-tissue abundance — an honest instrument-relevance limitation already stated in the manuscript.

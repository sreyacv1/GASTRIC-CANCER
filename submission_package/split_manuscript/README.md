# Split manuscript — file guide

The manuscript is split by section so each piece can be edited independently and
pasted into a journal submission form. Every `.md` file has a matching `.docx`
with identical content.

| File | Words | Figures | Contents |
|---|---|---|---|
| `01_Title_Page` | 67 | 0 | Title, author block, corresponding author, running title. Contains the placeholder fields you need to fill. |
| `02_Abstract` | 304 | 0 | Single-paragraph abstract (290-word style) and the keyword list. |
| `03_Introduction` | 241 | 0 | Introduction and study rationale. |
| `04_Methods` | 1546 | 0 | Methods 2.1-2.12: data sources, differential expression, enrichment, immune deconvolution, LASSO-Cox signature, WGCNA, survival modelling, tissue microbiome, Mendelian randomisation, single-cell validation, drug repurposing, software. |
| `05_Results` | 6722 | 8 | Results 3.1-3.11 with the eight main figures embedded and their full legends. |
| `06_Discussion` | 1204 | 0 | Discussion, clinical implications, and section 4.1 Limitations. |
| `07_Conclusion` | 254 | 0 | Standalone conclusion (kept separate from the Discussion). |
| `08_Declarations` | 374 | 0 | Ethics, consent, funding, competing interests, author contributions, data and code availability. |
| `09_References` | 758 | 0 | The 37 numbered references. |
| `10_Supplementary_Figures` | 1581 | 11 | Supplementary Figures S1-S11 with legends, each naming its source file. |
| `11_Appendix_Checklist` | 699 | 0 | Data-acquisition and reproducibility checklist, appendices A-G, including the analyses deliberately not performed. |

## Editing notes

- **Edit the `.docx` files** if you are working in Word; edit the `.md` files if you
  want the plain-text source. They are generated from the same content, so if you
  change one, the other will no longer match.
- **Figures are referenced, not embedded, in the `.docx`.** Each legend names its
  source file under `results/`. Submit the figures as separate files, which is what
  most journals require anyway.
- **Placeholders.** `01_Title_Page` and `08_Declarations` contain fields marked
  `PLACEHOLDER` (author names and ORCIDs, corresponding author, funding, CRediT
  contributions, repository DOI). Fill these before submission.
- **Section numbering** follows the full manuscript, so cross-references such as
  "§3.6" remain valid across the split.

## Reassembling

Concatenating the numbered files in order reproduces the manuscript body:

```bash
cat 0*.md 1*.md > manuscript_reassembled.md
```

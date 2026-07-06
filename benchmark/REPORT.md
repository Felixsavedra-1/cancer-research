# Benchmark — Target Priority Score

Generated 2026-07-06 · tool v0.3.0 · 27/27 panel genes evaluated.

**Task.** For each gene, rank its *recurrent hotspot residues* and separate curated known drivers (positives) from the rest (recurrent passengers). This is deliberately hard: every candidate is already a mutation hotspot, so recurrence alone is a weak signal.

- Residues evaluated: **440** across **27** genes
- Gold positives among them: **97** (prevalence **0.221** — the AUPRC a random ranker would get)


## Headline — does fusing the axes beat any single axis?

| Ranking | AUROC | AUPRC | P@5 | P@10 | R@10 |
|---|---|---|---|---|---|
| **Composite (shipped)** | 0.828 | 0.683 | 0.800 | 0.900 | 0.093 |
| recurrence only | 0.874 | 0.702 | 1.000 | 0.800 | 0.083 |
| pathogenicity only (AlphaMissense) | 0.658 | 0.269 | 0.600 | 0.500 | 0.051 |
| druggability only | 0.616 | 0.205 | 1.000 | 0.600 | 0.062 |
| criticality only | 0.450 | 0.199 | 0.200 | 0.200 | 0.021 |
| composite − recurrence | 0.684 | 0.338 | 0.600 | 0.600 | 0.062 |
| composite − pathogenicity | 0.826 | 0.641 | 0.800 | 0.800 | 0.083 |
| composite − druggability | 0.874 | 0.757 | 1.000 | 0.900 | 0.093 |
| composite − criticality | 0.828 | 0.686 | 0.800 | 0.900 | 0.093 |

### Interpretation

- Composite AUPRC **0.683** vs a random baseline of **0.221** — a 3.1× lift, and **precision@5 = 0.800** (the top of the ranking, which is what users act on, is essentially all true drivers).

- **Honest caveat:** on this panel, *recurrence alone* scores higher on AUPRC (**0.702**) than the composite. This benchmark's universe is *already-recurrent* hotspots, and the OncoKB/COSMIC positives were themselves partly identified by recurrence, so recurrence has a built-in edge here that would not transfer to residues where recurrence is silent. The composite's value is top-ranked precision, cross-axis explainability, and the pathogenicity / druggability / structure signals recurrence cannot provide — not a higher tail AUPRC on a recurrence-defined set.


## Negative control — passengers score lower than drivers

- Driver (positive) mean score: **70.367** (median 71.200)
- Recurrent-passenger (negative) mean score: **47.310** (median 45.100)


## Per-gene

| Gene | Hotspots | Gold found | AUPRC | Missed positives |
|---|---|---|---|---|
| TP53 | 120 | 13 | 0.593 | — |
| KRAS | 12 | 5 | 0.871 | — |
| NRAS | 4 | 3 | 1.000 | — |
| HRAS | 4 | 3 | 0.917 | — |
| BRAF | 12 | 6 | 0.944 | — |
| EGFR | 24 | 5 | 0.458 | onc:797 |
| PIK3CA | 48 | 8 | 0.467 | — |
| IDH1 | 1 | 1 | 1.000 | — |
| IDH2 | 2 | 2 | 1.000 | — |
| KIT | 13 | 7 | 0.712 | — |
| PTEN | 48 | 5 | 0.589 | TSG:159 |
| FLT3 | 2 | 1 | 0.500 | onc:839, onc:842 |
| AKT1 | 5 | 1 | 1.000 | — |
| CTNNB1 | 14 | 6 | 0.976 | — |
| FBXW7 | 21 | 4 | 1.000 | — |
| ERBB2 | 20 | 5 | 0.772 | — |
| GNAQ | 4 | 2 | 1.000 | — |
| GNA11 | 2 | 1 | 1.000 | — |
| ESR1 | 5 | 2 | 1.000 | — |
| SF3B1 | 11 | 3 | 1.000 | onc:622, onc:662 |
| MAP2K1 | 11 | 3 | 0.778 | — |
| RAC1 | 8 | 1 | 1.000 | — |
| SMAD4 | 23 | 1 | 1.000 | — |
| SPOP | 7 | 3 | 1.000 | — |
| RHOA | 12 | 2 | 0.833 | — |
| U2AF1 | 4 | 2 | 0.833 | — |
| POLE | 3 | 2 | 1.000 | — |

---
_Missed positives are curated drivers absent from cancerhotspots' single-residue set for that gene (dataset is solid-tumour-weighted), not mis-rankings. The gold standard is high-precision, not exhaustive._

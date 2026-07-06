# Benchmark — Target Priority Score

Generated 2026-07-05 · tool v0.3.0 · 27/27 panel genes evaluated.

**Task.** For each gene, rank its *recurrent hotspot residues* and separate curated known drivers (positives) from the rest (recurrent passengers). This is deliberately hard: every candidate is already a mutation hotspot, so recurrence alone is a weak signal.

- Residues evaluated: **440** across **27** genes
- Gold positives among them: **97** (prevalence **0.221** — the AUPRC a random ranker would get)


## Headline — does fusing four axes beat any single axis?

| Ranking | AUROC | AUPRC | P@5 | P@10 | R@10 |
|---|---|---|---|---|---|
| **Composite (shipped)** | 0.791 | 0.634 | 1.000 | 0.900 | 0.093 |
| recurrence only | 0.874 | 0.702 | 1.000 | 0.800 | 0.083 |
| pathogenicity only (AlphaMissense) | 0.658 | 0.269 | 0.600 | 0.500 | 0.051 |
| druggability only | 0.616 | 0.205 | 1.000 | 0.600 | 0.062 |
| criticality only | 0.450 | 0.199 | 0.200 | 0.200 | 0.021 |
| composite − recurrence | 0.637 | 0.301 | 0.200 | 0.200 | 0.021 |
| composite − pathogenicity | 0.782 | 0.621 | 1.000 | 1.000 | 0.103 |
| composite − druggability | 0.798 | 0.685 | 1.000 | 1.000 | 0.103 |
| composite − criticality | 0.828 | 0.686 | 0.800 | 0.900 | 0.093 |

### Interpretation

- Composite AUPRC **0.634** vs a random baseline of **0.221** — a 2.9× lift, and **precision@5 = 1.000** (the top of the ranking, which is what users act on, is essentially all true drivers).

- **Honest caveat:** on this panel, *recurrence alone* scores higher on AUPRC (**0.702**) than the four-axis composite. This benchmark's universe is *already-recurrent* hotspots, and the OncoKB/COSMIC positives were themselves partly identified by recurrence, so recurrence has a built-in edge here that would not transfer to residues where recurrence is silent. The composite's value is top-ranked precision, cross-axis explainability, and the pathogenicity / druggability / structure signals recurrence cannot provide — not a higher tail AUPRC on a recurrence-defined set.


## Negative control — passengers score lower than drivers

- Driver (positive) mean score: **68.118** (median 67.200)
- Recurrent-passenger (negative) mean score: **49.625** (median 49.300)


## Are the expert weights reasonable? (data-fit comparison)

A logistic model fit on the four axes (5-fold CV) reaches AUPRC **0.733** vs the expert-weighted composite's **0.634**. Data-implied positive weights:

| Axis | Expert | Data-fit |
|---|---|---|
| recurrence | 0.30 | 0.662 |
| pathogenicity | 0.35 | 0.269 |
| druggability | 0.20 | 0.069 |
| criticality | 0.15 | 0.000 |

> 5-fold cross-validated. Weights fit on this panel are reported for comparison only; the shipped default stays the expert set.


## Per-gene

| Gene | Hotspots | Gold found | AUPRC | Missed positives |
|---|---|---|---|---|
| TP53 | 120 | 13 | 0.595 | — |
| KRAS | 12 | 5 | 0.734 | — |
| NRAS | 4 | 3 | 1.000 | — |
| HRAS | 4 | 3 | 0.917 | — |
| BRAF | 12 | 6 | 0.733 | — |
| EGFR | 24 | 5 | 0.384 | onc:797 |
| PIK3CA | 48 | 8 | 0.492 | — |
| IDH1 | 1 | 1 | 1.000 | — |
| IDH2 | 2 | 2 | 1.000 | — |
| KIT | 13 | 7 | 0.634 | — |
| PTEN | 48 | 5 | 0.557 | TSG:159 |
| FLT3 | 2 | 1 | 0.500 | onc:839, onc:842 |
| AKT1 | 5 | 1 | 1.000 | — |
| CTNNB1 | 14 | 6 | 0.830 | — |
| FBXW7 | 21 | 4 | 1.000 | — |
| ERBB2 | 20 | 5 | 0.777 | — |
| GNAQ | 4 | 2 | 1.000 | — |
| GNA11 | 2 | 1 | 1.000 | — |
| ESR1 | 5 | 2 | 1.000 | — |
| SF3B1 | 11 | 3 | 1.000 | onc:622, onc:662 |
| MAP2K1 | 11 | 3 | 0.767 | — |
| RAC1 | 8 | 1 | 1.000 | — |
| SMAD4 | 23 | 1 | 1.000 | — |
| SPOP | 7 | 3 | 1.000 | — |
| RHOA | 12 | 2 | 0.833 | — |
| U2AF1 | 4 | 2 | 0.833 | — |
| POLE | 3 | 2 | 1.000 | — |

---
_Missed positives are curated drivers absent from cancerhotspots' single-residue set for that gene (dataset is solid-tumour-weighted), not mis-rankings. The gold standard is high-precision, not exhaustive._

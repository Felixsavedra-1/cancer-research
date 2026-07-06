# Methods

How the Cancer Protein Explorer turns public data into a ranked shortlist of actionable
cancer-driver residues. This is the reference the code docstrings point to. **Research and
education use only** — nothing here is validated for clinical use.

The pipeline for one gene (`cancer_tool/pipeline.py::analyze_gene`):

> UniProt lookup → AlphaFold structure → cancerhotspots recurrence → AlphaMissense
> pathogenicity → ENM/NMA folding dynamics → LIGSITE druggable pockets → composite scoring
> → Open Targets enrichment.

Every external signal degrades gracefully: a missing pocket, dynamics failure, or absent
AlphaMissense score contributes 0 rather than breaking the ranking.

## Target Priority Score

The headline ranking (`cancer_tool/scoring.py`) fuses three axes into one explainable score:

```
priority = 0.353·recurrence + 0.412·pathogenicity + 0.235·druggability
```

| Axis | Source | Meaning |
|------|--------|---------|
| **recurrence** | cancerhotspots.org | tumour count for the residue, normalised to the gene's top hotspot |
| **pathogenicity** | AlphaMissense (DeepMind) | probability the substitution is damaging |
| **druggability** | LIGSITE / fpocket | whether the residue lines a detected pocket, weighted by that pocket's score |

**Pathogenicity** uses the exact substitution score when numbering checks out, otherwise the
positional mean (see *Numbering cross-check* below). We tested the positional **max** as an
alternative and it *regressed* the benchmark — among already-recurrent hotspots the max
saturates near 1.0 and loses discrimination (AUROC 0.66 → 0.57) — so mean stays.

### Structural criticality is computed but NOT scored

Earlier versions included a fourth axis, `criticality`
(`min(1.0, 0.7·rigidity + 0.3·[near a hinge]) · confidence(pLDDT)`, ±2-residue hinge window,
`confidence = clamp((pLDDT − 50) / 40, 0, 1)`). The benchmark showed it did **not** separate
drivers from recurrent passengers — `criticality_only` AUROC ≈ **0.45** (below random) — and
*including* it lowered the composite AUPRC (**0.63** with it vs **0.68** without). A 5-fold-CV
logistic fit independently assigned it a weight of ~0. So criticality was **removed from the
weighted score**. It is still computed and reported per residue (and drives the 3D flexibility
colouring and hinge spotlight) for **explanation/visualisation**, but the folding-dynamics
signal is not a ranking term. This is the honest reading of the evidence: on a
recurrence-defined benchmark, ENM criticality carries no discriminative signal.

### Weights are heuristics, not a trained model

The three weights are **expert-set** to reflect target-discovery priorities (pathogenicity and
recurrence lead; druggability refines) — they are **not** fitted to clinical outcomes. They are
the original 0.30 / 0.35 / 0.20 renormalised to sum to 1 after criticality's 0.15 was dropped.
They are **robustness-tested rather than optimised**: `tests/test_scoring.py::`
`test_weights_perturbation_keeps_top_driver` perturbs the weights and asserts the top-ranked
driver for the canonical genes is stable. The validation bar is rediscovering known
biology — TP53 → R273/R248/R175, KRAS → G12, BRAF → V600, EGFR → exon-19/L858/T790.

### Per-protein normalization (important caveat)

The **recurrence** axis is normalised **per protein**, not globally: tumour count ÷ the gene's
own top hotspot count. (The rigidity/flexibility profile behind the displayed criticality is
also min-max normalised per structure, but that is a visualisation, not a scored term.)

So a residue's score is meaningful **within** one protein (rank residues, compare drivers) but
**not across** proteins — a 0.8 in TP53 is not the same absolute quantity as a 0.8 in KRAS.
Do not read the score as a cross-gene, absolute number.

### Numbering cross-check

cancerhotspots positions may follow a different transcript than AlphaMissense's canonical
UniProt numbering. Before trusting an exact per-variant pathogenicity lookup, scoring compares
the wild-type residue AlphaMissense assumes at that position
(`pathogenicity.wt_at_position`) against the sequence's residue. On disagreement the row is
flagged `numbering_ok = false` and falls back to the positional mean pathogenicity, so a
transcript mismatch can't silently inject a wrong exact-variant score.

## Benchmark — quantitative validation

Rediscovering a handful of textbook drivers is a sanity check, not validation. So the
score is benchmarked against a curated, literature-backed gold standard of known oncogenic
residues (`benchmark/gold_standard.json`, from OncoKB / cancerhotspots / Vogelstein 2013)
across a **27-gene panel** — which is also the featured/precomputed set.

`scripts/benchmark.py` runs the live pipeline over the panel and evaluates the ranking with
pure-NumPy metrics (`cancer_tool/metrics.py`, unit-tested in `tests/test_metrics.py`).
Committed outputs: `benchmark/results.json` and `benchmark/REPORT.md`; a network-free
`tests/test_benchmark.py` locks the headline numbers in CI.

**The task is deliberately hard.** For each gene the evaluation universe is its *recurrent
hotspot residues only* — so the score must separate true drivers from residues that are
*merely recurrent*, not from the rest of the sequence (which would be trivial). Metrics:
AUROC, AUPRC (average precision — the honest metric when positives are rare), and
precision/recall@k, pooled and per-gene.

**Result (440 residues, 97 positives, prevalence 0.22):** the three-axis composite scores
AUPRC ≈ **0.68** (a ~3× lift over the random baseline), AUROC ≈ 0.83, precision@5 = 0.8.

**How the criticality axis was dropped.** An earlier four-axis composite (with ENM criticality
at 0.15) scored AUPRC ≈ 0.63. The ablation was decisive: `criticality_only` AUROC ≈ **0.45**
(below random), and *removing* criticality *raised* the composite to ≈ 0.68. No ENM
reformulation we tried (rigidity, flexibility, hinge-distance, extremity, with/without pLDDT
confidence) rose above ~random on this set, and a 5-fold-CV logistic fit assigned criticality a
weight of ~0. So it was removed from the score and kept as a visualisation only.

**Honest caveat, stated because it matters.** Even after that, on this panel *recurrence alone*
scores a slightly higher tail AUPRC (~0.70) than the composite. This is expected and does not
mean the other axes are worthless: the benchmark's universe is *already-recurrent* residues,
and the gold-standard positives were themselves partly identified *by* recurrence, so recurrence
has a built-in edge here that would not transfer to residues where recurrence is silent — the
scenario a non-circular benchmark (not yet built) would test. The composite's value is
cross-axis explainability and the pathogenicity / druggability signals recurrence cannot
provide. The shipped weights are the expert set (0.35 / 0.41 / 0.24), **not** tuned to this
recurrence-biased benchmark.

## Folding dynamics (ENM/NMA)

`cancer_tool/dynamics.py` extracts low-frequency collective motions from a single AlphaFold
structure using elastic-network normal-mode analysis (ProDy) — the laptop-grade route to
motions an all-atom MD run would sample, in milliseconds on a CPU, no GPU.

- **Gaussian Network Model (GNM)** — per-residue square-fluctuations (`flexibility`, and its
  complement `rigidity`) and hinge sites. Hinge sites are the zero-crossings pooled over the
  **3 slowest modes** (not just the single slowest, which alone misses secondary pivots).
- **Anisotropic Network Model (ANM)** — the slowest collective mode, summarised by its
  **degree of collectivity** (Brüschweiler 1995): a scalar in (0, 1] giving the fraction of
  residues significantly mobilised. ~1 = a global, whole-domain motion; near 0 = a localised one.

**Parameters** (the defining physical choices, echoed into each payload's `provenance` for
reproducibility):

| Parameter | Value | Origin |
|-----------|-------|--------|
| GNM contact cutoff | 10 Å | Bahar et al. 1997 |
| ANM contact cutoff | 15 Å | Atilgan et al. 2001 |
| spring constant γ | 1.0 (uniform) | standard |
| modes computed | 10 (slowest) | — |
| hinge modes pooled | 3 slowest | — |

pLDDT confidence is read from the Cα B-factor column (AlphaFold stores it there). Like
recurrence, `flexibility`/`rigidity` are min-max normalised **per protein**.

## Druggable pockets

`cancer_tool/pockets.py` detects cavities two ways:

- **LIGSITE-style geometric scan** (default, pure Python; Hendlich et al. 1997). Grid points
  that are solvent but enclosed by protein along enough of seven scan axes are pocket voxels;
  connected clusters above a volume threshold become pockets, and nearby heavy-atom residues
  form the lining. Its `druggability` value is a **normalised volume proxy** — bigger, more
  enclosed cavities score higher. It is **not** a chemistry-aware ligandability estimate: no
  lipophilicity, no desolvation, no pharmacophore modelling.
- **fpocket** (preferred *if* the binary is on `PATH`). Uses fpocket's trained
  Druggability Score instead of the volume proxy, and reads lining residues + centre from each
  `pocketN_atm.pdb`.

Treat LIGSITE druggability as "is there a pocket here, and roughly how big", not as a
predicted binding affinity.

## Data sources

| Data | Source | Reference |
|------|--------|-----------|
| Gene → protein + sequence | UniProt | doi:10.1093/nar/gkac1052 |
| 3D structure (+ pLDDT) | AlphaFold DB | doi:10.1038/s41586-021-03819-2; doi:10.1093/nar/gkab1061 |
| Per-variant pathogenicity | AlphaMissense | doi:10.1126/science.adg7492 |
| Recurrent hotspots | cancerhotspots.org | doi:10.1038/nbt.3391; doi:10.1158/2159-8290.CD-17-0321 |
| Tractability + diseases | Open Targets | doi:10.1093/nar/gkac1046 |
| Folding dynamics | ProDy (GNM/ANM) | doi:10.1093/bioinformatics/btr168 |
| Druggable pockets | LIGSITE / fpocket | doi:10.1016/S1093-3263(98)00002-3 |

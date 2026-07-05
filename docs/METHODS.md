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

The headline ranking (`cancer_tool/scoring.py`) fuses four axes into one explainable score:

```
priority = 0.30·recurrence + 0.35·pathogenicity + 0.20·druggability + 0.15·criticality
```

| Axis | Source | Meaning |
|------|--------|---------|
| **recurrence** | cancerhotspots.org | tumour count for the residue, normalised to the gene's top hotspot |
| **pathogenicity** | AlphaMissense (DeepMind) | probability the substitution is damaging |
| **druggability** | LIGSITE / fpocket | whether the residue lines a detected pocket, weighted by that pocket's score |
| **criticality** | ENM/NMA dynamics | structural importance — rigidity + hinge proximity, down-weighted by pLDDT |

**Pathogenicity** uses the exact substitution score when numbering checks out, otherwise the
positional mean (see *Numbering cross-check* below).

**Criticality** is
`min(1.0, 0.7·rigidity + 0.3·[near a hinge]) · confidence(pLDDT)`, where hinge proximity is a
±2-residue window and `confidence = clamp((pLDDT − 50) / 40, 0, 1)`. A residue in a
low-confidence (low-pLDDT) region has its structural signal discounted, since ENM on an
uncertain fold is unreliable.

### Weights are heuristics, not a trained model

The four weights are **expert-set** to reflect target-discovery priorities (pathogenicity and
recurrence lead; druggability and structure refine) — they are **not** fitted to clinical
outcomes. They are **robustness-tested rather than optimised**: `tests/test_scoring.py::`
`test_weights_perturbation_keeps_top_driver` perturbs every weight by ±0.05 and asserts the
top-ranked driver for the canonical genes is stable. The validation bar is rediscovering known
biology — TP53 → R273/R248/R175, KRAS → G12, BRAF → V600, EGFR → exon-19/L858/T790.

### Per-protein normalization (important caveat)

Two inputs are normalised **per protein**, not globally:

- **recurrence** = tumour count ÷ the gene's own top hotspot count;
- the **rigidity** underlying criticality is min-max normalised across the single structure.

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

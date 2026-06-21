# Changelog

## 0.3.0

Credibility and reproducibility hardening.

- **Confidence-aware dynamics.** Structural criticality is now scaled by AlphaFold
  pLDDT confidence (read from the Cα B-factor column), so low-confidence disordered
  regions no longer register as structurally critical. `compute_dynamics` returns a
  per-residue `plddt` array.
- **Automated biology validation.** `tests/test_validation.py` asserts the committed
  rankings rediscover known drivers (TP53 R273/R248/R175, KRAS G12, BRAF V600,
  EGFR L858/T790, …) and validates the data schema — the headline claim is now in CI.
- **Weight-sensitivity test.** Confirms a genuine driver tops the ranking under ±25%
  perturbation of each scoring axis.
- **Data provenance.** Each `data/{GENE}.json` records a `provenance` block
  (tool version, resolved AlphaFold model version, fetch date, sources).
- **Pinned dependencies** for reproducible dynamics output.
- **Honest limitations** surfaced in the dashboard and README; added `CITATION.cff`.

## 0.2.0

Target-discovery engine: ENM/NMA dynamics, AlphaMissense pathogenicity, LIGSITE-style
pocket detection, and the composite Target Priority Score across two front-ends.

"""Locks in the committed benchmark result (network-free).

Loads benchmark/results.json (produced by scripts/benchmark.py against live APIs)
and asserts the scientifically meaningful, *honest* properties — not an
aspirational "composite beats everything." The benchmark showed recurrence alone
is hard to beat on a recurrence-defined universe, so we do NOT assert the
composite beats recurrence; we assert what is actually true and worth defending:

  1. the composite gives a large lift over the random baseline,
  2. its top-ranked residues are almost all real drivers (precision@5),
  3. fusing genuinely helps versus the pathogenicity / druggability / criticality
     axes taken alone.

Run `python scripts/benchmark.py` to regenerate results.json when the pipeline
or gold standard changes.
"""

import json
from pathlib import Path

import pytest

RESULTS = Path(__file__).resolve().parent.parent / "benchmark" / "results.json"


@pytest.fixture(scope="module")
def results():
    if not RESULTS.exists():
        pytest.skip("benchmark/results.json not present — run scripts/benchmark.py")
    return json.loads(RESULTS.read_text())


def test_benchmark_schema(results):
    for key in ("generated", "tool_version", "evaluation", "per_gene", "n_genes_evaluated"):
        assert key in results, f"results.json missing '{key}'"
    ev = results["evaluation"]
    assert "ablation" in ev and "composite" in ev["ablation"]
    for axis in ("recurrence", "pathogenicity", "druggability", "criticality"):
        assert f"{axis}_only" in ev["ablation"]


def test_panel_is_substantial(results):
    ev = results["evaluation"]
    assert results["n_genes_evaluated"] >= 20
    assert ev["n_positives"] >= 50
    assert ev["n_residues"] >= 200


def test_composite_beats_random_baseline(results):
    ev = results["evaluation"]
    comp = ev["ablation"]["composite"]
    # A large, unambiguous lift over the prevalence a random ranker would score.
    assert comp["auprc"] >= 2.0 * ev["prevalence"]
    assert comp["auprc"] >= 0.55
    assert comp["auroc"] >= 0.70


def test_top_of_ranking_is_almost_all_drivers(results):
    comp = results["evaluation"]["ablation"]["composite"]
    assert comp["precision@5"] >= 0.8


def test_fusion_beats_the_non_recurrence_axes(results):
    """The composite must genuinely add value over pathogenicity, druggability,
    and criticality taken alone (the axes fusion is supposed to combine)."""
    ab = results["evaluation"]["ablation"]
    comp = ab["composite"]["auprc"]
    for axis in ("pathogenicity", "druggability", "criticality"):
        assert comp > ab[f"{axis}_only"]["auprc"], f"composite should beat {axis}-only"


def test_passengers_score_below_drivers(results):
    sep = results["evaluation"]["separation"]
    assert sep["positive_mean_score"] > sep["negative_mean_score"]

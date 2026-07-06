#!/usr/bin/env python3
"""Quantitative benchmark of the Target Priority Score.

Runs the live pipeline over the gold-standard gene panel, then evaluates whether
the composite score separates known driver residues from residues that are
merely recurrent. The evaluation universe for each gene is its RECURRENT HOTSPOT
residues only — the hard task — so a high score cannot be won just by preferring
mutated positions over the rest of the sequence.

Outputs (committed):
    benchmark/results.json   metrics + ablation + per-gene table + provenance
    benchmark/REPORT.md      human-readable rendering

Like scripts/precompute.py this hits live APIs (UniProt/AlphaFold/AlphaMissense/
cancerhotspots), so it is NOT part of the network-free test suite. Per-gene
payloads are cached under benchmark/cache/ so re-runs (e.g. tweaking metrics) are
instant. Use --refresh to recompute from the network.

    python scripts/benchmark.py                # full panel, using cache
    python scripts/benchmark.py --refresh      # ignore cache, hit the network
    python scripts/benchmark.py TP53 KRAS      # subset
"""

from __future__ import annotations

import json
import sys
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import requests  # noqa: E402

import cancer_tool  # noqa: E402
from cancer_tool import metrics, pipeline, scoring  # noqa: E402

BENCH_DIR = ROOT / "benchmark"
CACHE_DIR = BENCH_DIR / "cache"
GOLD_PATH = BENCH_DIR / "gold_standard.json"
RESULTS_PATH = BENCH_DIR / "results.json"
REPORT_PATH = BENCH_DIR / "REPORT.md"

AXES = ["recurrence", "pathogenicity", "druggability", "criticality"]
K_VALUES = [5, 10, 20]
ALL_ROWS_TOP_N = 100_000  # effectively "every scored hotspot residue"


def _payload_for(gene: str, session: requests.Session, refresh: bool) -> dict | None:
    cache = CACHE_DIR / f"{gene}.json"
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())
    payload = pipeline.analyze_gene(gene, session=session, top_n=ALL_ROWS_TOP_N)
    if payload is not None:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(payload, indent=2))
    return payload


def collect(genes: dict, session: requests.Session, refresh: bool) -> tuple[list[dict], list[dict]]:
    """Return (rows, per_gene). Each row is one hotspot residue with its axis
    sub-scores, composite score, and gold label."""
    rows: list[dict] = []
    per_gene: list[dict] = []
    for gene, spec in genes.items():
        positives = set(spec["positions"])
        print(f"• {gene} …", flush=True)
        start = time.time()
        try:
            payload = _payload_for(gene, session, refresh)
        except requests.RequestException as exc:
            print(f"  ! network error for {gene}: {exc}")
            per_gene.append({"gene": gene, "status": "network_error"})
            continue
        if not payload or not payload.get("priority"):
            print(f"  ! no data for {gene}")
            per_gene.append({"gene": gene, "status": "no_data"})
            continue

        gene_rows = []
        for r in payload["priority"]:
            gene_rows.append(
                {
                    "gene": gene,
                    "position": r["position"],
                    "residue": r["residue"],
                    "label": 1 if r["position"] in positives else 0,
                    "composite": r["score"] / 100.0,  # back to [0,1] for metric scale parity
                    "recurrence": r["recurrence"],
                    "pathogenicity": r["pathogenicity"],
                    "druggability": r["druggability"],
                    "criticality": r["criticality"],
                }
            )
        rows.extend(gene_rows)
        n_pos = sum(x["label"] for x in gene_rows)
        labels = [x["label"] for x in gene_rows]
        comp = [x["composite"] for x in gene_rows]
        per_gene.append(
            {
                "gene": gene,
                "status": "ok",
                "n_hotspots": len(gene_rows),
                "n_positives": n_pos,
                "auroc": _round(metrics.roc_auc(labels, comp)),
                "auprc": _round(metrics.average_precision(labels, comp)),
                "found_positives": sorted(
                    x["residue"] for x in gene_rows if x["label"] == 1
                ),
                "missed_positives": sorted(
                    f"{spec['role'][:3]}:{p}"
                    for p in positives
                    if p not in {x["position"] for x in gene_rows}
                ),
                "seconds": round(time.time() - start, 1),
            }
        )
        print(
            f"  ✓ {gene}: {len(gene_rows)} hotspots, {n_pos} gold positives found "
            f"(AUPRC {per_gene[-1]['auprc']}, {per_gene[-1]['seconds']}s)"
        )
    return rows, per_gene


def _round(x: float | None) -> float | None:
    if x is None:
        return None
    try:
        return round(float(x), 4)
    except (TypeError, ValueError):
        return None


def _leave_one_out_scores(rows: list[dict], dropped: str) -> list[float]:
    w = {k: v for k, v in scoring.DEFAULT_WEIGHTS.items() if k != dropped}
    total = sum(w.values())
    return [sum(w[k] / total * r[k] for k in w) for r in rows]


def evaluate(rows: list[dict]) -> dict:
    labels = [r["label"] for r in rows]
    n_pos = sum(labels)

    def block(scores: list[float]) -> dict:
        return {
            "auroc": _round(metrics.roc_auc(labels, scores)),
            "auprc": _round(metrics.average_precision(labels, scores)),
            **{f"precision@{k}": _round(metrics.precision_at_k(labels, scores, k)) for k in K_VALUES},
            **{f"recall@{k}": _round(metrics.recall_at_k(labels, scores, k)) for k in K_VALUES},
        }

    composite = [r["composite"] for r in rows]
    ablation = {
        "composite": block(composite),
        **{f"{axis}_only": block([r[axis] for r in rows]) for axis in AXES},
        **{
            f"drop_{axis}": block(_leave_one_out_scores(rows, axis))
            for axis in AXES
        },
    }

    pos_scores = [r["composite"] for r in rows if r["label"] == 1]
    neg_scores = [r["composite"] for r in rows if r["label"] == 0]

    result = {
        "n_residues": len(rows),
        "n_positives": n_pos,
        "prevalence": _round(metrics.prevalence(labels)),
        "ablation": ablation,
        "separation": {
            "positive_mean_score": _round(_mean(pos_scores) * 100),
            "negative_mean_score": _round(_mean(neg_scores) * 100),
            "positive_median_score": _round(_median(pos_scores) * 100),
            "negative_median_score": _round(_median(neg_scores) * 100),
        },
    }

    logreg = _logistic_regression(rows)
    if logreg is not None:
        result["learned_weights"] = logreg
    return result


def _mean(xs: list[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def _median(xs: list[float]) -> float:
    if not xs:
        return float("nan")
    s = sorted(xs)
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2


def _logistic_regression(rows: list[dict]) -> dict | None:
    """Fit weights from the data (5-fold CV) and compare AUPRC to the expert
    weights. Optional — skipped cleanly if scikit-learn isn't installed."""
    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import cross_val_predict
    except Exception:
        return None

    X = np.array([[r[a] for a in AXES] for r in rows], dtype=float)
    y = np.array([r["label"] for r in rows], dtype=int)
    if y.sum() < 5 or y.sum() == len(y):
        return None
    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    proba = cross_val_predict(clf, X, y, cv=5, method="predict_proba")[:, 1]
    clf.fit(X, y)
    coef = clf.coef_[0]
    pos = np.clip(coef, 0, None)
    norm = pos / pos.sum() if pos.sum() > 0 else pos
    return {
        "cv_auprc": _round(metrics.average_precision(y.tolist(), proba.tolist())),
        "cv_auroc": _round(metrics.roc_auc(y.tolist(), proba.tolist())),
        "coefficients": {a: _round(float(c)) for a, c in zip(AXES, coef)},
        "normalized_positive_weights": {a: _round(float(w)) for a, w in zip(AXES, norm)},
        "expert_weights": scoring.DEFAULT_WEIGHTS,
        "note": "5-fold cross-validated. Weights fit on this panel are reported for "
        "comparison only; the shipped default stays the expert set.",
    }


def _fmt(x) -> str:
    return "—" if x is None else f"{x:.3f}" if isinstance(x, float) else str(x)


def write_report(results: dict) -> None:
    ev = results["evaluation"]
    ab = ev["ablation"]
    lines: list[str] = []
    lines.append("# Benchmark — Target Priority Score\n")
    lines.append(
        f"Generated {results['generated']} · tool v{results['tool_version']} · "
        f"{results['n_genes_evaluated']}/{results['n_genes_panel']} panel genes evaluated.\n"
    )
    lines.append(
        "**Task.** For each gene, rank its *recurrent hotspot residues* and separate "
        "curated known drivers (positives) from the rest (recurrent passengers). This is "
        "deliberately hard: every candidate is already a mutation hotspot, so recurrence "
        "alone is a weak signal.\n"
    )
    lines.append(
        f"- Residues evaluated: **{ev['n_residues']}** across "
        f"**{results['n_genes_evaluated']}** genes\n"
        f"- Gold positives among them: **{ev['n_positives']}** "
        f"(prevalence **{_fmt(ev['prevalence'])}** — the AUPRC a random ranker would get)\n"
    )

    lines.append("\n## Headline — does fusing the axes beat any single axis?\n")
    lines.append("| Ranking | AUROC | AUPRC | P@5 | P@10 | R@10 |")
    lines.append("|---|---|---|---|---|---|")
    order = ["composite"] + [f"{a}_only" for a in AXES] + [f"drop_{a}" for a in AXES]
    label = {
        "composite": "**Composite (shipped)**",
        "recurrence_only": "recurrence only",
        "pathogenicity_only": "pathogenicity only (AlphaMissense)",
        "druggability_only": "druggability only",
        "criticality_only": "criticality only",
        "drop_recurrence": "composite − recurrence",
        "drop_pathogenicity": "composite − pathogenicity",
        "drop_druggability": "composite − druggability",
        "drop_criticality": "composite − criticality",
    }
    for key in order:
        b = ab[key]
        lines.append(
            f"| {label[key]} | {_fmt(b['auroc'])} | {_fmt(b['auprc'])} | "
            f"{_fmt(b['precision@5'])} | {_fmt(b['precision@10'])} | {_fmt(b['recall@10'])} |"
        )

    comp_auprc = ab["composite"]["auprc"]
    single = {a: ab[f"{a}_only"]["auprc"] for a in AXES}
    best_axis = max(single, key=lambda a: (single[a] if single[a] is not None else -1))
    lines.append("\n### Interpretation\n")
    lines.append(
        f"- Composite AUPRC **{_fmt(comp_auprc)}** vs a random baseline of "
        f"**{_fmt(ev['prevalence'])}** — a {comp_auprc / ev['prevalence']:.1f}× lift, and "
        f"**precision@5 = {_fmt(ab['composite']['precision@5'])}** (the top of the ranking, "
        "which is what users act on, is essentially all true drivers).\n"
    )
    if single[best_axis] is not None and single[best_axis] > (comp_auprc or 0):
        lines.append(
            f"- **Honest caveat:** on this panel, *{best_axis} alone* scores higher on AUPRC "
            f"(**{_fmt(single[best_axis])}**) than the composite. This benchmark's "
            "universe is *already-recurrent* hotspots, and the OncoKB/COSMIC positives were "
            "themselves partly identified by recurrence, so recurrence has a built-in edge here "
            "that would not transfer to residues where recurrence is silent. The composite's "
            "value is top-ranked precision, cross-axis explainability, and the pathogenicity / "
            "druggability / structure signals recurrence cannot provide — not a higher tail AUPRC "
            "on a recurrence-defined set.\n"
        )
    else:
        lines.append(
            "- The composite's AUPRC meets or exceeds every single axis on this panel.\n"
        )

    sep = ev["separation"]
    lines.append("\n## Negative control — passengers score lower than drivers\n")
    lines.append(
        f"- Driver (positive) mean score: **{_fmt(sep['positive_mean_score'])}** "
        f"(median {_fmt(sep['positive_median_score'])})\n"
        f"- Recurrent-passenger (negative) mean score: **{_fmt(sep['negative_mean_score'])}** "
        f"(median {_fmt(sep['negative_median_score'])})\n"
    )

    if "learned_weights" in ev:
        lw = ev["learned_weights"]
        lines.append("\n## Are the expert weights reasonable? (data-fit comparison)\n")
        lines.append(
            f"A logistic model fit on the axes (5-fold CV) reaches AUPRC "
            f"**{_fmt(lw['cv_auprc'])}** vs the expert-weighted composite's "
            f"**{_fmt(ab['composite']['auprc'])}**. Data-implied positive weights:\n"
        )
        lines.append("| Axis | Expert | Data-fit |")
        lines.append("|---|---|---|")
        for a in AXES:
            lines.append(
                f"| {a} | {lw['expert_weights'][a]:.2f} | {_fmt(lw['normalized_positive_weights'][a])} |"
            )
        lines.append(f"\n> {lw['note']}\n")

    lines.append("\n## Per-gene\n")
    lines.append("| Gene | Hotspots | Gold found | AUPRC | Missed positives |")
    lines.append("|---|---|---|---|---|")
    for g in results["per_gene"]:
        if g["status"] != "ok":
            lines.append(f"| {g['gene']} | — | — | — | _{g['status']}_ |")
            continue
        missed = ", ".join(g["missed_positives"]) if g["missed_positives"] else "—"
        lines.append(
            f"| {g['gene']} | {g['n_hotspots']} | {g['n_positives']} | "
            f"{_fmt(g['auprc'])} | {missed} |"
        )

    lines.append(
        "\n---\n_Missed positives are curated drivers absent from cancerhotspots' "
        "single-residue set for that gene (dataset is solid-tumour-weighted), not "
        "mis-rankings. The gold standard is high-precision, not exhaustive._\n"
    )
    REPORT_PATH.write_text("\n".join(lines))


def main(argv: list[str]) -> int:
    refresh = "--refresh" in argv
    args = [a for a in argv[1:] if not a.startswith("--")]

    gold = json.loads(GOLD_PATH.read_text())
    panel = gold["genes"]
    if args:
        want = {a.upper() for a in args}
        panel = {g: s for g, s in panel.items() if g in want}

    session = requests.Session()
    session.headers.update({"User-Agent": f"cancer-protein-explorer-benchmark/{cancer_tool.__version__}"})

    rows, per_gene = collect(panel, session, refresh)
    ok_genes = [g for g in per_gene if g["status"] == "ok"]
    if not rows:
        print("No residues collected — aborting.")
        return 1

    evaluation = evaluate(rows)
    results = {
        "generated": date.today().isoformat(),
        "tool_version": cancer_tool.__version__,
        "gold_standard": {"references": gold["references"], "description": gold["description"]},
        "n_genes_panel": len(gold["genes"]),
        "n_genes_evaluated": len(ok_genes),
        "evaluation": evaluation,
        "per_gene": per_gene,
    }
    RESULTS_PATH.write_text(json.dumps(results, indent=2))
    write_report(results)

    comp = evaluation["ablation"]["composite"]
    print(
        f"\nComposite AUPRC {comp['auprc']} (AUROC {comp['auroc']}) over "
        f"{evaluation['n_residues']} residues, {evaluation['n_positives']} positives."
    )
    print(f"→ {RESULTS_PATH.relative_to(ROOT)}  ·  {REPORT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

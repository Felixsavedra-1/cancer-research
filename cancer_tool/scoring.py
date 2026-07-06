"""Target Priority Score — the composite target-discovery ranking.

Pure, network-free fusion of four axes into one explainable score:

    0.30·recurrence + 0.35·pathogenicity + 0.20·druggability + 0.15·criticality

Weights are expert-set heuristics, not trained. recurrence and criticality are
normalised per protein, so scores rank residues within a protein but do not
compare across proteins. See docs/METHODS.md.
"""

from __future__ import annotations

from . import dynamics as dyn
from . import pathogenicity as patho
from . import pockets as pock

# Mirrored in cancer-explorer.html's glossary text (the GLOSS "Target Priority
# Score" entry near its scoring comment) — update both.
DEFAULT_WEIGHTS = {
    "recurrence": 0.30,
    "pathogenicity": 0.35,
    "druggability": 0.20,
    "criticality": 0.15,
}

HINGE_WINDOW = 2

# The [0,100] scale in ``score_residues`` only holds if the weights sum to 1, so
# custom weights are validated rather than silently rescaling every score.
_REQUIRED_AXES = frozenset(DEFAULT_WEIGHTS)


def _validate_weights(weights: dict) -> dict:
    missing = _REQUIRED_AXES - weights.keys()
    if missing:
        raise ValueError(f"weights missing axes: {sorted(missing)}")
    total = sum(weights[k] for k in _REQUIRED_AXES)
    if abs(total - 1.0) > 1e-6:
        raise ValueError(
            f"weights must sum to 1.0 for the 0-100 scale to hold, got {total:.4f}"
        )
    return weights


def _most_common_variant(variants: dict) -> str | None:
    if not variants:
        return None
    return max(variants.items(), key=lambda kv: int(kv[1]))[0]


def _confidence(plddt: float | None) -> float:
    if plddt is None:
        return 1.0
    return min(1.0, max(0.0, (plddt - 50.0) / 40.0))


def _criticality(
    position: int,
    rigidity: dict[int, float],
    hinges: set[int],
    plddt: dict[int, float],
) -> float:
    rigid = rigidity.get(position, 0.0)
    near_hinge = any(abs(position - h) <= HINGE_WINDOW for h in hinges)
    base = min(1.0, 0.7 * rigid + (0.3 if near_hinge else 0.0))
    return base * _confidence(plddt.get(position))


def _rationale(components: dict, n_tumours: int, am_label: str | None) -> str:
    parts: list[str] = []
    if components["recurrence"] >= 0.5:
        parts.append(f"recurrently mutated ({n_tumours} tumours)")
    elif n_tumours:
        parts.append(f"mutated in {n_tumours} tumours")
    if components["pathogenicity"] >= 0.7:
        parts.append(f"predicted pathogenic ({am_label or 'AlphaMissense'})")
    elif components["pathogenicity"] >= 0.4:
        parts.append("moderately damaging")
    if components["druggability"] > 0:
        parts.append("lines a druggable pocket")
    if components["criticality"] >= 0.6:
        parts.append("in a structurally critical site")
    return "; ".join(parts) if parts else "weak signal across all axes"


def score_residues(
    hotspots: list[dict],
    pathogenicity: dict | None = None,
    dynamics: dict | None = None,
    pockets: list[dict] | None = None,
    sequence: str = "",
    weights: dict | None = None,
) -> list[dict]:
    """Rank hotspot residues by the composite Target Priority Score.

    Fuses recurrence, AlphaMissense pathogenicity, pocket druggability, and
    ENM/NMA criticality (see module docstring for the formula and weights). Every
    optional signal degrades gracefully to 0 when absent. Returns rows sorted by
    descending score, each carrying its sub-scores, a ``numbering_ok`` flag, and a
    plain-English ``rationale``.
    """
    weights = _validate_weights(weights or DEFAULT_WEIGHTS)
    pathogenicity = pathogenicity or {}
    pockets = pockets or []

    rigidity = dyn.rigidity_by_position(dynamics) if dynamics else {}
    hinges = set(dynamics.get("hinges", [])) if dynamics else set()
    plddt = dyn.plddt_by_position(dynamics) if dynamics else {}

    max_count = max((h.get("count", 0) for h in hotspots), default=0)

    rows: list[dict] = []
    for hot in hotspots:
        position = hot["position"]
        n_tumours = int(hot.get("count", 0))
        recurrence = (n_tumours / max_count) if max_count else 0.0

        wt = sequence[position - 1] if 0 < position <= len(sequence) else ""

        # cancerhotspots positions may use a different transcript than
        # AlphaMissense's canonical UniProt numbering; if the wild-type residues
        # disagree, flag the row and fall back to the positional mean.
        am_wt = patho.wt_at_position(pathogenicity, position)
        numbering_ok = (not wt) or (am_wt is None) or (am_wt == wt)

        am_label = None
        am_class = None
        path_score = patho.position_pathogenicity(pathogenicity, position)
        top_mut = _most_common_variant(hot.get("variants", {}))
        if wt and top_mut and numbering_ok:
            variant = patho.variant_score(pathogenicity, f"{wt}{position}{top_mut}")
            if variant:
                path_score = variant["score"]
                am_label = variant["class_label"]
                am_class = variant["class"]

        druggability = pock.pocket_proximity(position, pockets)
        criticality = _criticality(position, rigidity, hinges, plddt)

        components = {
            "recurrence": recurrence,
            "pathogenicity": path_score,
            "druggability": druggability,
            "criticality": criticality,
        }
        score = 100.0 * sum(weights[k] * components[k] for k in weights)

        rows.append(
            {
                "position": position,
                "residue": hot.get("residue", f"{wt}{position}"),
                "wt": wt,
                "score": round(score, 1),
                "recurrence": round(recurrence, 3),
                "pathogenicity": round(path_score, 3),
                "druggability": round(druggability, 3),
                "criticality": round(criticality, 3),
                "tumours": n_tumours,
                "am_class": am_class,
                "numbering_ok": numbering_ok,
                "rationale": _rationale(components, n_tumours, am_label),
            }
        )

    rows.sort(key=lambda r: r["score"], reverse=True)
    return rows

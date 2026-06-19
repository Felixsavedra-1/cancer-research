"""Target Priority Score — fuse four orthogonal axes into one ranked, explainable shortlist.

Axes: recurrence (cancerhotspots), pathogenicity (AlphaMissense), druggability (pocket
geometry), and structural criticality (ENM/NMA rigidity + hinge proximity), tuned to
surface druggable cancer drivers. Pure: same inputs → same ranking, no I/O.
"""

from __future__ import annotations

from . import dynamics as dyn
from . import pathogenicity as patho
from . import pockets as pock

# Weights for the druggable-cancer-driver objective. Pathogenicity and recurrence
# lead (is it a real driver?); druggability gates actionability; criticality breaks
# ties toward structurally pivotal sites. Overridable via ``score_residues(weights=)``.
DEFAULT_WEIGHTS = {
    "recurrence": 0.30,
    "pathogenicity": 0.35,
    "druggability": 0.20,
    "criticality": 0.15,
}

# Distance (in residues) within which a site counts as "at a hinge".
HINGE_WINDOW = 2


def _most_common_variant(variants: dict) -> str | None:
    """Return the single-letter substituted residue mutated in the most tumours."""
    if not variants:
        return None
    return max(variants.items(), key=lambda kv: int(kv[1]))[0]


def _criticality(position: int, rigidity: dict[int, float], hinges: set[int]) -> float:
    """Structural load-bearing signal in 0–1: rigid core + hinge proximity.

    NMA rigidity is a strong proxy for burial — buried core residues barely
    fluctuate — so a separate solvent-accessibility term would be largely redundant.
    Hinge proximity adds the functionally critical domain-pivot sites that aren't
    necessarily the most rigid.
    """
    rigid = rigidity.get(position, 0.0)
    near_hinge = any(abs(position - h) <= HINGE_WINDOW for h in hinges)
    return min(1.0, 0.7 * rigid + (0.3 if near_hinge else 0.0))


def _rationale(components: dict, n_tumours: int, am_label: str | None) -> str:
    """Human-readable 'why this ranks here', so the score is never a black box."""
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
    """Rank cancer-hotspot residues by druggable-driver priority.

    Candidates are the gene's recurrent-mutation sites (``hotspots``); each is scored
    on the four axes above and sorted high→low. Any missing signal (e.g. a gene with
    no AlphaMissense coverage, or pocket detection unavailable) simply contributes 0
    for that axis, so the tool degrades gracefully rather than failing.

    Each returned row::

        {"position", "residue", "wt", "score",          # 0-100 composite
         "recurrence", "pathogenicity", "druggability", "criticality",  # 0-1 axes
         "tumours", "am_class", "rationale"}
    """
    weights = weights or DEFAULT_WEIGHTS
    pathogenicity = pathogenicity or {}
    pockets = pockets or []

    rigidity = dyn.rigidity_by_position(dynamics) if dynamics else {}
    hinges = set(dynamics.get("hinges", [])) if dynamics else set()

    max_count = max((h.get("count", 0) for h in hotspots), default=0)

    rows: list[dict] = []
    for hot in hotspots:
        position = hot["position"]
        n_tumours = int(hot.get("count", 0))
        recurrence = (n_tumours / max_count) if max_count else 0.0

        wt = sequence[position - 1] if 0 < position <= len(sequence) else ""

        # Prefer the exact cancer variant's pathogenicity; fall back to the
        # residue's mean intolerance across all substitutions.
        am_label = None
        am_class = None
        path_score = patho.position_pathogenicity(pathogenicity, position)
        top_mut = _most_common_variant(hot.get("variants", {}))
        if wt and top_mut:
            variant = patho.variant_score(pathogenicity, f"{wt}{position}{top_mut}")
            if variant:
                path_score = variant["score"]
                am_label = variant["class_label"]
                am_class = variant["class"]

        druggability = pock.pocket_proximity(position, pockets)
        criticality = _criticality(position, rigidity, hinges)

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
                "rationale": _rationale(components, n_tumours, am_label),
            }
        )

    rows.sort(key=lambda r: r["score"], reverse=True)
    return rows

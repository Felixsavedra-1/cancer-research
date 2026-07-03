from __future__ import annotations

import re

import requests

_MUTATION_RE = re.compile(r"^\s*([A-Za-z])\s*(\d+)\s*([A-Za-z*])\s*$")
CANCER_HOTSPOTS = "https://www.cancerhotspots.org/api/hotspots/single"


def parse_mutation(text: str) -> dict | None:
    match = _MUTATION_RE.match(text or "")
    if not match:
        return None
    wt, pos, mut = match.group(1).upper(), int(match.group(2)), match.group(3).upper()
    return {"wt": wt, "position": pos, "mut": mut, "label": f"{wt}{pos}{mut}"}


def parse_mutations(text: str) -> tuple[list[dict], list[str]]:
    tokens = [t for t in re.split(r"[,\s]+", text or "") if t]
    parsed, bad = [], []
    for token in tokens:
        mutation = parse_mutation(token)
        (parsed if mutation else bad).append(mutation or token)
    return parsed, bad


def validate_mutation(mutation: dict, sequence: str) -> tuple[bool, str]:
    pos = mutation["position"]
    if pos < 1 or pos > len(sequence):
        return False, f"Position {pos} is outside the protein (length {len(sequence)})."
    actual = sequence[pos - 1]
    if actual != mutation["wt"]:
        return (
            False,
            f"Sequence has {actual} at position {pos}, not {mutation['wt']} — "
            "check the residue letter or numbering.",
        )
    return True, "ok"


def fetch_hotspots(gene: str, session: requests.Session | None = None) -> list[dict]:
    http = session or requests
    resp = http.get(CANCER_HOTSPOTS, headers={"Accept": "application/json"}, timeout=60)
    if not resp.ok:
        return []
    try:
        data = resp.json()
    except ValueError:
        return []

    gene = gene.upper()
    hotspots = []
    for entry in data or []:
        if (entry.get("hugoSymbol") or "").upper() != gene:
            continue
        residue = entry.get("residue", "")
        match = re.search(r"(\d+)", residue)
        if not match:
            continue
        variants = entry.get("variantAminoAcid", {}) or {}
        count = entry.get("tumorCount")
        if count is None:
            count = sum(int(v) for v in variants.values()) if variants else 0
        hotspots.append(
            {
                "residue": residue,
                "position": int(match.group(1)),
                "count": int(count),
                "variants": variants,
            }
        )
    hotspots.sort(key=lambda h: h["count"], reverse=True)
    return hotspots

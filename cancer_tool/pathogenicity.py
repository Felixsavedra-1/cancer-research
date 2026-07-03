from __future__ import annotations

import csv
import re

import requests

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

_VARIANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z*])$")

CLASS_LABELS = {
    "LPath": "Likely pathogenic",
    "Amb": "Ambiguous",
    "LBen": "Likely benign",
}


def parse_alphamissense_csv(text: str) -> dict:
    by_variant: dict[str, dict] = {}
    acc: dict[int, list[float]] = {}

    reader = csv.reader(text.splitlines())
    header = next(reader, None)
    if header and header[0].strip().lower() != "protein_variant":
        rows = [header, *reader]
    else:
        rows = reader

    for row in rows:
        if len(row) < 3:
            continue
        variant, score_str, cls = row[0].strip(), row[1].strip(), row[2].strip()
        match = _VARIANT_RE.match(variant)
        if not match:
            continue
        try:
            score = float(score_str)
        except ValueError:
            continue
        pos = int(match.group(2))
        by_variant[variant] = {
            "score": score,
            "class": cls,
            "class_label": CLASS_LABELS.get(cls, cls),
        }
        acc.setdefault(pos, []).append(score)

    by_position = {
        pos: {
            "mean": round(sum(scores) / len(scores), 4),
            "max": round(max(scores), 4),
            "n": len(scores),
        }
        for pos, scores in acc.items()
    }
    return {"by_position": by_position, "by_variant": by_variant}


def fetch_alphamissense(
    accession: str, session: requests.Session | None = None
) -> dict:
    http = session or requests
    empty = {"by_position": {}, "by_variant": {}}
    try:
        meta = http.get(ALPHAFOLD_API.format(accession=accession), timeout=30)
        if not meta.ok:
            return empty
        payload = meta.json()
        if not payload:
            return empty
        url = payload[0].get("amAnnotationsUrl")
        if not url:
            return empty
        csv_resp = http.get(url, timeout=60)
        if not csv_resp.ok:
            return empty
    except (requests.RequestException, ValueError, KeyError, IndexError):
        return empty

    return parse_alphamissense_csv(csv_resp.text)


def variant_score(pathogenicity: dict, label: str) -> dict | None:
    return pathogenicity.get("by_variant", {}).get(label.upper())


def position_pathogenicity(pathogenicity: dict, position: int) -> float:
    rec = pathogenicity.get("by_position", {}).get(position)
    return float(rec["mean"]) if rec else 0.0

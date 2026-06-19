"""Target context from the Open Targets Platform GraphQL API: clinical tractability
(modalities) and top disease associations. Best-effort enrichment — any failure yields
``None`` and never blocks the page.
"""

from __future__ import annotations

import requests

OPEN_TARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"

_SEARCH_QUERY = """
query Search($q: String!) {
  search(queryString: $q, entityNames: ["target"]) {
    hits { id name }
  }
}
"""

_TARGET_QUERY = """
query Target($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    tractability { label modality value }
    associatedDiseases(page: { index: 0, size: 5 }) {
      rows { disease { name } score }
    }
  }
}
"""

# Modality codes Open Targets uses in tractability buckets.
_MODALITY_LABELS = {"SM": "Small molecule", "AB": "Antibody", "PR": "PROTAC", "OC": "Other"}


def _post(query: str, variables: dict, session: requests.Session | None) -> dict | None:
    http = session or requests
    resp = http.post(
        OPEN_TARGETS_API, json={"query": query, "variables": variables}, timeout=30
    )
    if not resp.ok:
        return None
    payload = resp.json()
    if payload.get("errors"):
        return None
    return payload.get("data")


def fetch_target_context(
    gene: str, session: requests.Session | None = None
) -> dict | None:
    """Return tractability + top disease associations for a gene, or ``None``.

    Result shape::

        {"ensembl_id": "ENSG...",
         "symbol": "TP53",
         "tractability": [{"modality": "Small molecule", "label": "Advanced Clinical"}, ...],
         "top_diseases": [{"name": "Li-Fraumeni syndrome", "score": 0.88}, ...]}

    Only tractability buckets that are *true* (the assessment passed) are returned.
    """
    try:
        found = _post(_SEARCH_QUERY, {"q": gene}, session)
        hits = (found or {}).get("search", {}).get("hits", []) if found else []
        if not hits:
            return None
        ensembl_id = hits[0]["id"]

        data = _post(_TARGET_QUERY, {"id": ensembl_id}, session)
        target = (data or {}).get("target") if data else None
        if not target:
            return None
    except (requests.RequestException, ValueError, KeyError, IndexError):
        return None

    tractability = [
        {
            "modality": _MODALITY_LABELS.get(t["modality"], t["modality"]),
            "label": t["label"],
        }
        for t in target.get("tractability", [])
        if t.get("value")
    ]
    top_diseases = [
        {"name": row["disease"]["name"], "score": round(float(row["score"]), 3)}
        for row in target.get("associatedDiseases", {}).get("rows", [])
    ]
    return {
        "ensembl_id": ensembl_id,
        "symbol": target.get("approvedSymbol", gene.upper()),
        "tractability": tractability,
        "top_diseases": top_diseases,
    }

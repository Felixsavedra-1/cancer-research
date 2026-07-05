"""Live analysis API for the Cancer Protein Explorer.

Runs the ``cancer_tool`` pipeline for any human gene and returns the same schema
as the committed ``data/{GENE}.json`` files, so the static HTML can fall back to
it for genes that aren't precomputed. Committed genes are served from disk;
everything else is computed once and disk-cached under ``data/cache/``.

    uvicorn api.main:app --reload --port 8100
"""

from __future__ import annotations

import json
import threading
from pathlib import Path

import requests
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

import cancer_tool
from cancer_tool import pipeline

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
CACHE_DIR = DATA_DIR / "cache"

app = FastAPI(
    title="Cancer Protein Explorer API",
    version=cancer_tool.__version__,
    summary="Live target-discovery analysis for any human gene.",
)
# The static HTML is served from another origin (file://, Pages, http.server).
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)

# One shared HTTP session; the lock serialises the CPU-heavy ProDy/pocket work.
_session = requests.Session()
_session.headers.update({"User-Agent": f"cancer-protein-explorer/{cancer_tool.__version__}"})
_compute_lock = threading.Lock()


def _cached_path(gene: str) -> Path | None:
    """Return a committed or previously-cached JSON path for ``gene``, if any."""
    committed = DATA_DIR / f"{gene}.json"
    if committed.exists():
        return committed
    cached = CACHE_DIR / f"{gene}.json"
    if cached.exists():
        return cached
    return None


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "version": cancer_tool.__version__}


@app.get("/analyze/{gene}")
def analyze(gene: str) -> dict:
    """Full analysis payload for ``gene`` (committed → cache → live compute)."""
    gene = gene.strip().upper()
    if not gene.isalnum():
        raise HTTPException(status_code=400, detail="Gene symbol must be alphanumeric.")

    existing = _cached_path(gene)
    if existing is not None:
        return json.loads(existing.read_text())

    with _compute_lock:
        # Re-check under the lock in case a concurrent request just computed it.
        existing = _cached_path(gene)
        if existing is not None:
            return json.loads(existing.read_text())
        try:
            payload = pipeline.analyze_gene(gene, session=_session)
        except requests.RequestException as exc:
            raise HTTPException(status_code=502, detail=f"Upstream data error: {exc}") from exc
        if payload is None:
            raise HTTPException(status_code=404, detail=f"No reviewed human protein for '{gene}'.")
        pipeline.validate_payload(payload)
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        (CACHE_DIR / f"{gene}.json").write_text(json.dumps(payload, indent=2))
        return payload

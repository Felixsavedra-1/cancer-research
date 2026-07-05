#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import requests  # noqa: E402

from cancer_tool import pipeline  # noqa: E402

FEATURED_GENES = ["TP53", "KRAS", "BRAF", "KIT", "FLT3", "IDH2", "PTEN", "EGFR"]
TOP_N = pipeline.DEFAULT_TOP_N
DATA_DIR = ROOT / "data"


def build_gene(gene: str, session: requests.Session) -> dict | None:
    """Run the shared pipeline for a gene, noting missing records and numbering
    mismatches on the way."""
    payload = pipeline.analyze_gene(gene, session=session, top_n=TOP_N)
    if payload is None:
        print(f"  ! no UniProt record for {gene}")
        return None
    mismatches = [r["residue"] for r in payload["priority"] if r.get("numbering_ok") is False]
    if mismatches:
        print(f"  ! {gene}: {len(mismatches)} numbering mismatch(es): {', '.join(mismatches[:8])}")
    return payload


def main(argv: list[str]) -> int:
    genes = [g.upper() for g in argv[1:]] or FEATURED_GENES
    DATA_DIR.mkdir(exist_ok=True)
    session = requests.Session()
    session.headers.update({"User-Agent": "cancer-protein-explorer-precompute"})

    index = []
    for gene in genes:
        print(f"• {gene} …", flush=True)
        start = time.time()
        try:
            payload = build_gene(gene, session)
        except requests.RequestException as exc:
            print(f"  ! network error for {gene}: {exc}")
            continue
        if not payload:
            continue
        pipeline.validate_payload(payload)
        out = DATA_DIR / f"{gene}.json"
        out.write_text(json.dumps(payload, indent=2))
        top_residue = payload["priority"][0]["residue"] if payload["priority"] else None
        index.append({"gene": gene, "name": payload["name"], "top": top_residue})
        print(
            f"  ✓ {gene}: {len(payload['priority'])} ranked residues "
            f"({time.time() - start:.1f}s) → {out.relative_to(ROOT)}"
        )

    (DATA_DIR / "index.json").write_text(
        json.dumps({"genes": index, "generated": date.today().isoformat()}, indent=2)
    )
    print(f"\nWrote {len(index)} gene file(s) to {DATA_DIR.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

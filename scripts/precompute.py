#!/usr/bin/env python3
"""Precompute the full engine for featured genes → ``data/{GENE}.json`` for the
browser-only ``cancer-explorer.html`` (which can't run the Python science live).

Usage::

    python scripts/precompute.py                 # all featured genes
    python scripts/precompute.py TP53 KRAS BRAF   # specific genes
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

from cancer_tool import (  # noqa: E402
    dynamics,
    mutations,
    pathogenicity,
    pockets,
    scoring,
    structures,
    targets,
    uniprot,
)

FEATURED_GENES = ["TP53", "KRAS", "BRAF", "KIT", "FLT3", "IDH2", "PTEN", "EGFR"]
TOP_N = 30
DATA_DIR = ROOT / "data"


def build_gene(gene: str, session: requests.Session) -> dict | None:
    protein = uniprot.get_protein(gene, session=session)
    if not protein:
        print(f"  ! no UniProt record for {gene}")
        return None

    pdb = structures.fetch_alphafold_pdb(protein["accession"], session=session)
    hotspots = mutations.fetch_hotspots(gene, session=session)
    pathos = pathogenicity.fetch_alphamissense(protein["accession"], session=session)
    try:
        dyn = dynamics.compute_dynamics(pdb)
    except dynamics.DynamicsError:
        dyn = None
    pocket_list = pockets.detect_pockets(pdb)
    rows = scoring.score_residues(hotspots, pathos, dyn, pocket_list, protein["sequence"])
    context = targets.fetch_target_context(gene, session=session)

    return {
        "gene": gene,
        "accession": protein["accession"],
        "name": protein["name"],
        "length": protein["length"],
        "generated": date.today().isoformat(),
        "priority": rows[:TOP_N],
        "dynamics": (
            {
                "residue_numbers": dyn["residue_numbers"],
                "flexibility": dyn["flexibility"],
                "hinges": dyn["hinges"],
                "n_modes": dyn["n_modes"],
            }
            if dyn
            else None
        ),
        "pockets": [
            {k: p[k] for k in ("druggability", "volume", "residues", "source")}
            for p in pocket_list[:3]
        ],
        "targets": context,
    }


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

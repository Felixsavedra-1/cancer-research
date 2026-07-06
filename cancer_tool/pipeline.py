"""End-to-end target-discovery pipeline for one gene.

Single source of truth for a gene's analysis payload: both
``scripts/precompute.py`` and the live API (``api/main.py``) call
:func:`analyze_gene`, so the committed ``data/{GENE}.json`` and live responses
share one schema. See docs/METHODS.md for the pipeline stages.
"""

from __future__ import annotations

import re
import subprocess
from datetime import date
from functools import lru_cache
from importlib import metadata

import requests

import cancer_tool
from cancer_tool import (
    dynamics,
    mutations,
    pathogenicity,
    pockets,
    scoring,
    structures,
    targets,
    uniprot,
)

DEFAULT_TOP_N = 30
SOURCES = ["UniProt", "AlphaFold DB", "AlphaMissense", "cancerhotspots.org", "Open Targets"]

# Packages that determine the numerical output — pinned in provenance so a
# ranking can be traced to the exact scientific stack that produced it.
_PROVENANCE_PACKAGES = ("prody", "numpy", "scipy")


@lru_cache(maxsize=1)
def _git_commit() -> str:
    """Short hash of the code that generated a payload, or 'unknown' outside a repo."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=cancer_tool.__path__[0],
        )
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


@lru_cache(maxsize=1)
def _package_versions() -> dict[str, str]:
    versions = {}
    for pkg in _PROVENANCE_PACKAGES:
        try:
            versions[pkg] = metadata.version(pkg)
        except Exception:
            versions[pkg] = "unknown"
    return versions


def analyze_gene(
    gene: str,
    session: requests.Session | None = None,
    top_n: int = DEFAULT_TOP_N,
) -> dict | None:
    """Run the full pipeline for ``gene`` and return its analysis payload.

    Returns the same schema as a committed ``data/{GENE}.json`` (validated by
    :func:`validate_payload`), or ``None`` if the gene has no reviewed UniProt
    record. Network failures propagate to the caller.
    """
    gene = gene.upper()
    session = session or requests.Session()

    protein = uniprot.get_protein(gene, session=session)
    if not protein:
        return None

    model_url = structures.resolve_pdb_url(protein["accession"], session=session)
    pdb = structures.fetch_alphafold_pdb(protein["accession"], session=session)
    version_match = re.search(r"-model_(v\d+)", model_url)
    model_version = version_match.group(1) if version_match else "unknown"

    hotspots = mutations.fetch_hotspots(gene, session=session)
    pathos = pathogenicity.fetch_alphamissense(protein["accession"], session=session)
    try:
        dyn = dynamics.compute_dynamics(pdb)
    except dynamics.DynamicsError:
        dyn = None
    # Like dynamics, pocket detection must degrade to nothing rather than break
    # the ranking (missing ProDy/scipy or an unparseable structure raise here).
    try:
        pocket_list = pockets.detect_pockets(pdb)
    except Exception:
        pocket_list = []
    rows = scoring.score_residues(hotspots, pathos, dyn, pocket_list, protein["sequence"])
    context = targets.fetch_target_context(gene, session=session)

    return {
        "gene": gene,
        "accession": protein["accession"],
        "name": protein["name"],
        "length": protein["length"],
        "provenance": {
            "tool_version": cancer_tool.__version__,
            "git_commit": _git_commit(),
            "packages": _package_versions(),
            "alphafold_model": model_version,
            "fetched": date.today().isoformat(),
            "sources": SOURCES,
            "enm_params": dyn["params"] if dyn else None,
            "weights": scoring.DEFAULT_WEIGHTS,
        },
        "priority": rows[:top_n],
        "dynamics": (
            {
                "residue_numbers": dyn["residue_numbers"],
                "plddt": dyn["plddt"],
                "flexibility": dyn["flexibility"],
                "hinges": dyn["hinges"],
                "collectivity": dyn["collectivity"],
                "n_modes": dyn["n_modes"],
                "params": dyn["params"],
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


def validate_payload(payload: dict) -> None:
    """Raise ``ValueError`` unless ``payload`` has the required analysis schema.

    Provenance is validated for content, not just presence, so a payload can
    always be traced to the code + scientific-stack that produced it.
    """
    for key in ("gene", "accession", "name", "length", "provenance", "priority"):
        if key not in payload:
            raise ValueError(f"{payload.get('gene', '?')}: missing '{key}'")
    if not isinstance(payload["priority"], list) or not payload["priority"]:
        raise ValueError(f"{payload['gene']}: empty priority list")
    for row in payload["priority"]:
        missing = {"position", "residue", "score", "rationale"} - row.keys()
        if missing:
            raise ValueError(f"{payload['gene']}: priority row missing {missing}")

    prov = payload["provenance"]
    if not isinstance(prov, dict):
        raise ValueError(f"{payload['gene']}: provenance must be an object")
    prov_missing = {
        "tool_version", "git_commit", "packages", "alphafold_model", "fetched", "weights"
    } - prov.keys()
    if prov_missing:
        raise ValueError(f"{payload['gene']}: provenance missing {prov_missing}")
    pkg_missing = set(_PROVENANCE_PACKAGES) - (prov.get("packages") or {}).keys()
    if pkg_missing:
        raise ValueError(f"{payload['gene']}: provenance.packages missing {pkg_missing}")

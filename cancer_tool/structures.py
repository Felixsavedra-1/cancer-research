"""Fetch a finished AlphaFold structure as PDB text (never runs AlphaFold).

pLDDT confidence is carried in the B-factor column. The model file version bumps
over time (v4 → v6 → …), so the URL is resolved via the API, not hardcoded.
"""

from __future__ import annotations

import requests

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb"


def resolve_pdb_url(accession: str, session: requests.Session | None = None) -> str:
    """Resolve the current AlphaFold model URL via the API, falling back to the
    versioned file URL if the API is unreachable."""
    http = session or requests
    try:
        meta = http.get(ALPHAFOLD_API.format(accession=accession), timeout=30)
        if meta.ok and meta.json():
            pdb_url = meta.json()[0].get("pdbUrl")
            if pdb_url:
                return pdb_url
    except (requests.RequestException, ValueError, KeyError, IndexError):
        pass
    return ALPHAFOLD_FILE.format(accession=accession)


def fetch_alphafold_pdb(accession: str, session: requests.Session | None = None) -> str:
    """Download the AlphaFold model for an accession as PDB text (pLDDT in the
    B-factor column). Raises on HTTP failure."""
    http = session or requests
    pdb = http.get(resolve_pdb_url(accession, session), timeout=60)
    pdb.raise_for_status()
    return pdb.text

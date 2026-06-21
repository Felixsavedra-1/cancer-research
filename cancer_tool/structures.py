"""Fetch precomputed AlphaFold structures as PDB text."""

from __future__ import annotations

import requests

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb"


def resolve_pdb_url(accession: str, session: requests.Session | None = None) -> str:
    """Return the current AlphaFold model PDB URL for an accession.

    Asks the API (the file version bumps over time — v4 → v6 → …) and falls back to
    the direct versioned file URL if the API is unavailable.
    """
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
    """Return AlphaFold model PDB text for a UniProt accession.

    Raises ``requests.HTTPError`` if the structure cannot be retrieved.
    """
    http = session or requests
    pdb = http.get(resolve_pdb_url(accession, session), timeout=60)
    pdb.raise_for_status()
    return pdb.text

"""Gene symbol → UniProt accession, name, and canonical sequence.

Queries UniProtKB for the single reviewed (Swiss-Prot) human entry matching an
exact gene symbol.
"""

from __future__ import annotations

import requests

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def get_protein(gene: str, session: requests.Session | None = None) -> dict | None:
    """Resolve a human gene symbol to its reviewed UniProt record.

    Returns ``{accession, name, gene, sequence, length}`` or ``None`` if no
    reviewed human entry matches. Raises on HTTP failure (caller decides).
    """
    http = session or requests
    params = {
        "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,protein_name,sequence,gene_names",
        "format": "json",
        "size": 1,
    }
    resp = http.get(UNIPROT_SEARCH, params=params, timeout=30)
    resp.raise_for_status()
    results = resp.json().get("results", [])
    if not results:
        return None

    rec = results[0]
    accession = rec["primaryAccession"]
    name = (
        rec.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value")
        or accession
    )
    sequence = rec.get("sequence", {}).get("value", "")
    return {
        "accession": accession,
        "name": name,
        "gene": gene.upper(),
        "sequence": sequence,
        "length": len(sequence),
    }

"""End-to-end pipeline reproducibility, network-free.

Replays real captured KRAS API responses (tests/fixtures/kras/) through the
*actual* ``pipeline.analyze_gene`` — UniProt parse → AlphaFold structure → ENM
dynamics → AlphaMissense → LIGSITE pockets → composite scoring — and asserts the
ranking re-derives the known driver (KRAS → G12) from first principles.

Unlike tests/test_validation.py (which checks the committed data/*.json artifacts
still contain the expected drivers), this recomputes the ranking, so it catches a
regression anywhere in the fetch/parse/science path. It uses a fake session, so
it stays in the network-free suite.
"""

import json
from pathlib import Path

import pytest
import requests

from cancer_tool import pipeline

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "kras"


class _FakeResponse:
    def __init__(self, text: str, ok: bool = True):
        self.text = text
        self.ok = ok
        self.status_code = 200 if ok else 500

    def json(self):
        return json.loads(self.text)

    def raise_for_status(self):
        if not self.ok:
            raise requests.HTTPError(f"status {self.status_code}")


class _FakeSession:
    """Serves the captured KRAS fixtures by URL, so no network is touched."""

    def __init__(self):
        self._uniprot = (FIXTURES / "uniprot_search.json").read_text()
        self._af = (FIXTURES / "af_prediction.json").read_text()
        self._pdb = (FIXTURES / "model.pdb").read_text()
        self._am = (FIXTURES / "alphamissense.csv").read_text()
        self._hotspots = (FIXTURES / "hotspots_kras.json").read_text()

    def get(self, url, **kwargs):
        if "rest.uniprot.org" in url:
            return _FakeResponse(self._uniprot)
        if "/api/prediction/" in url:  # AlphaFold prediction metadata
            return _FakeResponse(self._af)
        if url.endswith(".pdb"):
            return _FakeResponse(self._pdb)
        if url.endswith(".csv"):
            return _FakeResponse(self._am)
        if "cancerhotspots.org" in url:
            return _FakeResponse(self._hotspots)
        raise AssertionError(f"unexpected GET in offline test: {url}")

    def post(self, url, **kwargs):
        # Open Targets enrichment is best-effort and unused by scoring; fail it.
        raise requests.ConnectionError("offline")


@pytest.fixture(scope="module")
def kras_payload():
    return pipeline.analyze_gene("KRAS", session=_FakeSession())


def test_pipeline_reproduces_kras_g12_from_first_principles(kras_payload):
    top = kras_payload["priority"][0]
    assert top["position"] == 12, f"expected G12 to top KRAS, got {top['residue']}"


def test_offline_payload_passes_schema(kras_payload):
    pipeline.validate_payload(kras_payload)


def test_offline_pipeline_runs_full_science(kras_payload):
    # Dynamics and pockets really ran (not degraded to None/empty), proving the
    # whole scientific path executed on the captured structure.
    assert kras_payload["dynamics"] is not None
    assert kras_payload["dynamics"]["collectivity"] > 0
    top = kras_payload["priority"][0]
    assert top["pathogenicity"] > 0  # AlphaMissense joined
    assert top["rationale"]

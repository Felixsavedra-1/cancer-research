import json

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from api import main as api_main  # noqa: E402


@pytest.fixture()
def client(tmp_path, monkeypatch):
    # Redirect the committed-data and cache dirs at a temp location so the test
    # never reads/writes the real repo data.
    monkeypatch.setattr(api_main, "DATA_DIR", tmp_path)
    monkeypatch.setattr(api_main, "CACHE_DIR", tmp_path / "cache")
    return TestClient(api_main.app)


def _fake_payload(gene: str) -> dict:
    return {
        "gene": gene,
        "accession": "P00000",
        "name": "Fake protein",
        "length": 3,
        "provenance": {
            "tool_version": "test",
            "git_commit": "abc1234",
            "packages": {"prody": "0", "numpy": "0", "scipy": "0"},
            "alphafold_model": "v6",
            "fetched": "2026-01-01",
            "weights": {"recurrence": 0.3, "pathogenicity": 0.35, "druggability": 0.2, "criticality": 0.15},
        },
        "priority": [{"position": 1, "residue": "M1", "score": 42.0, "rationale": "test"}],
    }


def test_health_reports_ok(client):
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_rejects_non_alphanumeric_gene(client):
    assert client.get("/analyze/TP53;DROP").status_code == 400


def test_serves_committed_file_without_computing(client, tmp_path, monkeypatch):
    (tmp_path / "TP53.json").write_text(json.dumps(_fake_payload("TP53")))

    def _boom(*a, **k):  # must not be called when a committed file exists
        raise AssertionError("pipeline should not run for a committed gene")

    monkeypatch.setattr(api_main.pipeline, "analyze_gene", _boom)
    r = client.get("/analyze/TP53")
    assert r.status_code == 200
    assert r.json()["gene"] == "TP53"


def test_computes_caches_then_reuses(client, tmp_path, monkeypatch):
    calls = {"n": 0}

    def _fake_analyze(gene, session=None, top_n=30):
        calls["n"] += 1
        return _fake_payload(gene)

    monkeypatch.setattr(api_main.pipeline, "analyze_gene", _fake_analyze)
    first = client.get("/analyze/pik3ca")
    assert first.status_code == 200
    assert first.json()["gene"] == "PIK3CA"
    assert (tmp_path / "cache" / "PIK3CA.json").exists()
    # Second call is served from cache — pipeline runs only once.
    client.get("/analyze/PIK3CA")
    assert calls["n"] == 1


def test_unknown_gene_returns_404(client, monkeypatch):
    monkeypatch.setattr(api_main.pipeline, "analyze_gene", lambda *a, **k: None)
    assert client.get("/analyze/NOTAGENE").status_code == 404

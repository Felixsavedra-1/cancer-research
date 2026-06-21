"""The rankings must rediscover known cancer biology. These guard the committed
``data/*.json`` against scoring/parsing regressions — network-free, run in CI.

If a change breaks these, the science is wrong, not the test.
"""

import json
from pathlib import Path

import pytest

from scripts.precompute import validate_payload

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

# Textbook drivers each gene's ranking should surface near the top.
KNOWN_DRIVERS = {
    "TP53": ({"R273", "R248", "R175"}, 5),
    "KRAS": ({"G12"}, 1),
    "BRAF": ({"V600"}, 1),
    "EGFR": ({"L858", "T790"}, 5),
    "IDH2": ({"R172", "R140"}, 2),
    # KIT has several co-equal activating drivers; D816 sits in a low-pLDDT region so its
    # structural-criticality bonus is (correctly) down-weighted, but it stays a top driver.
    "KIT": ({"D816"}, 3),
}


def _load(gene: str) -> dict:
    return json.loads((DATA_DIR / f"{gene}.json").read_text())


@pytest.mark.parametrize("gene,expected", KNOWN_DRIVERS.items())
def test_known_drivers_rank_at_top(gene, expected):
    drivers, top_n = expected
    top = {row["residue"] for row in _load(gene)["priority"][:top_n]}
    assert drivers <= top, f"{gene}: expected {drivers} within top {top_n}, got {top}"


@pytest.mark.parametrize("path", sorted(DATA_DIR.glob("*.json")))
def test_committed_data_matches_schema(path):
    if path.name == "index.json":
        return
    validate_payload(json.loads(path.read_text()))

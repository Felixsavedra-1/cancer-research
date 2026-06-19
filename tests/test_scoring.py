"""Unit tests for the Target Priority composite score (pure, network-free)."""

from cancer_tool import scoring

# A short sequence where position 5 is R and position 8 is G (1-indexed).
SEQUENCE = "MKTLR" "AAG" "AA"  # positions: 1 M,2 K,3 T,4 L,5 R,6 A,7 A,8 G,9 A,10 A


def _hotspots():
    return [
        {"residue": "R5", "position": 5, "count": 100, "variants": {"H": 90, "C": 10}},
        {"residue": "G8", "position": 8, "count": 10, "variants": {"D": 10}},
    ]


def _pathogenicity():
    return {
        "by_position": {5: {"mean": 0.5, "max": 0.9, "n": 19}, 8: {"mean": 0.2, "max": 0.3, "n": 19}},
        "by_variant": {
            "R5H": {"score": 0.95, "class": "LPath", "class_label": "Likely pathogenic"},
            "G8D": {"score": 0.20, "class": "LBen", "class_label": "Likely benign"},
        },
    }


def _dynamics():
    # Position 5 rigid (low flex), position 8 flexible; hinge near 8.
    return {
        "residue_numbers": [5, 8],
        "flexibility": [0.1, 0.9],
        "rigidity": [0.9, 0.1],
        "hinges": [8],
        "collective_motion": [0.2, 0.8],
        "n_modes": 5,
    }


def _pockets():
    return [{"residues": [5, 6, 7], "druggability": 0.8, "volume": 300.0, "center": [0, 0, 0]}]


def test_ranks_driver_above_passenger():
    rows = scoring.score_residues(
        _hotspots(), _pathogenicity(), _dynamics(), _pockets(), SEQUENCE
    )
    assert [r["position"] for r in rows] == [5, 8]
    assert rows[0]["score"] > rows[1]["score"]


def test_uses_exact_variant_pathogenicity():
    rows = scoring.score_residues(
        _hotspots(), _pathogenicity(), _dynamics(), _pockets(), SEQUENCE
    )
    top = rows[0]
    # R5H exact score (0.95) is used, not the position mean (0.5).
    assert top["pathogenicity"] == 0.95
    assert top["am_class"] == "LPath"
    assert top["wt"] == "R"


def test_score_formula_matches_weights():
    rows = scoring.score_residues(
        _hotspots(), _pathogenicity(), _dynamics(), _pockets(), SEQUENCE
    )
    top = rows[0]  # position 5
    w = scoring.DEFAULT_WEIGHTS
    # recurrence 1.0, pathogenicity 0.95, druggability 0.8, criticality=0.7*0.9=0.63
    expected = 100 * (
        w["recurrence"] * 1.0
        + w["pathogenicity"] * 0.95
        + w["druggability"] * 0.8
        + w["criticality"] * 0.63
    )
    assert abs(top["score"] - round(expected, 1)) < 0.2


def test_rationale_is_explained():
    rows = scoring.score_residues(
        _hotspots(), _pathogenicity(), _dynamics(), _pockets(), SEQUENCE
    )
    assert "recurrently mutated" in rows[0]["rationale"]
    assert "druggable pocket" in rows[0]["rationale"]


def test_degrades_without_optional_signals():
    # Only hotspots available — no AM, dynamics, or pockets.
    rows = scoring.score_residues(_hotspots(), None, None, None, SEQUENCE)
    assert rows[0]["position"] == 5  # recurrence alone still ranks it first
    assert rows[0]["pathogenicity"] == 0.0
    assert rows[0]["druggability"] == 0.0
    assert rows[0]["criticality"] == 0.0


def test_empty_hotspots_returns_empty():
    assert scoring.score_residues([], _pathogenicity(), _dynamics(), _pockets(), SEQUENCE) == []


def test_custom_weights_override():
    weights = {"recurrence": 1.0, "pathogenicity": 0.0, "druggability": 0.0, "criticality": 0.0}
    rows = scoring.score_residues(
        _hotspots(), _pathogenicity(), _dynamics(), _pockets(), SEQUENCE, weights=weights
    )
    assert rows[0]["score"] == 100.0  # recurrence 1.0 * 100

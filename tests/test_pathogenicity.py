from cancer_tool import pathogenicity as patho

SAMPLE_CSV = """protein_variant,am_pathogenicity,am_class
M1A,0.80,LPath
M1C,0.60,LPath
R175H,0.99,LPath
R175G,0.10,LBen
junkrow
G2,bad,Amb
"""


def test_parse_builds_variant_index():
    parsed = patho.parse_alphamissense_csv(SAMPLE_CSV)
    assert parsed["by_variant"]["R175H"] == {
        "score": 0.99,
        "class": "LPath",
        "class_label": "Likely pathogenic",
    }


def test_parse_aggregates_by_position():
    parsed = patho.parse_alphamissense_csv(SAMPLE_CSV)
    pos1 = parsed["by_position"][1]
    assert pos1["n"] == 2
    assert pos1["mean"] == 0.70
    assert pos1["max"] == 0.80


def test_parse_skips_malformed_rows():
    parsed = patho.parse_alphamissense_csv(SAMPLE_CSV)
    assert "junk" not in str(parsed["by_variant"])
    assert 2 not in parsed["by_position"]


def test_parse_tolerates_missing_header():
    headerless = "M1A,0.80,LPath\nR175H,0.99,LPath\n"
    parsed = patho.parse_alphamissense_csv(headerless)
    assert "M1A" in parsed["by_variant"]
    assert "R175H" in parsed["by_variant"]


def test_lookup_helpers():
    parsed = patho.parse_alphamissense_csv(SAMPLE_CSV)
    assert patho.variant_score(parsed, "r175h")["score"] == 0.99
    assert patho.variant_score(parsed, "X999Y") is None
    assert patho.position_pathogenicity(parsed, 1) == 0.70
    assert patho.position_pathogenicity(parsed, 9999) == 0.0

from cancer_tool import mutations as mut

TP53_FRAGMENT_170 = "TEVVRRCPHHE"


def _tp53_like_sequence() -> str:
    return "A" * 169 + TP53_FRAGMENT_170 + "A" * 200


def test_parse_basic_mutation():
    assert mut.parse_mutation("R175H") == {
        "wt": "R",
        "position": 175,
        "mut": "H",
        "label": "R175H",
    }


def test_parse_is_case_insensitive_and_trims():
    assert mut.parse_mutation("  r175h ")["label"] == "R175H"


def test_parse_nonsense_returns_none():
    assert mut.parse_mutation("hello") is None
    assert mut.parse_mutation("175") is None
    assert mut.parse_mutation("") is None


def test_parse_multiple_separators():
    parsed, bad = mut.parse_mutations("R175H, G245S nonsense\nR248Q")
    assert [m["label"] for m in parsed] == ["R175H", "G245S", "R248Q"]
    assert bad == ["nonsense"]


def test_validate_matches_wild_type():
    seq = _tp53_like_sequence()
    assert seq[174] == "R"
    ok, _ = mut.validate_mutation(mut.parse_mutation("R175H"), seq)
    assert ok


def test_validate_rejects_wrong_wild_type():
    seq = _tp53_like_sequence()
    ok, message = mut.validate_mutation(mut.parse_mutation("A175H"), seq)
    assert not ok
    assert "175" in message


def test_validate_rejects_out_of_range():
    ok, message = mut.validate_mutation(mut.parse_mutation("R9999H"), "MRT")
    assert not ok
    assert "outside" in message

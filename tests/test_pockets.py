import math

import numpy as np

from cancer_tool import pockets


def test_shift_moves_voxels_with_zero_fill():
    arr = np.zeros((3, 1, 1), dtype=bool)
    arr[0, 0, 0] = True
    shifted = pockets._shift(arr, (1, 0, 0))
    assert shifted[1, 0, 0]
    assert not shifted[0, 0, 0] and not shifted[2, 0, 0]


def test_ray_presence_sees_protein_behind():
    protein = np.zeros((3, 1, 1), dtype=bool)
    protein[0, 0, 0] = True
    acc = pockets._ray_presence(protein, (1, 0, 0), steps=3)
    assert list(acc[:, 0, 0]) == [False, True, True]


def test_pocket_proximity_returns_best_lining_druggability():
    found = [
        {"residues": [5, 6, 7], "druggability": 0.8},
        {"residues": [7, 8], "druggability": 0.3},
    ]
    assert pockets.pocket_proximity(6, found) == 0.8
    assert pockets.pocket_proximity(7, found) == 0.8
    assert pockets.pocket_proximity(99, found) == 0.0
    assert pockets.pocket_proximity(5, []) == 0.0


def _hollow_sphere_pdb(n_points: int = 260, radius: float = 8.0) -> str:
    golden = math.pi * (3.0 - math.sqrt(5.0))
    lines = []
    for i in range(n_points):
        y = 1.0 - 2.0 * (i + 0.5) / n_points
        r = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden * i
        x, z = math.cos(theta) * r, math.sin(theta) * r
        lines.append(
            f"ATOM  {i + 1:>5}  CA  ALA A{i + 1:>4}    "
            f"{x * radius:8.3f}{y * radius:8.3f}{z * radius:8.3f}  1.00 50.00           C"
        )
    lines.append("END")
    return "\n".join(lines)


def test_ligsite_finds_enclosed_cavity():
    found = pockets._detect_ligsite(_hollow_sphere_pdb())
    assert found, "the hollow centre should register as a buried pocket"
    top = found[0]
    assert top["source"] == "ligsite"
    assert top["volume"] > pockets.MIN_POCKET_VOLUME
    assert 0.0 < top["druggability"] <= 1.0
    assert top["residues"]

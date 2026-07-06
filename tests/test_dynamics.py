import math

import numpy as np

from cancer_tool import dynamics


def _alpha_helix_pdb(n: int = 16, bfactor: float = 50.0) -> str:
    radius, rise, turn = 2.3, 1.5, math.radians(100.0)
    lines = []
    for i in range(n):
        x = radius * math.cos(i * turn)
        y = radius * math.sin(i * turn)
        z = i * rise
        lines.append(
            f"ATOM  {i + 1:>5}  CA  ALA A{i + 1:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{bfactor:6.2f}           C"
        )
    lines.append("END")
    return "\n".join(lines)


def test_returns_per_residue_arrays():
    result = dynamics.compute_dynamics(_alpha_helix_pdb(16))
    assert len(result["flexibility"]) == 16
    assert len(result["rigidity"]) == 16
    assert len(result["mode_amplitude"]) == 16
    assert result["n_modes"] >= 1


def test_records_enm_parameters_for_provenance():
    result = dynamics.compute_dynamics(_alpha_helix_pdb(16), gnm_cutoff=9.0, gamma=1.0)
    assert result["params"]["gnm_cutoff"] == 9.0
    assert result["params"]["anm_cutoff"] == dynamics.DEFAULT_ANM_CUTOFF
    assert 0.0 <= result["collectivity"] <= 1.0


def test_collectivity_high_for_global_mode():
    # A uniform whole-body amplitude mobilises every residue → collectivity ~1.
    uniform = np.ones(20)
    assert dynamics._degree_of_collectivity(uniform) > 0.99
    # A single mobile residue → collectivity ~1/N.
    localised = np.zeros(20)
    localised[0] = 1.0
    assert dynamics._degree_of_collectivity(localised) < 0.1


def test_flexibility_is_normalised():
    flex = dynamics.compute_dynamics(_alpha_helix_pdb(16))["flexibility"]
    assert min(flex) >= 0.0
    assert max(flex) <= 1.0
    assert math.isclose(max(flex), 1.0, abs_tol=1e-6)


def test_termini_more_flexible_than_core():
    flex = dynamics.compute_dynamics(_alpha_helix_pdb(20))["flexibility"]
    core_mean = float(np.mean(flex[8:12]))
    termini_mean = float(np.mean([flex[0], flex[1], flex[-2], flex[-1]]))
    assert termini_mean > core_mean


def test_rigidity_is_complement_of_flexibility():
    result = dynamics.compute_dynamics(_alpha_helix_pdb(16))
    for flex, rigid in zip(result["flexibility"], result["rigidity"]):
        assert math.isclose(flex + rigid, 1.0, abs_tol=1e-3)


def test_is_deterministic():
    pdb = _alpha_helix_pdb(16)
    assert dynamics.compute_dynamics(pdb)["flexibility"] == dynamics.compute_dynamics(pdb)["flexibility"]


def test_lookup_helpers_align_positions():
    result = dynamics.compute_dynamics(_alpha_helix_pdb(12))
    flex_map = dynamics.flexibility_by_position(result)
    assert set(flex_map) == set(result["residue_numbers"])
    rigid_map = dynamics.rigidity_by_position(result)
    assert math.isclose(flex_map[6] + rigid_map[6], 1.0, abs_tol=1e-3)


def test_plddt_is_read_from_b_factors():
    result = dynamics.compute_dynamics(_alpha_helix_pdb(12))
    assert len(result["plddt"]) == 12
    assert all(b == 50.0 for b in result["plddt"])
    plddt_map = dynamics.plddt_by_position(result)
    assert set(plddt_map) == set(result["residue_numbers"])


def test_plddt_is_confidence_flag_for_alphafold_range():
    # AlphaFold pLDDT lives in [0, 100] → treated as confidence.
    result = dynamics.compute_dynamics(_alpha_helix_pdb(12, bfactor=80.0))
    assert result["plddt_is_confidence"] is True


def test_plddt_is_confidence_false_for_experimental_bfactors():
    # A real crystal structure stores B-factors that can exceed 100 → NOT pLDDT,
    # so the flag must be False and scoring must not read them as confidence.
    result = dynamics.compute_dynamics(_alpha_helix_pdb(12, bfactor=120.0))
    assert result["plddt_is_confidence"] is False

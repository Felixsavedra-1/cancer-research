"""Folding dynamics via elastic-network normal-mode analysis (ProDy).

Recovers the low-frequency collective motions an all-atom MD run would sample —
per-residue flexibility and hinge sites (GNM), plus the slowest collective mode
(ANM) — in milliseconds on a CPU from a single structure. See docs/METHODS.md
for parameters and references.
"""

from __future__ import annotations

import io

import numpy as np

try:
    import prody

    prody.confProDy(verbosity="none")
    from prody import ANM, GNM, calcSqFlucts, parsePDBStream
except Exception:
    prody = None


DEFAULT_MODES = 10

# Standard elastic-network cutoffs and a uniform spring constant, kept explicit
# for reproducibility rather than left to ProDy's defaults.
DEFAULT_GNM_CUTOFF = 10.0
DEFAULT_ANM_CUTOFF = 15.0
DEFAULT_GAMMA = 1.0
DEFAULT_HINGE_MODES = 3  # slowest GNM modes pooled for hinge zero-crossings


class DynamicsError(RuntimeError):
    pass


def _parse_calphas(pdb_text: str):
    if prody is None:
        raise DynamicsError(
            "ProDy is not installed. Install it with `pip install prody`."
        )
    structure = parsePDBStream(io.StringIO(pdb_text))
    if structure is None:
        raise DynamicsError("Could not parse the PDB text.")
    calphas = structure.select("calpha")
    if calphas is None or calphas.numAtoms() < 4:
        raise DynamicsError("Too few Cα atoms for elastic network analysis.")
    return calphas


def _zero_crossings(vector: np.ndarray) -> list[int]:
    signs = np.sign(vector)
    for i in range(1, len(signs)):
        if signs[i] == 0:
            signs[i] = signs[i - 1]
    return [i for i in range(1, len(signs)) if signs[i] != signs[i - 1]]


def _normalize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return values
    lo, hi = float(values.min()), float(values.max())
    if hi - lo < 1e-12:
        return np.zeros_like(values)
    return (values - lo) / (hi - lo)


def _degree_of_collectivity(amplitudes: np.ndarray) -> float:
    """Fraction of residues significantly mobilised by a mode, from per-residue
    squared displacement amplitudes: ~1 is a global motion, near 0 a localised one."""
    a2 = np.asarray(amplitudes, dtype=float) ** 2
    total = a2.sum()
    if total <= 0 or a2.size == 0:
        return 0.0
    a2 = a2 / total
    nonzero = a2[a2 > 0]
    entropy = -np.sum(nonzero * np.log(nonzero))
    return round(float(np.exp(entropy) / a2.size), 4)


def compute_dynamics(
    pdb_text: str,
    n_modes: int = DEFAULT_MODES,
    gnm_cutoff: float = DEFAULT_GNM_CUTOFF,
    anm_cutoff: float = DEFAULT_ANM_CUTOFF,
    gamma: float = DEFAULT_GAMMA,
    hinge_modes: int = DEFAULT_HINGE_MODES,
) -> dict:
    """Compute folding dynamics from a single structure via elastic-network NMA.

    Returns per-residue flexibility/rigidity, hinge sites, the slowest ANM mode's
    collectivity, and pLDDT (read from the Cα B-factor column). ``plddt_is_confidence``
    flags whether that column is really AlphaFold pLDDT (all values in [0, 100]) vs
    experimental B-factors, so scoring can avoid misreading crystallographic B-factors
    as confidence. Cutoffs and gamma are echoed into ``params`` for provenance.
    flexibility/rigidity are min-max normalised per protein — comparable within a
    structure, not across.
    """
    calphas = _parse_calphas(pdb_text)
    n_atoms = calphas.numAtoms()
    resnums = [int(n) for n in calphas.getResnums()]
    betas = [float(b) for b in calphas.getBetas()]
    plddt = [round(b, 1) for b in betas]
    plddt_is_confidence = bool(betas) and all(0.0 <= b <= 100.0 for b in betas)
    n_modes = max(1, min(n_modes, n_atoms - 1))

    gnm = GNM("enm")
    gnm.buildKirchhoff(calphas, cutoff=gnm_cutoff, gamma=gamma)
    gnm.calcModes(n_modes=n_modes)

    sqflucts = np.asarray(calcSqFlucts(gnm), dtype=float)
    flexibility = _normalize(sqflucts)

    # Pool hinges over the slowest few modes, not just the slowest, which alone
    # misses secondary pivots.
    n_hinge = max(1, min(hinge_modes, gnm.numModes()))
    hinge_set: set[int] = set()
    for m in range(n_hinge):
        vec = np.asarray(gnm[m].getEigvec(), dtype=float).ravel()
        for i in _zero_crossings(vec):
            if 0 <= i < len(resnums):
                hinge_set.add(resnums[i])
    hinges = sorted(hinge_set)

    anm = ANM("enm")
    anm.buildHessian(calphas, cutoff=anm_cutoff, gamma=gamma)
    anm.calcModes(n_modes=min(n_modes, 3 * n_atoms - 6))
    slow_anm = np.asarray(anm[0].getEigvec(), dtype=float).reshape(-1, 3)
    amplitude = np.linalg.norm(slow_anm, axis=1)
    collectivity = _degree_of_collectivity(amplitude)

    return {
        "residue_numbers": resnums,
        "plddt": plddt,
        "plddt_is_confidence": plddt_is_confidence,
        "flexibility": [round(float(x), 4) for x in flexibility],
        "rigidity": [round(float(1.0 - x), 4) for x in flexibility],
        "hinges": hinges,
        "mode_amplitude": [round(float(x), 4) for x in _normalize(amplitude)],
        "collectivity": collectivity,
        "n_modes": int(gnm.numModes()),
        "params": {
            "gnm_cutoff": gnm_cutoff,
            "anm_cutoff": anm_cutoff,
            "gamma": gamma,
            "n_modes_requested": n_modes,
            "hinge_modes": n_hinge,
        },
    }


def flexibility_by_position(dynamics: dict) -> dict[int, float]:
    return dict(zip(dynamics["residue_numbers"], dynamics["flexibility"]))


def rigidity_by_position(dynamics: dict) -> dict[int, float]:
    return dict(zip(dynamics["residue_numbers"], dynamics["rigidity"]))


def plddt_by_position(dynamics: dict) -> dict[int, float]:
    return dict(zip(dynamics["residue_numbers"], dynamics.get("plddt", [])))

"""Intrinsic protein dynamics from one structure via Elastic Network Model analysis.

GNM square fluctuations give per-residue flexibility and hinge sites; ANM gives the
slowest collective ("breathing") mode — the laptop-grade alternative to molecular
dynamics. References: Bahar et al. 1997; Atilgan et al. 2001; ProDy (Bakan et al. 2011).
"""

from __future__ import annotations

import io

import numpy as np

# ProDy is chatty by default; keep it quiet so it doesn't pollute app/CLI output.
try:  # pragma: no cover - import-time configuration
    import prody

    prody.confProDy(verbosity="none")
    from prody import ANM, GNM, calcSqFlucts, parsePDBStream
except Exception:  # pragma: no cover - surfaced at call time instead
    prody = None


# Number of low-frequency modes to retain. The slowest handful carry the
# functionally relevant collective motions; more just adds high-frequency noise.
DEFAULT_MODES = 10


class DynamicsError(RuntimeError):
    """Raised when ENM analysis cannot be performed on the given structure."""


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
    """Indices where a mode shape crosses zero — the network's hinge points."""
    signs = np.sign(vector)
    # Treat exact zeros as continuations of the previous sign so a single zero
    # residue isn't double-counted as two crossings.
    for i in range(1, len(signs)):
        if signs[i] == 0:
            signs[i] = signs[i - 1]
    return [i for i in range(1, len(signs)) if signs[i] != signs[i - 1]]


def _normalize(values: np.ndarray) -> np.ndarray:
    """Min–max scale to 0–1; returns all-zeros for a flat/empty input."""
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return values
    lo, hi = float(values.min()), float(values.max())
    if hi - lo < 1e-12:
        return np.zeros_like(values)
    return (values - lo) / (hi - lo)


def compute_dynamics(pdb_text: str, n_modes: int = DEFAULT_MODES) -> dict:
    """Run ENM/NMA on a structure and summarise its intrinsic dynamics.

    Returns (all per-residue lists aligned to ``residue_numbers``)::

        {"residue_numbers": [int, ...],
         "plddt": [0-100, ...],            # AlphaFold per-residue confidence
         "flexibility": [0-1, ...],        # GNM square fluctuations, 0 rigid → 1 mobile
         "rigidity": [0-1, ...],           # 1 - flexibility
         "hinges": [int, ...],             # slowest-mode zero-crossings (domain pivots)
         "collective_motion": [0-1, ...],  # slowest ANM mode magnitude
         "n_modes": int}

    Deterministic. Raises :class:`DynamicsError` if the structure can't be analysed.
    """
    calphas = _parse_calphas(pdb_text)
    n_atoms = calphas.numAtoms()
    resnums = [int(n) for n in calphas.getResnums()]
    # AlphaFold stores per-residue pLDDT confidence (0–100) in the B-factor column.
    plddt = [round(float(b), 1) for b in calphas.getBetas()]
    # Can't request more non-trivial modes than the system supports.
    n_modes = max(1, min(n_modes, n_atoms - 1))

    gnm = GNM("enm")
    gnm.buildKirchhoff(calphas)
    gnm.calcModes(n_modes=n_modes)

    sqflucts = np.asarray(calcSqFlucts(gnm), dtype=float)
    flexibility = _normalize(sqflucts)

    # Slowest non-trivial GNM mode → hinge sites (domain boundaries).
    slowest = np.asarray(gnm[0].getEigvec(), dtype=float).ravel()
    hinges = [resnums[i] for i in _zero_crossings(slowest) if 0 <= i < len(resnums)]

    # ANM gives 3D mode vectors; the slowest one's per-residue magnitude is the
    # dominant collective ("breathing") motion.
    anm = ANM("enm")
    anm.buildHessian(calphas)
    anm.calcModes(n_modes=min(n_modes, 3 * n_atoms - 6))
    slow_anm = np.asarray(anm[0].getEigvec(), dtype=float).reshape(-1, 3)
    collective = _normalize(np.linalg.norm(slow_anm, axis=1))

    return {
        "residue_numbers": resnums,
        "plddt": plddt,
        "flexibility": [round(float(x), 4) for x in flexibility],
        "rigidity": [round(float(1.0 - x), 4) for x in flexibility],
        "hinges": hinges,
        "collective_motion": [round(float(x), 4) for x in collective],
        "n_modes": int(gnm.numModes()),
    }


def flexibility_by_position(dynamics: dict) -> dict[int, float]:
    """Map residue number → normalised flexibility for O(1) lookups during scoring."""
    return dict(zip(dynamics["residue_numbers"], dynamics["flexibility"]))


def rigidity_by_position(dynamics: dict) -> dict[int, float]:
    """Map residue number → normalised rigidity (1 − flexibility)."""
    return dict(zip(dynamics["residue_numbers"], dynamics["rigidity"]))


def plddt_by_position(dynamics: dict) -> dict[int, float]:
    """Map residue number → AlphaFold pLDDT confidence (0–100)."""
    return dict(zip(dynamics["residue_numbers"], dynamics.get("plddt", [])))

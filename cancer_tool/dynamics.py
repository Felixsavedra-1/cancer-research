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


def compute_dynamics(pdb_text: str, n_modes: int = DEFAULT_MODES) -> dict:
    calphas = _parse_calphas(pdb_text)
    n_atoms = calphas.numAtoms()
    resnums = [int(n) for n in calphas.getResnums()]
    plddt = [round(float(b), 1) for b in calphas.getBetas()]
    n_modes = max(1, min(n_modes, n_atoms - 1))

    gnm = GNM("enm")
    gnm.buildKirchhoff(calphas)
    gnm.calcModes(n_modes=n_modes)

    sqflucts = np.asarray(calcSqFlucts(gnm), dtype=float)
    flexibility = _normalize(sqflucts)

    slowest = np.asarray(gnm[0].getEigvec(), dtype=float).ravel()
    hinges = [resnums[i] for i in _zero_crossings(slowest) if 0 <= i < len(resnums)]

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
    return dict(zip(dynamics["residue_numbers"], dynamics["flexibility"]))


def rigidity_by_position(dynamics: dict) -> dict[int, float]:
    return dict(zip(dynamics["residue_numbers"], dynamics["rigidity"]))


def plddt_by_position(dynamics: dict) -> dict[int, float]:
    return dict(zip(dynamics["residue_numbers"], dynamics.get("plddt", [])))

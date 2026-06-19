"""Druggable-pocket detection via a LIGSITE-style geometric scan (Hendlich et al. 1997).

Grid the protein, flag solvent voxels enclosed between protein along several scan lines,
cluster them into pockets, rank by enclosed volume. Uses the `fpocket` binary when present,
else the always-available pure-Python scan. No GPU, no docking.
"""

from __future__ import annotations

import io
import shutil

import numpy as np

try:  # pragma: no cover - import-time configuration
    import prody

    prody.confProDy(verbosity="none")
    from prody import parsePDBStream
except Exception:  # pragma: no cover
    prody = None

# Grid scan parameters. Spacing trades resolution for speed; 1.0 Å is plenty for
# laptop-scale proteins. A voxel counts as "protein" within this radius of a heavy
# atom (~heavy-atom vdW + small probe).
GRID_SPACING = 1.0
PROTEIN_RADIUS = 3.0
# Min protein-solvent-protein scan lines (of 7) for a voxel to be "buried".
MIN_PSP = 5
# Smallest pocket worth reporting, in Å³.
MIN_POCKET_VOLUME = 80.0
# Volume (Å³) at which a pocket is considered fully druggable; smaller scales down.
DRUGGABLE_VOLUME_REF = 350.0
# A residue lines a pocket if a heavy atom sits within this of a pocket voxel.
LINING_DISTANCE = 4.5

# Seven scan-line directions (3 axes + 4 cube diagonals).
_DIRECTIONS = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 1),
    (1, 1, -1),
    (1, -1, 1),
    (-1, 1, 1),
]


def _shift(arr: np.ndarray, d: tuple[int, int, int]) -> np.ndarray:
    """Shift a 3D boolean array by one voxel along ``d`` with zero (False) fill."""
    out = np.zeros_like(arr)
    src = [slice(None)] * 3
    dst = [slice(None)] * 3
    for axis, step in enumerate(d):
        if step > 0:
            dst[axis], src[axis] = slice(step, None), slice(None, -step)
        elif step < 0:
            dst[axis], src[axis] = slice(None, step), slice(-step, None)
    out[tuple(dst)] = arr[tuple(src)]
    return out


def _ray_presence(protein: np.ndarray, d: tuple[int, int, int], steps: int) -> np.ndarray:
    """For each voxel, is there a protein voxel somewhere along direction ``d``?"""
    acc = np.zeros_like(protein)
    shifted = protein
    for _ in range(steps):
        shifted = _shift(shifted, d)
        acc |= shifted
    return acc


def detect_pockets(pdb_text: str) -> list[dict]:
    """Find candidate druggable pockets in a structure.

    Returns a list (most druggable first) of::

        {"residues": [int, ...],   # residue numbers lining the pocket
         "druggability": 0.0-1.0,  # from enclosed volume
         "volume": float,          # Å³
         "center": [x, y, z],
         "source": "fpocket" | "ligsite"}
    """
    if shutil.which("fpocket"):
        try:
            return _detect_with_fpocket(pdb_text)
        except Exception:
            pass  # fall through to the always-available geometric scan
    return _detect_ligsite(pdb_text)


def _heavy_atoms(pdb_text: str):
    if prody is None:
        raise RuntimeError("ProDy is not installed; cannot detect pockets.")
    structure = parsePDBStream(io.StringIO(pdb_text))
    if structure is None:
        raise RuntimeError("Could not parse the PDB text.")
    heavy = structure.select("heavy") or structure
    return heavy


def _detect_ligsite(pdb_text: str) -> list[dict]:
    from scipy import ndimage
    from scipy.spatial import cKDTree

    heavy = _heavy_atoms(pdb_text)
    coords = np.asarray(heavy.getCoords(), dtype=float)
    resnums = np.asarray(heavy.getResnums())

    # Grid spanning the protein plus a margin so surface pockets aren't clipped.
    margin = PROTEIN_RADIUS + GRID_SPACING
    lo = coords.min(axis=0) - margin
    hi = coords.max(axis=0) + margin
    dims = np.ceil((hi - lo) / GRID_SPACING).astype(int) + 1

    # Classify every voxel as protein (near an atom) or solvent, via KD-tree.
    axes = [lo[i] + np.arange(dims[i]) * GRID_SPACING for i in range(3)]
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    grid_pts = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
    tree = cKDTree(coords)
    near = tree.query_ball_point(grid_pts, r=PROTEIN_RADIUS, return_length=True)
    protein = (np.asarray(near) > 0).reshape(dims)
    solvent = ~protein

    # Count protein-solvent-protein scan lines through each solvent voxel.
    steps = int(max(dims))
    psp = np.zeros(dims, dtype=np.int8)
    for d in _DIRECTIONS:
        fwd = _ray_presence(protein, d, steps)
        bwd = _ray_presence(protein, tuple(-x for x in d), steps)
        psp += (fwd & bwd).astype(np.int8)

    pocket_voxels = solvent & (psp >= MIN_PSP)
    if not pocket_voxels.any():
        return []

    # Cluster contiguous buried voxels into discrete pockets.
    labels, n = ndimage.label(pocket_voxels)
    voxel_volume = GRID_SPACING**3
    pockets: list[dict] = []
    for label in range(1, n + 1):
        mask = labels == label
        count = int(mask.sum())
        volume = count * voxel_volume
        if volume < MIN_POCKET_VOLUME:
            continue
        idx = np.argwhere(mask)
        pts = lo + idx * GRID_SPACING
        center = pts.mean(axis=0)

        # Residues with a heavy atom lining the pocket cavity.
        lining_idx = tree.query_ball_point(pts, r=LINING_DISTANCE)
        lining_atoms = {a for sub in lining_idx for a in sub}
        residues = sorted({int(resnums[a]) for a in lining_atoms})

        druggability = min(1.0, volume / DRUGGABLE_VOLUME_REF)
        pockets.append(
            {
                "residues": residues,
                "druggability": round(float(druggability), 4),
                "volume": round(float(volume), 1),
                "center": [round(float(c), 2) for c in center],
                "source": "ligsite",
            }
        )

    pockets.sort(key=lambda p: p["druggability"], reverse=True)
    return pockets


def _detect_with_fpocket(pdb_text: str) -> list[dict]:  # pragma: no cover - needs binary
    """Run fpocket on the structure and parse its ranked pockets."""
    import os
    import re
    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        pdb_path = os.path.join(tmp, "model.pdb")
        with open(pdb_path, "w") as fh:
            fh.write(pdb_text)
        subprocess.run(
            ["fpocket", "-f", pdb_path],
            check=True,
            capture_output=True,
            timeout=300,
        )
        info_path = os.path.join(tmp, "model_out", "model_info.txt")
        pockets: list[dict] = []
        if not os.path.exists(info_path):
            return []
        with open(info_path) as fh:
            blocks = fh.read().split("Pocket")
        for block in blocks[1:]:
            score = re.search(r"Druggability Score\s*:\s*([\d.]+)", block)
            volume = re.search(r"Volume\s*:\s*([\d.]+)", block)
            if not score:
                continue
            pockets.append(
                {
                    "residues": [],
                    "druggability": round(float(score.group(1)), 4),
                    "volume": round(float(volume.group(1)), 1) if volume else 0.0,
                    "center": [0.0, 0.0, 0.0],
                    "source": "fpocket",
                }
            )
        pockets.sort(key=lambda p: p["druggability"], reverse=True)
        return pockets


def pocket_proximity(position: int, pockets: list[dict]) -> float:
    """Druggability of the best pocket this residue lines, else 0.0."""
    best = 0.0
    for pocket in pockets:
        if position in pocket["residues"]:
            best = max(best, float(pocket["druggability"]))
    return best

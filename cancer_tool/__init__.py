"""Cancer Protein Explorer — AlphaFold structures, protein dynamics, AlphaMissense
pathogenicity and pocket detection fused into a ranked shortlist of druggable
cancer-driver residues.

Research and education use only; not a clinical or diagnostic device.
"""

from . import (  # noqa: F401
    dynamics,
    mutations,
    pathogenicity,
    pockets,
    scoring,
    structures,
    targets,
    uniprot,
    viewer,
)

__version__ = "0.2.0"

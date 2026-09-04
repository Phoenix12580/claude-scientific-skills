"""Dual-engine CNV malignant-cell calling.

Phases stop after QC + labels. No 3CA / NMF / psbulk.
"""

from .run import run_malignant_id

__all__ = ["run_malignant_id"]
__version__ = "1.0.0"

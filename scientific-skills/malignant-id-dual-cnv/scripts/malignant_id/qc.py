"""Phase 04: fail-closed QC. Do not write Cancer cell if this fails."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .io import QC_FAILED, PASSED, write_json


def run_qc(
    labels: pd.DataFrame,
    phase01: dict[str, Any],
    phase02: dict[str, Any],
    contract: dict[str, Any],
    adata_obs: pd.DataFrame | None,
    *,
    outdir: Path,
) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    a1_fp = {pid: rec.get("spike_a1_fp", 0.0) for pid, rec in phase01.items()}
    b_fp = {pid: rec.get("ref_b_fp", 0.0) for pid, rec in phase02.items()}
    max_a1 = max(a1_fp.values()) if a1_fp else 0.0
    max_b = max(b_fp.values()) if b_fp else 0.0
    worst_a1 = max(a1_fp, key=a1_fp.get) if a1_fp else None
    worst_b = max(b_fp, key=b_fp.get) if b_fp else None
    if max_a1 > 0.02:
        errors.append(f"spike A1 FP {max_a1:.4f} > 0.02 (patient {worst_a1})")
    if max_b > 0.005:
        errors.append(f"spike/ref B FP {max_b:.4f} > 0.005 (patient {worst_b})")

    n_normal_malignant = None
    flag = contract.get("normal_sample_flag")
    if flag and flag != "none" and adata_obs is not None:
        col, values = flag["column"], set(map(str, flag["values"]))
        normal_cells = set(adata_obs.index[adata_obs[col].astype(str).isin(values)])
        n_normal_malignant = int(labels.loc[labels["cell"].isin(normal_cells), "is_cancer_cell"].sum())
        if n_normal_malignant != 0:
            errors.append(f"known normal samples have {n_normal_malignant} Cancer cells")

    status = PASSED if not errors else QC_FAILED
    report = {
        "status": status,
        "max_spike_a1_fp": max_a1,
        "max_spike_b_fp": max_b,
        "worst_a1_patient": worst_a1,
        "worst_b_patient": worst_b,
        "n_normal_malignant": n_normal_malignant,
        "n_hc": int((labels["malignant_label"] == "malignant_high_conf").sum()),
        "n_probable": int((labels["malignant_label"] == "probable_malignant").sum()),
        "n_cancer_cell": int(labels["is_cancer_cell"].sum()),
        "errors": errors,
        "warnings": warnings,
    }
    write_json(Path(outdir) / "phase04_qc.json", report)
    return report

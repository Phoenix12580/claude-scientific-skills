from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad

from .io import ContractError, SCHEMA_VERSION, sha256_obj, write_json

REQUIRED = (
    "h5ad_path",
    "organism",
    "genome",
    "sample_key",
    "patient_key",
    "lineage_key",
    "target_categories",
    "tme_categories",
    "spike_source",
    "normal_sample_flag",
    "keep_probable",
)

OBS_CANDIDATES = {
    "sample_key": ("sample_uid", "sample", "orig.ident", "library_id"),
    "patient_key": ("patient_uid", "patient", "donor", "donor_id"),
    "lineage_key": ("cell_type_plot", "celltype_major", "ct.main", "lineage"),
}


def _present(adata: ad.AnnData, col: str) -> bool:
    return col in adata.obs.columns


def suggest_columns(adata: ad.AnnData) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for key, cands in OBS_CANDIDATES.items():
        out[key] = [c for c in cands if _present(adata, c)]
    return out


def validate_contract(contract: dict[str, Any], adata: ad.AnnData | None = None) -> list[str]:
    errors: list[str] = []
    if contract.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"schema_version must be {SCHEMA_VERSION}")
    for k in REQUIRED:
        if k not in contract or contract[k] in (None, "", []):
            errors.append(f"missing required field: {k}")
    if contract.get("organism") not in ("human", "mouse"):
        errors.append("organism must be human or mouse")
    if contract.get("genome") not in ("hg38", "mm10", "mm39"):
        errors.append("genome must be hg38, mm10, or mm39")
    spike = contract.get("spike_source") or {}
    if not isinstance(spike, dict) or spike.get("mode") not in (
        "in_adata_normals",
        "external_atlas",
        "user_barcodes",
    ):
        errors.append("spike_source.mode must be in_adata_normals | external_atlas | user_barcodes")
    if adata is not None:
        for k in ("sample_key", "patient_key", "lineage_key"):
            col = contract.get(k)
            if col and col not in adata.obs.columns:
                errors.append(f"{k}={col!r} not in adata.obs")
        targets = set(contract.get("target_categories") or [])
        if targets:
            present = set(map(str, adata.obs[contract["lineage_key"]].astype(str).unique())) if contract.get("lineage_key") in adata.obs else set()
            missing = targets - present
            if missing:
                errors.append(f"target_categories not in {contract.get('lineage_key')}: {sorted(missing)}")
    return errors


def resolve_contract(
    adata: ad.AnnData,
    overrides: dict[str, Any],
    *,
    outdir: str | Path,
) -> dict[str, Any]:
    """Build a contract from user overrides. Do not invent spike donors."""
    suggestions = suggest_columns(adata)
    contract = dict(overrides)
    contract.setdefault("schema_version", SCHEMA_VERSION)
    for key in ("sample_key", "patient_key", "lineage_key"):
        if not contract.get(key) and len(suggestions.get(key, [])) == 1:
            contract[key] = suggestions[key][0]
    errors = validate_contract(contract, adata)
    if errors:
        payload = {
            "status": "awaiting_user_input",
            "errors": errors,
            "suggested_columns": suggestions,
            "obs_columns": list(map(str, adata.obs.columns)),
        }
        write_json(Path(outdir) / "awaiting_user_input.json", payload)
        raise ContractError("; ".join(errors))
    contract["contract_sha256"] = sha256_obj({k: v for k, v in contract.items() if k != "contract_sha256"})
    write_json(Path(outdir) / "dataset_contract.json", contract)
    return contract

"""Phase 00: split per-patient epi / tme / spike parquet packs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .io import ContractError, write_json


def _is_target(obs: pd.Series, categories: list[str]) -> pd.Series:
    return obs.astype(str).isin(categories)


def _counts_gene_by_cell(adata: ad.AnnData) -> pd.DataFrame:
    x = adata.X
    if hasattr(x, "toarray"):
        mat = x.T.toarray()
    else:
        mat = np.asarray(x).T
    return pd.DataFrame(mat, index=adata.var_names.astype(str), columns=adata.obs_names.astype(str))


def _draw_spike(
    adata: ad.AnnData,
    contract: dict[str, Any],
    query_patient: str,
) -> ad.AnnData:
    spike = contract["spike_source"]
    mode = spike["mode"]
    lineage_key = contract["lineage_key"]
    patient_key = contract["patient_key"]
    targets = contract["target_categories"]
    max_n = int(spike.get("max_cells_per_donor", 500))
    rng = np.random.default_rng(0)

    if mode == "in_adata_normals":
        mask = _is_target(adata.obs[lineage_key], targets)
        mask &= adata.obs[patient_key].astype(str) != str(query_patient)
        normal_pats = spike.get("normal_patient_values") or []
        normal_samps = spike.get("normal_sample_values") or []
        if normal_pats:
            mask &= adata.obs[patient_key].astype(str).isin(normal_pats)
        if normal_samps:
            mask &= adata.obs[contract["sample_key"]].astype(str).isin(normal_samps)
        if not mask.any():
            raise ContractError(f"no legal in-adata spike for patient {query_patient}")
        sub = adata[mask].copy()
        parts = []
        for donor, idx in sub.obs.groupby(patient_key).groups.items():
            if str(donor) == str(query_patient):
                continue
            take = list(idx)
            if len(take) > max_n:
                take = list(rng.choice(take, size=max_n, replace=False))
            parts.append(sub[take])
        if len(parts) < int(spike.get("min_donors", 1)):
            raise ContractError("spike min_donors not met")
        return ad.concat(parts, merge="same") if len(parts) > 1 else parts[0]

    if mode == "user_barcodes":
        barcodes = Path(spike["barcode_file"]).read_text().split()
        keep = [b for b in barcodes if b in adata.obs_names]
        if not keep:
            raise ContractError("user spike barcodes not in adata")
        return adata[keep].copy()

    if mode == "external_atlas":
        atlas = ad.read_h5ad(spike["atlas_h5ad"], backed="r")
        akey = spike.get("atlas_patient_key") or patient_key
        lkey = spike.get("atlas_lineage_key") or lineage_key
        cats = spike.get("atlas_target_categories") or targets
        mask = atlas.obs[lkey].astype(str).isin(cats)
        if akey in atlas.obs:
            mask &= atlas.obs[akey].astype(str) != str(query_patient)
        sub = atlas[mask].to_memory()
        genes = adata.var_names.intersection(sub.var_names)
        if len(genes) < 1000:
            raise ContractError(f"external atlas gene overlap too small: {len(genes)}")
        return sub[:, genes].copy()

    raise ContractError(f"unknown spike mode {mode}")


def pack_patient_runs(
    adata: ad.AnnData,
    contract: dict[str, Any],
    *,
    outdir: str | Path,
) -> dict[str, Any]:
    outdir = Path(outdir)
    pack_root = outdir / "packs"
    pack_root.mkdir(parents=True, exist_ok=True)
    patient_key = contract["patient_key"]
    lineage_key = contract["lineage_key"]
    targets = list(contract["target_categories"])
    tme_cats = list(contract.get("tme_categories") or [])
    min_n = int(contract.get("min_target_cells_per_patient", 200))

    target_mask = _is_target(adata.obs[lineage_key], targets)
    if tme_cats:
        tme_mask = adata.obs[lineage_key].astype(str).isin(tme_cats)
    else:
        tme_mask = ~target_mask

    manifest: dict[str, Any] = {"patients": {}, "skipped": []}
    for pid, idx in adata.obs.groupby(patient_key).groups.items():
        pid = str(pid)
        cells = adata.obs_names[adata.obs_names.isin(idx)]
        sub = adata[cells]
        epi = sub[target_mask.reindex(sub.obs_names).fillna(False).astype(bool)]
        tme = sub[tme_mask.reindex(sub.obs_names).fillna(False).astype(bool)]
        if epi.n_obs < min_n:
            manifest["skipped"].append({"patient": pid, "n_target": int(epi.n_obs), "reason": "below_min_target_cells"})
            continue
        try:
            spike = _draw_spike(adata, contract, pid)
        except ContractError as exc:
            manifest["skipped"].append({"patient": pid, "reason": str(exc)})
            continue
        dst = pack_root / pid
        dst.mkdir(parents=True, exist_ok=True)
        _counts_gene_by_cell(epi).to_parquet(dst / "epi.parquet")
        _counts_gene_by_cell(tme).to_parquet(dst / "tme.parquet")
        _counts_gene_by_cell(spike).to_parquet(dst / "spike.parquet")
        manifest["patients"][pid] = {
            "dir": str(dst),
            "n_epi": int(epi.n_obs),
            "n_tme": int(tme.n_obs),
            "n_spike": int(spike.n_obs),
        }
    if not manifest["patients"]:
        raise ContractError("no patient packs written; check min cells and spike_source")
    write_json(outdir / "pack_manifest.json", manifest)
    return manifest

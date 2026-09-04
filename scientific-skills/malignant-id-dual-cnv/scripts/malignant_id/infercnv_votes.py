"""Phase 02: pyinferCNV B1 burden + B2 ACC."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .io import write_json


def run_one_patient_infercnv(pack_dir: Path, *, genome: str, q_b1: float, q_b2: float) -> dict[str, Any]:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from pyinfercnv import InferCNVConfig, infercnv
    from pyinfercnv.io.genome import load_gene_positions

    epi = pd.read_parquet(pack_dir / "epi.parquet")
    tme = pd.read_parquet(pack_dir / "tme.parquet")
    spike = pd.read_parquet(pack_dir / "spike.parquet")
    mtx = pd.concat([epi, tme, spike], axis=1)
    adata = ad.AnnData(
        X=sp.csr_matrix(mtx.values.T),
        obs=pd.DataFrame(index=mtx.columns),
        var=pd.DataFrame(index=mtx.index),
    )
    epi_n, tme_n, spike_n = epi.shape[1], tme.shape[1], spike.shape[1]
    adata.obs["ref_group"] = "query"
    adata.obs.iloc[epi_n : epi_n + tme_n + spike_n, adata.obs.columns.get_loc("ref_group")] = "ref"
    adata.obs["role"] = ["epi"] * epi_n + ["tme"] * tme_n + ["spike"] * spike_n

    gp = load_gene_positions(genome if genome == "hg38" else genome).set_index("gene_symbol")
    var = adata.var.join(gp, how="left")
    adata = adata[:, var["chromosome"].notna()].copy()
    adata.var = var.loc[adata.var_names]

    cfg = InferCNVConfig(cutoff=0.1, num_threads=1, counts_layer=None)
    res = infercnv(adata, config=cfg, reference_key="ref_group", reference_cat="ref", inplace=False)
    Z = np.asarray(res.cnv_matrix_f64)
    role = adata.obs["role"].reindex(adata.obs_names)
    ref_mask = (role != "epi").to_numpy()
    epi_mask = (role == "epi").to_numpy()

    burden = np.sqrt((Z ** 2).sum(axis=1))
    thr_b1 = float(np.quantile(burden[ref_mask], q_b1))
    pass_b1 = burden > thr_b1

    k = max(3, int(np.ceil(0.01 * epi_mask.sum())))
    top_idx = np.argsort(burden[epi_mask])[-k:]
    epi_pos = np.where(epi_mask)[0][top_idx]
    T = Z[epi_pos].mean(axis=0)
    Zc = Z - Z.mean(axis=1, keepdims=True)
    Tc = T - T.mean()
    acc = (Zc @ Tc) / (np.linalg.norm(Zc, axis=1) * np.linalg.norm(Tc) + 1e-12)
    thr_b2 = float(np.quantile(acc[ref_mask], q_b2))
    pass_b2 = acc > thr_b2

    out = pd.DataFrame(
        {
            "cell": list(adata.obs_names),
            "role": role.to_numpy(),
            "burden": burden,
            "acc": acc,
            "thr_b1": thr_b1,
            "thr_b2": thr_b2,
            "pass_b1": pass_b1,
            "pass_b2": pass_b2,
        }
    )
    out.to_csv(pack_dir / "b_scores.tsv", sep="\t", index=False)
    ref_fp = float(((pass_b1 & pass_b2)[ref_mask]).mean()) if ref_mask.any() else 0.0
    return {"n_cells": int(Z.shape[0]), "ref_b_fp": ref_fp, "epi_b_both_frac": float((pass_b1 & pass_b2)[epi_mask].mean()) if epi_mask.any() else 0.0}


def run_pyinfercnv(pack_manifest: dict[str, Any], contract: dict[str, Any], *, outdir: Path) -> dict[str, Any]:
    q = contract.get("quantiles") or {}
    q_b1 = float(q.get("b1", 0.99))
    q_b2 = float(q.get("b2", 0.99))
    rows = {}
    for pid, rec in pack_manifest["patients"].items():
        rows[pid] = run_one_patient_infercnv(Path(rec["dir"]), genome=contract["genome"], q_b1=q_b1, q_b2=q_b2)
    write_json(Path(outdir) / "phase02_infercnv.json", rows)
    return rows

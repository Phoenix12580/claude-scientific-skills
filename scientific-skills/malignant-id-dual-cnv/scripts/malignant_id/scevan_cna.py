"""Phase 01: SCEVAN CNA → A2 distance + A1 Leiden co-clustering."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from .io import write_json


def _ensure_threads() -> None:
    for v in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMBA_NUM_THREADS",
        "RAYON_NUM_THREADS",
    ):
        os.environ.setdefault(v, "1")


def _import_scevan(scevan_root: Path):
    sys.path.insert(0, str(scevan_root / "tasks" / "scevan" / "src" / "py"))
    from scevan.core import breaks_from_segdf, compute_cna_matrix, nonlinear_smooth, vega_segment
    from scevan.preproc import load_annotations, preprocessing

    return {
        "breaks_from_segdf": breaks_from_segdf,
        "compute_cna_matrix": compute_cna_matrix,
        "nonlinear_smooth": nonlinear_smooth,
        "vega_segment": vega_segment,
        "load_annotations": load_annotations,
        "preprocessing": preprocessing,
    }


def run_one_patient_cna(
    pack_dir: Path,
    *,
    organism: str,
    quantile: float,
    scevan_root: Path,
) -> dict[str, Any]:
    _ensure_threads()
    sv = _import_scevan(scevan_root)
    epi = pd.read_parquet(pack_dir / "epi.parquet")
    spike = pd.read_parquet(pack_dir / "spike.parquet")
    mtx = pd.concat([epi, spike], axis=1)
    annotations = sv["load_annotations"]("human" if organism == "human" else "mouse")
    res = sv["preprocessing"](mtx, annotations, find_confident=False)
    cn, annot = res["count_norm"], res["annot"]
    ref_cols = [c for c in spike.columns if c in cn.columns]
    X = cn.values
    basel = np.median(X[:, [cn.columns.get_loc(c) for c in ref_cols]], axis=1)
    relat = sv["nonlinear_smooth"](X - basel[:, None])
    chr_arr = annot["seqnames"].values.astype(int)
    segdf = sv["vega_segment"](relat, chr_arr, annot["end"].values.astype(int), beta=0.5)
    breaks = sv["breaks_from_segdf"](segdf, relat.shape[0])
    cna = sv["compute_cna_matrix"](relat, breaks, np.ones(len(breaks) - 1, bool))
    cna = cna - cna.mean(axis=0, keepdims=True)

    col_list = list(cn.columns)
    ref_idx = np.array([col_list.index(c) for c in ref_cols])
    center = np.median(cna[:, ref_idx], axis=1)
    dist = np.linalg.norm(cna - center[:, None], axis=0)
    thr = float(np.quantile(dist[ref_idx], quantile))
    in_spike = np.asarray(pd.Index(col_list).isin(ref_cols))
    a2 = dist > thr

    # A1: Leiden on cells × segments, spike-depleted clusters = aneuploid
    cna_cells = cna.T.astype(np.float32)
    adata = ad.AnnData(cna_cells)
    adata.obs_names = col_list
    adata.obs["is_spike"] = in_spike
    n_pcs = min(30, max(2, cna_cells.shape[1] - 1), max(2, cna_cells.shape[0] - 1))
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1), use_rep="X_pca")
    sc.tl.leiden(adata, flavor="igraph", directed=False, n_iterations=2, key_added="cnv_leiden")
    spike_frac = adata.obs.groupby("cnv_leiden")["is_spike"].mean()
    global_spike = float(adata.obs["is_spike"].mean())
    aneuploid = set(spike_frac[spike_frac < global_spike].index.astype(str))
    a1 = adata.obs["cnv_leiden"].astype(str).isin(aneuploid).to_numpy()

    out = pd.DataFrame(
        {
            "cell": col_list,
            "is_spike": in_spike,
            "cnv_dist": dist,
            "is_tumor_a2": a2,
            "thr_a2": thr,
            "cnv_leiden": adata.obs["cnv_leiden"].astype(str).to_numpy(),
            "is_tumor_a1": a1,
        }
    )
    keep = np.asarray(pd.Index(col_list).isin(epi.columns)) | in_spike
    out = out.loc[keep]
    dst = pack_dir
    out.to_csv(dst / "a_scores.tsv", sep="\t", index=False)
    pd.DataFrame(cna[:, keep].T.astype(np.float32), index=out["cell"].tolist()).to_parquet(dst / "cna_matrix.parquet")
    spike_fp = float(out.loc[out.is_spike, "is_tumor_a1"].mean()) if out.is_spike.any() else 0.0
    return {
        "n_cells": int(out.shape[0]),
        "spike_a1_fp": spike_fp,
        "epi_a1_frac": float(out.loc[~out.is_spike, "is_tumor_a1"].mean()) if (~out.is_spike).any() else 0.0,
        "n_aneuploid_clusters": len(aneuploid),
    }


def run_scevan_cna(pack_manifest: dict[str, Any], contract: dict[str, Any], env: dict[str, Any], *, outdir: Path) -> dict[str, Any]:
    root = Path(env["scevan"]["root"])
    q = float((contract.get("quantiles") or {}).get("a2", 0.99))
    rows = {}
    for pid, rec in pack_manifest["patients"].items():
        rows[pid] = run_one_patient_cna(Path(rec["dir"]), organism=contract["organism"], quantile=q, scevan_root=root)
    write_json(Path(outdir) / "phase01_scevan.json", rows)
    return rows

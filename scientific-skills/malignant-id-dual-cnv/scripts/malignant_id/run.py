"""Orchestrate phases 00–04. Stop after QC. Never overwrite the source h5ad."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad

from .contract import ContractError, resolve_contract
from .environment import capture_environment, resolve_scevan_root
from .infercnv_votes import run_pyinfercnv
from .io import AWAITING, FAILED, PASSED, QC_FAILED, SCHEMA_VERSION, sha256_file, utc_now, write_json
from .labels import synthesize_labels
from .pack import pack_patient_runs
from .qc import run_qc
from .scevan_cna import run_scevan_cna


def _write_h5ad(adata: ad.AnnData, labels, contract: dict[str, Any], outdir: Path) -> Path:
    out = adata.copy()
    lab = labels.set_index("cell")
    out.obs["malignant_label"] = "nonmalignant_tme"
    shared = out.obs_names.intersection(lab.index)
    out.obs.loc[shared, "malignant_label"] = lab.loc[shared, "malignant_label"].astype(str)
    for col in ("is_tumor_a1", "is_tumor_a2", "pass_b1", "pass_b2", "cnv_dist", "burden", "acc"):
        if col in lab.columns:
            dest = {
                "is_tumor_a1": "cnv_a1",
                "is_tumor_a2": "cnv_a2",
                "pass_b1": "cnv_b1",
                "pass_b2": "cnv_b2",
                "cnv_dist": "cnv_dist",
                "burden": "cnv_burden",
                "acc": "cnv_acc",
            }[col]
            out.obs[dest] = False if col.startswith(("is_", "pass_")) else float("nan")
            out.obs.loc[shared, dest] = lab.loc[shared, col]
    plot_key = contract.get("lineage_key")
    if plot_key and plot_key in out.obs:
        cancer = shared[lab.loc[shared, "is_cancer_cell"].astype(bool)]
        out.obs[plot_key] = out.obs[plot_key].astype(str)
        out.obs.loc[cancer, plot_key] = "Cancer cell"
    dest = Path(outdir) / "annotated_malignant.h5ad"
    out.write_h5ad(dest)
    return dest


def run_malignant_id(
    *,
    h5ad_path: str,
    outdir: str,
    contract_overrides: dict[str, Any] | None = None,
    scevan_root: str | None = None,
) -> dict[str, Any]:
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)
    started = utc_now()
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "running",
        "started_at": started,
        "finished_at": None,
        "phases": {},
        "errors": [],
        "warnings": [],
    }
    try:
        adata = ad.read_h5ad(h5ad_path)
        overrides = dict(contract_overrides or {})
        overrides["h5ad_path"] = str(Path(h5ad_path).resolve())
        contract = resolve_contract(adata, overrides, outdir=outdir_p)
        root = resolve_scevan_root(scevan_root)
        if root is None:
            raise ContractError("SCEVAN_ROOT not found; set it or pass scevan_root")
        env = capture_environment(outdir=outdir_p, scevan_root=root, h5ad_path=h5ad_path)
        manifest["environment_path"] = str(outdir_p / "environment.json")
        manifest["contract_sha256"] = contract["contract_sha256"]
        manifest["input_h5ad_sha256"] = sha256_file(h5ad_path)

        packs = pack_patient_runs(adata, contract, outdir=outdir_p)
        manifest["phases"]["phase00_pack"] = "ok"
        manifest["n_patients"] = len(packs["patients"])
        manifest["n_target"] = int(sum(p["n_epi"] for p in packs["patients"].values()))
        manifest["n_tme"] = int(sum(p["n_tme"] for p in packs["patients"].values()))
        manifest["n_spike"] = int(sum(p["n_spike"] for p in packs["patients"].values()))

        p01 = run_scevan_cna(packs, contract, env, outdir=outdir_p)
        manifest["phases"]["phase01_scevan"] = "ok"
        p02 = run_pyinfercnv(packs, contract, outdir=outdir_p)
        manifest["phases"]["phase02_infercnv"] = "ok"
        labels = synthesize_labels(packs, contract, outdir=outdir_p)
        manifest["phases"]["phase03_labels"] = "ok"
        qc = run_qc(labels, p01, p02, contract, adata.obs, outdir=outdir_p)
        manifest["phases"]["phase04_qc"] = qc["status"]
        manifest["n_hc"] = qc["n_hc"]
        manifest["n_probable"] = qc["n_probable"]
        manifest["n_cancer_cell"] = qc["n_cancer_cell"]
        if qc["status"] != PASSED:
            manifest["status"] = QC_FAILED
            manifest["errors"] = qc["errors"]
        else:
            dest = _write_h5ad(adata, labels, contract, outdir_p)
            manifest["status"] = PASSED
            manifest["annotated_h5ad"] = str(dest)
        manifest["finished_at"] = utc_now()
        write_json(outdir_p / "run_manifest.json", manifest)
        return manifest
    except ContractError as exc:
        manifest["status"] = AWAITING
        manifest["errors"] = [str(exc)]
        manifest["finished_at"] = utc_now()
        write_json(outdir_p / "run_manifest.json", manifest)
        return manifest
    except Exception as exc:
        manifest["status"] = FAILED
        manifest["errors"] = [repr(exc)]
        manifest["finished_at"] = utc_now()
        write_json(outdir_p / "run_manifest.json", manifest)
        raise

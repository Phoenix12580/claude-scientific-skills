from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from .io import sha256_file, utc_now, write_json


def _pkg_version(name: str) -> str | None:
    try:
        mod = __import__(name)
        return str(getattr(mod, "__version__", None) or "unknown")
    except Exception:
        return None


def _git_commit(root: Path) -> str | None:
    git = shutil.which("git")
    if git is None or not (root / ".git").exists():
        return None
    try:
        out = subprocess.check_output(
            [git, "-C", str(root), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip()
    except Exception:
        return None


def resolve_scevan_root(explicit: str | None = None) -> Path | None:
    if explicit:
        p = Path(explicit).expanduser().resolve()
        return p if p.exists() else None
    env = os.environ.get("SCEVAN_ROOT")
    if env:
        p = Path(env).expanduser().resolve()
        return p if p.exists() else None
    fallback = Path.home() / "scevan-001"
    return fallback if fallback.exists() else None


def capture_environment(
    *,
    outdir: str | Path,
    scevan_root: str | Path | None,
    h5ad_path: str | Path,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    root = Path(scevan_root) if scevan_root else None
    pixi_lock = (root / "pixi.lock") if root else None
    record: dict[str, Any] = {
        "schema_version": "1.0.0",
        "captured_at": utc_now(),
        "python": {
            "version": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
        },
        "packages": {
            "numpy": _pkg_version("numpy"),
            "pandas": _pkg_version("pandas"),
            "anndata": _pkg_version("anndata"),
            "scanpy": _pkg_version("scanpy"),
            "scipy": _pkg_version("scipy"),
            "pyinfercnv": _pkg_version("pyinfercnv"),
            "scevan": _pkg_version("scevan"),
        },
        "threads": {
            k: os.environ.get(k)
            for k in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMBA_NUM_THREADS",
                "RAYON_NUM_THREADS",
                "PYTHONNOUSERSITE",
            )
        },
        "scevan": {
            "root": str(root) if root else None,
            "git_commit": _git_commit(root) if root else None,
            "pixi_lock_sha256": sha256_file(pixi_lock) if pixi_lock and pixi_lock.exists() else None,
            "run_v2_cna": str(root / "tasks/scevan/src/run_v2_cna.py") if root else None,
            "run_v2_infercnv": str(root / "tasks/scevan/src/run_v2_infercnv.py") if root else None,
        },
        "input_h5ad": str(Path(h5ad_path).resolve()),
        "hostname": platform.node(),
        "cpu_count": os.cpu_count(),
    }
    if extra:
        record["extra"] = extra
    write_json(Path(outdir) / "environment.json", record)
    return record

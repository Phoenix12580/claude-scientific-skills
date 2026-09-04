"""Snapshot the interpreter that is actually running this process.

The environment is whatever python launched the CLI / orchestrator — not a
guessed conda name and not a later reconstruction. Capture it once at the
start of the run and write environment.json + requirements.lock.txt.
"""

from __future__ import annotations

import os
import platform
import shutil
import site
import subprocess
import sys
from pathlib import Path
from typing import Any

from .io import sha256_file, utc_now, write_json

PINNED_IMPORTS = (
    "numpy",
    "pandas",
    "anndata",
    "scanpy",
    "scipy",
    "sklearn",
    "pyinfercnv",
    "scevan",
    "numba",
    "igraph",
    "leidenalg",
)


def _module_record(name: str) -> dict[str, Any]:
    rec: dict[str, Any] = {"name": name, "imported": False, "version": None, "file": None}
    try:
        mod = __import__(name)
    except Exception as exc:
        rec["import_error"] = type(exc).__name__
        return rec
    rec["imported"] = True
    rec["version"] = getattr(mod, "__version__", None)
    rec["file"] = getattr(mod, "__file__", None)
    return rec


def _pip_freeze(python: str) -> list[str]:
    """Freeze THIS executable's site-packages. Never call a different python."""
    try:
        out = subprocess.check_output(
            [python, "-m", "pip", "freeze"],
            stderr=subprocess.STDOUT,
            text=True,
            timeout=120,
        )
        return [ln.strip() for ln in out.splitlines() if ln.strip() and not ln.startswith("#")]
    except Exception as exc:
        return [f"# pip freeze failed: {type(exc).__name__}: {exc}"]


def _git_commit(root: Path) -> str | None:
    git = shutil.which("git")
    if git is None or not (root / ".git").exists():
        return None
    try:
        return subprocess.check_output(
            [git, "-C", str(root), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
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
    """Record the live runtime. Call this from the same process that will run CNV."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    root = Path(scevan_root).resolve() if scevan_root else None
    pixi_lock = (root / "pixi.lock") if root else None
    exe = Path(sys.executable).resolve()
    freeze = _pip_freeze(str(exe))
    lock_path = outdir / "requirements.lock.txt"
    lock_path.write_text("\n".join(freeze) + "\n")

    record: dict[str, Any] = {
        "schema_version": "1.1.0",
        "captured_at": utc_now(),
        "note": "This is a snapshot of the running interpreter (sys.executable), not a reconstructed env.",
        "python": {
            "version": sys.version,
            "version_info": list(sys.version_info[:3]),
            "executable": str(exe),
            "executable_sha256": sha256_file(exe, max_bytes=1024 * 1024) if exe.exists() else None,
            "prefix": sys.prefix,
            "base_prefix": getattr(sys, "base_prefix", None),
            "executable_real": str(Path(os.path.realpath(exe))),
            "platform": platform.platform(),
            "implementation": platform.python_implementation(),
        },
        "site": {
            "usersite": site.getusersitepackages() if hasattr(site, "getusersitepackages") else None,
            "sitepackages": site.getsitepackages() if hasattr(site, "getsitepackages") else None,
            "PYTHONNOUSERSITE": os.environ.get("PYTHONNOUSERSITE"),
            "CONDA_PREFIX": os.environ.get("CONDA_PREFIX"),
            "CONDA_DEFAULT_ENV": os.environ.get("CONDA_DEFAULT_ENV"),
            "VIRTUAL_ENV": os.environ.get("VIRTUAL_ENV"),
            "PIXI_ENVIRONMENT_NAME": os.environ.get("PIXI_ENVIRONMENT_NAME"),
        },
        "modules": {name: _module_record(name) for name in PINNED_IMPORTS},
        "pip_freeze_file": str(lock_path),
        "pip_freeze_sha256": sha256_file(lock_path),
        "pip_freeze_n_lines": len([ln for ln in freeze if not ln.startswith("#")]),
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
        "cwd": os.getcwd(),
        "argv": list(sys.argv),
    }
    if extra:
        record["extra"] = extra
    write_json(outdir / "environment.json", record)
    return record

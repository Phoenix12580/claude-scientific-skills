from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
AWAITING = "awaiting_user_input"
FAILED = "failed"
PASSED = "passed_qc"
QC_FAILED = "qc_failed"


class ContractError(RuntimeError):
    """Missing or illegal dataset contract field."""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: str | Path, max_bytes: int | None = None) -> str:
    p = Path(path)
    h = hashlib.sha256()
    with p.open("rb") as fh:
        if max_bytes is None:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        else:
            h.update(fh.read(max_bytes))
    return h.hexdigest()


def sha256_obj(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def write_json(path: str | Path, obj: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")
    tmp.replace(path)
    return path


def read_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text())

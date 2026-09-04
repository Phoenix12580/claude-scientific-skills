#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SKILL_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from malignant_id.run import run_malignant_id  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Dual-engine CNV malignant identification (stops after QC).")
    p.add_argument("--h5ad", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--contract-json", default=None, help="Optional pre-filled dataset_contract.json")
    p.add_argument("--scevan-root", default=None)
    p.add_argument("--organism", choices=["human", "mouse"])
    p.add_argument("--genome", choices=["hg38", "mm10", "mm39"])
    p.add_argument("--sample-key")
    p.add_argument("--patient-key")
    p.add_argument("--lineage-key")
    p.add_argument("--target-categories", nargs="+")
    p.add_argument("--keep-probable", action="store_true", default=None)
    p.add_argument("--no-probable", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    overrides = {}
    if args.contract_json:
        overrides.update(json.loads(Path(args.contract_json).read_text()))
    for key in ("organism", "genome", "sample_key", "patient_key", "lineage_key"):
        val = getattr(args, key.replace("-", "_") if False else key)
        # argparse dest uses underscore
    mapping = {
        "organism": args.organism,
        "genome": args.genome,
        "sample_key": args.sample_key,
        "patient_key": args.patient_key,
        "lineage_key": args.lineage_key,
        "target_categories": args.target_categories,
    }
    for k, v in mapping.items():
        if v:
            overrides[k] = v
    if args.no_probable:
        overrides["keep_probable"] = False
    elif args.keep_probable:
        overrides["keep_probable"] = True
    manifest = run_malignant_id(
        h5ad_path=args.h5ad,
        outdir=args.outdir,
        contract_overrides=overrides,
        scevan_root=args.scevan_root,
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if manifest.get("status") in {"passed_qc", "awaiting_user_input"} else 2


if __name__ == "__main__":
    raise SystemExit(main())

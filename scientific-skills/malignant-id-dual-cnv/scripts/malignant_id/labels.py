"""Phase 03: synthesize hc / probable / Cancer cell labels."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .io import write_json


def synthesize_labels(
    pack_manifest: dict[str, Any],
    contract: dict[str, Any],
    *,
    outdir: Path,
) -> pd.DataFrame:
    keep_probable = bool(contract.get("keep_probable", True))
    frames = []
    for pid, rec in pack_manifest["patients"].items():
        d = Path(rec["dir"])
        a = pd.read_csv(d / "a_scores.tsv", sep="\t")
        b = pd.read_csv(d / "b_scores.tsv", sep="\t")
        a = a[~a["is_spike"]].copy()
        b_epi = b[b["role"] == "epi"].copy()
        m = a.merge(b_epi, on="cell", how="inner")
        m["patient"] = pid
        a1 = m["is_tumor_a1"].astype(bool)
        b1 = m["pass_b1"].astype(bool)
        b2 = m["pass_b2"].astype(bool)
        hc = a1 & b1 & b2
        probable = a1 & (b1 ^ b2)
        if not keep_probable:
            probable = pd.Series(False, index=m.index)
        label = pd.Series("nonmalignant_lineage", index=m.index)
        label[a1 & ~b1 & ~b2] = "uncertain"
        label[probable] = "probable_malignant"
        label[hc] = "malignant_high_conf"
        m["malignant_label"] = label
        m["is_cancer_cell"] = label.isin(["malignant_high_conf", "probable_malignant"])
        frames.append(m)
    labels = pd.concat(frames, ignore_index=True)
    labels.to_csv(Path(outdir) / "cell_labels.tsv", sep="\t", index=False)
    write_json(
        Path(outdir) / "phase03_counts.json",
        labels["malignant_label"].value_counts().to_dict() | {"n_cancer_cell": int(labels["is_cancer_cell"].sum())},
    )
    return labels

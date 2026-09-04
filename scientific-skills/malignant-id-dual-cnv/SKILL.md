---
name: malignant-id-dual-cnv
description: Dual-engine CNV malignant-cell calling for scRNA-seq on a NEW dataset. Epithelium is necessary but not sufficient. SCEVAN A1 Leiden co-clustering plus pyinferCNV B1 burden and B2 ACC votes define high-confidence and probable malignant cells. Use when the user asks to identify tumor/malignant cells from annotated scRNA, run SCEVAN plus inferCNV/pyinferCNV gates, apply a CNV triple gate, or write Cancer cell labels that are CNV-defined rather than marker-defined. Dataset-specific names (Henry D17, PRAD paths, KLK3) are examples only — resolve all inputs from the current AnnData. Do not use for bulk CNV, scATAC, CopyKAT-only scoring, marker-only tumor calling, or 3CA/NMF/psbulk.
---

# Dual-engine CNV malignant identification

## Goal

Turn the **current** annotated scRNA AnnData into cell-level malignant labels:

- `malignant_high_conf` (hc)
- `probable_malignant`
- `Cancer cell = hc ∪ probable`

Do **not** call epithelial/lineage cells malignant because they express lineage markers.

This skill is dataset-agnostic. Methods and gates are fixed. Paths, column names, spike donors, and species are **inputs to resolve**, not constants.

Run the packaged functions; do not reimplement the engines in a notebook.

```text
scripts/run_malignant_id.py
scripts/malignant_id/
  contract.py        # resolve_contract / validate_contract
  environment.py     # capture_environment → environment.json
  pack.py            # Phase 00 pack_patient_runs
  scevan_cna.py      # Phase 01 A1/A2
  infercnv_votes.py  # Phase 02 B1/B2
  labels.py          # Phase 03 hc/probable
  qc.py              # Phase 04 fail-closed QC
  run.py             # orchestrator, stops after QC
```

```bash
python scripts/run_malignant_id.py --h5ad PATH.h5ad --outdir RUNDIR \
  --organism human --genome hg38 --patient-key patient_uid \
  --lineage-key cell_type_plot --target-categories EPI \
  --scevan-root "$SCEVAN_ROOT"
```

Missing contract fields write `awaiting_user_input.json` and stop. QC failure does not write `Cancer cell` into a new h5ad.

## When to use

Use when the user wants:

- malignant / tumor cell identification from scRNA
- SCEVAN + inferCNV / pyinferCNV dual votes
- high-conf vs probable malignant tiers
- a CNV-defined `Cancer cell` label on the epithelial/target lineage

Do **not** use for:

- bulk CNV or WGS/WES
- scATAC
- CopyKAT-only workflows (`sc_cnv_infercnvpy_copykat`)
- marker-only tumor calling
- 3CA / NMF / metaprogram / psbulk (stop at malignant labels)

## Dataset contract — fill this before any compute

Write `dataset_contract.json` from the current h5ad. If any required field is missing, write `awaiting_user_input` and stop. **Do not reuse PRAD contract values.**

Required:

| field | meaning |
|---|---|
| `h5ad_path` | current counts AnnData (do not overwrite) |
| `organism` | `human` or `mouse` |
| `genome` | `hg38` or `mm10`/`mm39` |
| `sample_key` | obs column for sample |
| `patient_key` | obs column for biological donor/patient |
| `lineage_key` | obs column with target vs TME |
| `target_categories` | categories to test (usually epithelium; other solid-tumor lineages allowed if user names them) |
| `tme_categories` | all remaining non-target categories |
| `spike_source` | how external normal target-lineage cells are obtained (see Spike) |
| `normal_sample_flag` | obs column/values for known normal samples, or `none` |
| `keep_probable` | `true` unless the user forbids the probable tier |

Optional: `condition_key`, `dataset_key`, `min_target_cells_per_patient` (default 200).

Example only (PRAD v2, do not copy blindly):

- human / hg38 / `sample_uid` / `patient_uid` / `cell_type_plot` / target=`EPI`
- spike = Henry D17/D27/D35 epithelium, 500 cells each, self-exclusion
- proven yield on that dataset: 256,024 → 144,592 epithelium → 63,731 Cancer cell

## Non-negotiable gates

Stop and ask, do not invent:

1. Target-lineage vs TME categories exist. TME never receives a malignant call.
2. Analysis unit is `patient_key` (merge samples of the same patient first).
3. Spike cells are **external** target-lineage normals. Self-as-reference is forbidden.
4. `X` is integer raw counts. Do not feed lognorm to either engine.
5. Final labels require target-lineage ∧ A1. B1/B2 only decide hc vs probable.
6. A2 is QC / heatmap only. It does **not** enter the final gate.
7. Thresholds come from **spike/ref quantiles**, never from query self-quantiles.
8. If spike FP or known-normal leakage fails QC, do not write `Cancer cell`.

## Spike (must be re-chosen per dataset)

Spike = diploid **target-lineage** cells that are not from the query patient.

Acceptable, in order:

1. Normal samples already in this AnnData (same tissue lineage, different patients)
2. An external normal atlas of the same species/tissue lineage, intersected to the query gene space
3. User-supplied spike barcodes

Forbidden:

- using the query patient's own cells as spike
- using TME as the only SCEVAN baseline (TME is inferCNV `ref`, not SCEVAN spike)
- mixing species
- copying PRAD Henry D17/D27/D35 when those donors are not in the current data

Default draw: up to 500 cells per normal donor, ≥2 donors if available. If no legal spike exists, stop.

## Execution environment

Prefer a dedicated SCEVAN env if present:

```bash
cd "$SCEVAN_ROOT"   # resolve on this machine; do not assume ~/scevan-001
export PYTHONNOUSERSITE=1
```

If that env is absent, install/run the same two engines elsewhere, but keep the vote definitions identical.

Hard limits when using the pixi box:

- ≤ 12 processes
- `mp.get_context("spawn")`
- each process `OMP/OPENBLAS/MKL/NUMBA/RAYON_NUM_THREADS=1`
- long jobs: `nohup ... > log 2>&1 &`

The Python package wraps those engines. Call `malignant_id.run.run_malignant_id(...)` or the CLI above. Do not hard-code `/home/y413007/...` in a new-dataset run. Each run writes `environment.json` (python, packages, SCEVAN git, pixi.lock hash, thread caps) and `run_manifest.json`.

## Phase 00 — inspect and pack

1. Backed-read the current h5ad. Confirm counts, unique cell IDs, contract fields.
2. Split one directory per patient: `epi.parquet` (target lineage), `tme.parquet`, `spike.parquet`.
3. Write `manifest.json` and `dataset_contract.json`.
4. Independent run if target-lineage cells ≥ `min_target_cells_per_patient`; else merge same-condition neighbors or emit `uncertain` only.

Do not start Stage1 if spike includes the query patient.

## Phase 01 — SCEVAN CNA (A1 / A2)

Input: `epi + spike` only.

1. `preprocessing(mtx, annot_<organism>, find_confident=False)`
2. Baseline = median of spike cells
3. `relat = nonlinear_smooth(X - baseline)`
4. `vega_segment(..., beta=0.5)`
5. `compute_cna_matrix` then center **per cell**

A2 (auxiliary):

- center = median CNA vector of spike
- `cnv_dist` = Euclidean distance to center
- `thr` = 99th percentile of spike distances
- `A2 = dist > thr`

A1 (enters the final gate):

1. PCA + neighbors + Leiden (`flavor="igraph", directed=False`) on cells × segments
2. Score clusters by spike fraction
3. Aneuploid cluster = spike-depleted; diploid/reference = spike-enriched
4. `A1 = target-lineage cell in an aneuploid cluster`

Outputs: `cna_matrix.parquet`, `a2_scores.tsv`, A1 labels.

## Phase 02 — pyinferCNV (B1 / B2)

Input: `epi + tme + spike`.

1. AnnData `X = csr counts`. `counts_layer=None`.
2. `ref` = TME+spike, `query` = target lineage
3. Keep genes with coordinates for `genome`
4. `InferCNVConfig(cutoff=0.1, num_threads=1)` → `Z`

B1: `burden = sqrt(sum(Z^2))`; `thr_b1` = 99th percentile of **ref**; `B1 = burden > thr_b1`.

B2: template `T` = mean Z of query cells with top 1% burden (minimum 3); `ACC` = Pearson(cell Z, T); `thr_b2` = 99th percentile of **ref**; `B2 = ACC > thr_b2`.

Output: `b_scores.tsv`. Default pool size 3.

## Phase 03 — label synthesis

Only target-lineage cells get malignant labels. TME is always `nonmalignant_tme`.

| label | rule |
|---|---|
| `malignant_high_conf` | target ∧ A1 ∧ B1 ∧ B2 |
| `probable_malignant` | target ∧ A1 ∧ (B1 XOR B2) |
| `uncertain` | target ∧ A1 ∧ ¬B1 ∧ ¬B2 |
| `nonmalignant_lineage` | target ∧ ¬A1 |
| `nonmalignant_tme` | not target |

`Cancer cell = hc ∪ probable` (keep this name unless the user wants another).

Keep probable by default for low-CNA tumors. Set `keep_probable=false` only if the user says so.

Write a **new** h5ad:

- `obs['malignant_label']`
- `obs['cnv_a1']`, `obs['cnv_a2']`, `obs['cnv_b1']`, `obs['cnv_b2']`
- `obs['cnv_dist']`, `obs['cnv_burden']`, `obs['cnv_acc']`
- overwrite a plot column to `Cancer cell` **only** for hc+probable target cells

Never overwrite the source file.

## Phase 04 — QC

Fail closed:

| check | pass bar |
|---|---|
| Spike A1 false positive | ≤ 2% |
| Spike/ref B false positive | ≤ 0.5% |
| Known normal samples | 0 malignant cells, if `normal_sample_flag` exists |
| Self-as-spike | forbidden |

Do not relax query self-quantiles to rescue yield.

## Endpoint

Stop after labels + QC + new h5ad + summary tables. Do not start 3CA/NMF/psbulk.

## What the agent must report

1. Contract fields actually used (columns, spike source, genome)
2. n target, n TME, n spike, n patients/runs
3. Spike A1 FP and B FP (max and which patient)
4. Known-normal malignant count (0 if normals exist)
5. hc / probable / Cancer cell counts
6. Statement: Cancer cell is CNV-defined, not marker-defined; Cancer cell ⊂ target lineage

## What the agent must never do

- Copy PRAD donor IDs, paths, or counts onto a new dataset
- Call the target lineage malignant from markers alone
- Use A2 in the final gate
- Use query-derived thresholds
- Use a patient as its own spike
- Feed lognorm to either engine
- Continue into 3CA/NMF/psbulk

## Minimal claim language

Correct: "Cancer cells are target-lineage cells that passed SCEVAN A1 and at least one pyinferCNV vote (B1 XOR/AND B2)."

Incorrect: "Epithelial cells were annotated as tumor."

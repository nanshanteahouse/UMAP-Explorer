#!/usr/bin/env python3
"""
build_integrated.py - Integrate all datasets into shared UMAP space via Harmony.

Loads all 04_clustered.h5ad files, finds shared highly variable genes,
concatenates into a single AnnData, runs PCA + Harmony integration + UMAP,
and exports parquet files for the integrated cross-dataset view.

Usage:
    python build_integrated.py \
        --data-dir /path/to/neurobiology \
        --output-dir /path/to/vis_website/data/integrated
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Cell type unification mapping
# ---------------------------------------------------------------------------

CELL_TYPE_UNIFICATION: dict[str, str] = {
    # Retinal cell types — standardize across datasets
    "RPC": "RPC",
    "RPCs": "RPC",
    "Retinal Progenitor": "RPC",
    "RPC_1": "RPC",
    "RPC_2": "RPC",
    "RPC_3": "RPC",
    "RPC_4": "RPC",
    "RPC_5": "RPC",
    "RPCs_1": "RPC",
    "RPCs_2": "RPC",
    "RPCs_3": "RPC",
    "Neurogenic": "Neurogenic",
    "Neurogenic Cells": "Neurogenic",
    "Neurogenic_Cells": "Neurogenic",
    "RGC": "RGC",
    "RGCs": "RGC",
    "Retinal Ganglion Cells": "RGC",
    "Cone": "Cone",
    "Cones": "Cone",
    "Cone Photoreceptor": "Cone",
    "Cone_Photoreceptor": "Cone",
    "Rod": "Rod",
    "Rods": "Rod",
    "Rod Photoreceptor": "Rod",
    "Rod_Photoreceptor": "Rod",
    "Bipolar": "Bipolar",
    "Bipolar Cell": "Bipolar",
    "Bipolar_Cell": "Bipolar",
    "Amacrine": "Amacrine",
    "Amacrine Cell": "Amacrine",
    "Amacrine_Cell": "Amacrine",
    "Muller": "Muller Glia",
    "Muller Glia": "Muller Glia",
    "Muller_Glia": "Muller Glia",
    "MullerGlia": "Muller Glia",
    "Horizontal": "Horizontal",
    "Horizontal Cell": "Horizontal",
    "Horizontal_Cell": "Horizontal",
    "Microglia": "Microglia",
    "Astrocyte": "Astrocyte",
    "Astrocytes": "Astrocyte",
    "Vascular": "Vascular",
    "Endothelial": "Vascular",
    "RPE": "RPE",
    "Retinal Pigment Epithelium": "RPE",
    # Brain cell types (from GSE140877 / GSE268630)
    "Radial Glia": "Radial Glia",
    "Radial_Glia": "Radial Glia",
    "Excitatory Neuron": "Excitatory Neuron",
    "Excitatory_Neuron": "Excitatory Neuron",
    "Inhibitory Neuron": "Inhibitory Neuron",
    "Inhibitory_Neuron": "Inhibitory Neuron",
    "Oligodendrocyte": "Oligodendrocyte",
    "Oligodendrocytes": "Oligodendrocyte",
    "OPC": "OPC",
    "Pericyte": "Pericyte",
    "IP/Neuroprogenitor": "Neuroprogenitor",
    "Neuroprogenitor": "Neuroprogenitor",
    "Neural Progenitors": "Neuroprogenitor",
    "Neural_Progenitors": "Neuroprogenitor",
    "Ventricular Zone": "Neuroprogenitor",
    "IP/Neuroprogenitors": "Neuroprogenitor",
    "Interneurons": "Interneuron",
    "Motor Neurons": "Motor Neuron",
    "Motor_Neurons": "Motor Neuron",
    "Neurons": "Neuron",
}


def _unify_cell_type(label: str) -> str:
    """Map a raw cell-type label to its unified name."""
    if not isinstance(label, str):
        return "Unknown"
    label_clean = label.strip()
    return CELL_TYPE_UNIFICATION.get(label_clean, label_clean)


# ---------------------------------------------------------------------------
# Path resolution (same strategy as build_data_cache.py)
# ---------------------------------------------------------------------------


def find_dataset_paths(data_dir: str) -> dict[str, str]:
    """
    Scan for all 04_clustered.h5ad files under *data_dir*.

    Returns dict: {gse_id: absolute_path}
    """
    dataset_paths: dict[str, str] = {}

    patterns = [
        lambda gse: os.path.join(gse, "results", "h5ad", "04_clustered.h5ad"),
        lambda gse: os.path.join(gse, f"{gse}_results", "h5ad", "04_clustered.h5ad"),
        lambda gse: os.path.join(gse, "pipeline", "results", "h5ad", "04_clustered.h5ad"),
    ]

    if not os.path.isdir(data_dir):
        print(f"[ERROR] data_dir does not exist: {data_dir}")
        return dataset_paths

    for entry in sorted(os.listdir(data_dir)):
        entry_path = os.path.join(data_dir, entry)
        if not os.path.isdir(entry_path) or not entry.startswith("GSE"):
            continue

        found = False
        for pattern_fn in patterns:
            candidate = os.path.join(data_dir, pattern_fn(entry))
            if os.path.isfile(candidate):
                dataset_paths[entry] = os.path.abspath(candidate)
                found = True
                break

        if not found:
            print(f"[INFO]  Skipping {entry}: no 04_clustered.h5ad found")

    return dataset_paths


# ---------------------------------------------------------------------------
# Shared HVG discovery
# ---------------------------------------------------------------------------


def find_shared_hvgs(adatas: list[sc.AnnData], n_top: int = 3000) -> list[str]:
    """Find the intersection of highly variable genes across all datasets.

    Returns up to *n_top* genes, ranked by average HVG rank across datasets.
    """
    hvg_sets: list[set[str]] = []
    for ad in adatas:
        if "highly_variable" not in ad.var.columns:
            print("  [WARN] Dataset missing 'highly_variable' column — skipping HVG filter")
            return list(ad.var_names)
        hvg_genes = set(ad.var_names[ad.var["highly_variable"].values])
        hvg_sets.append(hvg_genes)

    if not hvg_sets:
        return []

    # Intersection of HVG sets
    shared: set[str] = hvg_sets[0]
    for s in hvg_sets[1:]:
        shared &= s

    shared_list = sorted(shared)
    print(f"  Shared HVGs across all datasets: {len(shared_list)}")

    if not shared_list:
        print("  [WARN] No shared HVGs found — falling back to intersection of all genes")
        # Use all gene intersection as fallback
        all_gene_sets = [set(ad.var_names) for ad in adatas]
        shared_genes = all_gene_sets[0]
        for s in all_gene_sets[1:]:
            shared_genes &= s
        shared_list = sorted(shared_genes)
        print(f"  Shared genes (all): {len(shared_list)}")

    return shared_list[:n_top]


# ---------------------------------------------------------------------------
# PCA + Harmony + UMAP integration
# ---------------------------------------------------------------------------


def run_integration(
    adatas: list[sc.AnnData],
    dataset_ids: list[str],
    n_top_genes: int = 3000,
    n_comps_pca: int = 50,
) -> sc.AnnData:
    """
    Run the full integration pipeline:
      1. Concatenate with inner join on shared genes
      2. Normalize + log1p
      3. HVG selection on concatenated data
      4. Scale
      5. PCA
      6. Harmony integration
      7. UMAP
    """
    # ── 1. Concatenate ────────────────────────────────────────────────
    print("  Concatenating datasets (inner join on shared genes) ...")
    t0 = time.time()
    concat = sc.AnnData.concatenate(
        *adatas,
        batch_key="dataset_id",
        batch_categories=dataset_ids,
        join="inner",
        index_unique="_",
    )
    print(f"    → {concat.n_obs} cells × {concat.n_vars} genes ({time.time() - t0:.1f}s)")

    # ── 2. Normalize + log1p ──────────────────────────────────────────
    print("  Normalizing ...")
    t1 = time.time()
    sc.pp.normalize_total(concat, target_sum=1e4)
    sc.pp.log1p(concat)
    print(f"    Done ({time.time() - t1:.1f}s)")

    # ── 3. HVG selection on concatenated ──────────────────────────────
    print(f"  Selecting top {n_top_genes} HVGs on concatenated data ...")
    t2 = time.time()
    n_hvg = min(n_top_genes, concat.n_vars - 1)
    sc.pp.highly_variable_genes(concat, n_top_genes=n_hvg)
    concat = concat[:, concat.var["highly_variable"]].copy()
    print(f"    → {concat.n_vars} genes retained ({time.time() - t2:.1f}s)")

    # ── 4. Scale ──────────────────────────────────────────────────────
    print("  Scaling ...")
    t3 = time.time()
    sc.pp.scale(concat, max_value=10)
    print(f"    Done ({time.time() - t3:.1f}s)")

    # ── 5. PCA ────────────────────────────────────────────────────────
    print(f"  Running PCA ({n_comps_pca} components) ...")
    t4 = time.time()
    sc.tl.pca(concat, n_comps=n_comps_pca, svd_solver="arpack")
    print(f"    Done ({time.time() - t4:.1f}s)")

    # ── 6. Harmony integration ────────────────────────────────────────
    print("  Running Harmony integration ...")
    t5 = time.time()
    try:
        sc.external.pp.harmony_integrate(concat, key="dataset_id")
        harmony_rep = "X_pca_harmony"
        print(f"    Done ({time.time() - t5:.1f}s)")
    except Exception as exc:
        print(f"    [WARN] harmony_integrate failed ({exc}) — trying harmonypy fallback ...")
        try:
            import harmonypy

            pca_df = pd.DataFrame(
                concat.obsm["X_pca"],
                index=concat.obs_names,
            )
            meta_df = concat.obs[["dataset_id"]]
            ho = harmonypy.run_harmony(pca_df, meta_df, "dataset_id")
            concat.obsm["X_pca_harmony"] = ho.Z_corr.T
            print(f"    Done ({time.time() - t5:.1f}s)")
            harmony_rep = "X_pca_harmony"
        except ImportError:
            print("    [WARN] harmonypy not available — falling back to PCA (no Harmony)")
            harmony_rep = "X_pca"

    # ── 7. UMAP on Harmony-corrected PCA ──────────────────────────────
    print("  Running UMAP ...")
    t6 = time.time()
    sc.pp.neighbors(concat, use_rep=harmony_rep, n_neighbors=30)
    sc.tl.umap(concat, min_dist=0.3)
    print(f"    Done ({time.time() - t6:.1f}s)")

    return concat


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------


def write_outputs(concat: sc.AnnData, output_dir: str) -> dict:
    """Write integrated parquet files and info JSON.

    Returns the info dict.
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n  Writing outputs to {output_dir} ...")

    # ── UMAP + metadata ───────────────────────────────────────────────
    t1 = time.time()
    umap_coords = concat.obsm["X_umap"]
    n_dims = umap_coords.shape[1]

    meta_df = pd.DataFrame(index=concat.obs_names)
    for i in range(n_dims):
        meta_df[f"UMAP_{i + 1}"] = umap_coords[:, i].astype(np.float32)

    # Copy all obs columns
    for col in concat.obs.columns:
        s = concat.obs[col]
        if s.dtype.name == "category":
            meta_df[col] = s
        elif pd.api.types.is_string_dtype(s) or pd.api.types.is_object_dtype(s):
            meta_df[col] = s.astype("category")
        elif pd.api.types.is_float_dtype(s):
            meta_df[col] = s.astype(np.float32)
        else:
            meta_df[col] = s

    # Add unified cell type column
    cell_type_col = _detect_cell_type_column(concat)
    if cell_type_col and cell_type_col in meta_df.columns:
        meta_df["unified_cell_type"] = meta_df[cell_type_col].astype(str).apply(_unify_cell_type).astype("category")

    meta_path = os.path.join(output_dir, "integrated_umap.parquet")
    meta_df.to_parquet(meta_path)
    print(f"    → integrated_umap.parquet  ({meta_df.shape[0]} rows × {meta_df.shape[1]} cols, {time.time() - t1:.1f}s)")

    # ── Shared gene expression (top 1000 HVGs, float16) ──────────────
    t2 = time.time()
    hvg_mask = concat.var["highly_variable"].values if "highly_variable" in concat.var.columns else slice(None)
    hvg_indices = np.where(hvg_mask)[0]
    n_expr = min(len(hvg_indices), 1000)
    expr_indices = hvg_indices[:n_expr]
    expr_genes = concat.var_names[expr_indices].tolist()

    X_subset = concat[:, expr_indices].X
    if hasattr(X_subset, "toarray"):
        X_subset = X_subset.toarray()

    expr_df = pd.DataFrame(
        X_subset.astype(np.float16),
        index=concat.obs_names,
        columns=expr_genes,
    )
    expr_df.index.name = "cell"
    expr_df["dataset_id"] = concat.obs["dataset_id"].values

    expr_path = os.path.join(output_dir, "integrated_expr.parquet")
    expr_df.to_parquet(expr_path)
    print(f"    → integrated_expr.parquet  ({expr_df.shape[0]} rows × {expr_df.shape[1]} cols, float16, {time.time() - t2:.1f}s)")

    # ── Info JSON ─────────────────────────────────────────────────────
    n_cell_types = int(meta_df["unified_cell_type"].nunique()) if "unified_cell_type" in meta_df.columns else 0

    info = {
        "n_cells": int(concat.n_obs),
        "n_genes": int(concat.n_vars),
        "n_shared_hvgs": int(concat.var["highly_variable"].sum()) if "highly_variable" in concat.var.columns else int(concat.n_vars),
        "n_datasets": int(concat.obs["dataset_id"].nunique()),
        "datasets": sorted(concat.obs["dataset_id"].unique().tolist()),
        "n_unified_cell_types": n_cell_types,
        "cell_type_column": cell_type_col,
        "available_columns": sorted(
            c for c in meta_df.columns if not c.startswith("UMAP_") and c != "dataset_id"
        ),
        "files": {
            "umap": "integrated_umap.parquet",
            "expression": "integrated_expr.parquet",
        },
    }

    info_path = os.path.join(output_dir, "integrated_info.json")
    with open(info_path, "w") as f:
        json.dump(info, f, indent=2, default=str)
    print(f"    → integrated_info.json")

    return info


def _detect_cell_type_column(adata: sc.AnnData) -> str | None:
    """Find the best cell-type annotation column in the concatenated object."""
    for candidate in ("cell_type", "predicted_celltype", "Cell_type", "celltype"):
        if candidate in adata.obs.columns:
            return candidate
    return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Integrate all datasets into shared UMAP space via Harmony."
    )
    parser.add_argument(
        "--data-dir",
        default=".",
        help="Root directory containing GSE* subdirectories (default: current dir)",
    )
    parser.add_argument(
        "--output-dir",
        default="vis_website/data/integrated",
        help="Output directory for integrated parquet files (default: vis_website/data/integrated)",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=3000,
        help="Number of top highly variable genes to use (default: 3000)",
    )
    parser.add_argument(
        "--n-comps-pca",
        type=int,
        default=50,
        help="Number of PCA components (default: 50)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    data_dir = os.path.abspath(args.data_dir)
    output_dir = os.path.abspath(args.output_dir)

    print(f"Data dir:   {data_dir}")
    print(f"Output dir: {output_dir}")
    print()

    # ── 1. Find dataset paths ─────────────────────────────────────────────
    dataset_paths = find_dataset_paths(data_dir)

    if not dataset_paths:
        print("[ERROR] No 04_clustered.h5ad files found.")
        sys.exit(1)

    print(f"Found {len(dataset_paths)} dataset(s):\n")
    for gse, path in dataset_paths.items():
        fsize = os.path.getsize(path) / (1024**3)
        print(f"  {gse}: {path} ({fsize:.2f} GB)")

    # ── 2. Load all datasets ──────────────────────────────────────────────
    print("\nLoading datasets ...")
    adatas: list[sc.AnnData] = []
    dataset_ids: list[str] = []
    for gse_id, h5ad_path in dataset_paths.items():
        t0 = time.time()
        ad = sc.read_h5ad(h5ad_path)
        print(f"  {gse_id}: {ad.n_obs} cells × {ad.n_vars} genes ({time.time() - t0:.1f}s)")
        adatas.append(ad)
        dataset_ids.append(gse_id)

    total_cells = sum(a.n_obs for a in adatas)
    print(f"\nTotal: {total_cells:,} cells across {len(adatas)} datasets")

    # ── 3. Run integration ────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Integration pipeline")
    print(f"{'=' * 60}")

    t_start = time.time()
    concat = run_integration(
        adatas,
        dataset_ids,
        n_top_genes=args.n_top_genes,
        n_comps_pca=args.n_comps_pca,
    )

    # ── 4. Write outputs ──────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Writing outputs")
    print(f"{'=' * 60}")
    info = write_outputs(concat, output_dir)

    elapsed = time.time() - t_start
    print(f"\n{'=' * 60}")
    print(f"[DONE] Integration complete in {elapsed:.1f}s")
    print(f"  {info['n_cells']:,} cells, {info['n_genes']} genes, {info['n_datasets']} datasets")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()

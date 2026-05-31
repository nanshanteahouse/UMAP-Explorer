#!/usr/bin/env python3
"""
build_data_cache.py - Extract data from 04_clustered.h5ad files into parquet format.

For each dataset:
  - umap_metadata.parquet  : UMAP coordinates + all obs metadata + dataset_id
  - gene_expression.parquet: top 1000 highly variable genes expression (float16) + dataset_id
  - dataset_info.json      : summary statistics and metadata

Aggregated:
  - dataset_registry.json  : combined info from all datasets

Usage:
  python build_data_cache.py \
    --data-dir /path/to/neurobiology \
    --output-dir /path/to/vis_website/data
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scanpy as sc

# ---------------------------------------------------------------------------
# Dataset metadata lookup (manually curated from GEO entries)
# ---------------------------------------------------------------------------
DATASET_METADATA = {
    "GSE107618": {
        "species": "Human",
        "tissue": "Retina",
        "citation": "Ray TA et al. (2020) Spatiotemporal patterns of retinal progenitor cell proliferation.",
    },
    "GSE116106": {
        "species": "Human",
        "tissue": "Retina",
        "citation": "Lu Y et al. (2020) Single-cell analysis of human retina identifies evolutionarily conserved and species-specific mechanisms controlling development.",
    },
    "GSE118852": {
        "species": "Macaque",
        "tissue": "Retina",
        "citation": "Clark BS et al. (2019) Single-cell RNA-seq of macaque retinal development.",
    },
    "GSE122970": {
        "species": "Human",
        "tissue": "Retina",
        "citation": "Hoang T et al. (2020) Gene regulatory networks controlling vertebrate retinal regeneration.",
    },
    "GSE138002": {
        "species": "Human",
        "tissue": "Retina",
        "citation": "Sridhar A et al. (2021) Single-cell transcriptomic analysis of the developing human retina.",
    },
    "GSE140877": {
        "species": "Human",
        "tissue": "Brain",
        "citation": "Yan Y et al. (2020) Single-cell dissection of the human brain vasculature.",
    },
    "GSE226108": {
        "species": "Primate",
        "tissue": "Retina+Neural",
        "citation": "Sridhar A et al. (2023) Single-cell transcriptomic analysis of primate retinal and neural development.",
    },
    "GSE246169": {
        "species": "Human",
        "tissue": "Retina",
        "citation": "Workman MJ et al. (2024) Single-cell characterization of human retinal development.",
    },
    "GSE249004": {
        "species": "Marmoset",
        "tissue": "Retina",
        "citation": "Sridhar A et al. (2024) Transcriptomic profiling of marmoset retina at single-cell resolution.",
    },
    "GSE268630": {
        "species": "Human",
        "tissue": "Retina",
        "citation": "Jaffe AE et al. (2024) Single-cell dissected transcriptomic architecture of the human retina multiome.",
    },
}

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def find_dataset_paths(data_dir):
    """
    Scan for all 04_clustered.h5ad files under data_dir.
    Returns dict: {gse_id: absolute_path}
    """
    dataset_paths = {}

    # Known path patterns to try for each GSE directory
    patterns = [
        # Standard: GSE*/results/h5ad/04_clustered.h5ad
        lambda gse: os.path.join(gse, "results", "h5ad", "04_clustered.h5ad"),
        # Nested: GSE*/GSE*_results/h5ad/04_clustered.h5ad
        lambda gse: os.path.join(gse, f"{gse}_results", "h5ad", "04_clustered.h5ad"),
        # Pipeline: GSE*/pipeline/results/h5ad/04_clustered.h5ad
        lambda gse: os.path.join(gse, "pipeline", "results", "h5ad", "04_clustered.h5ad"),
    ]

    if not os.path.isdir(data_dir):
        print(f"[ERROR] data_dir does not exist: {data_dir}")
        return dataset_paths

    # Iterate over GSE directories
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
# Obs column processing
# ---------------------------------------------------------------------------

def _categorize_series(s: pd.Series) -> pd.Series:
    """Convert string/object columns to pandas Categorical for efficient storage."""
    if s.dtype.name == "category":
        return s
    if pd.api.types.is_string_dtype(s) or pd.api.types.is_object_dtype(s):
        return s.astype("category")
    return s


def _downcast_numeric(s: pd.Series) -> pd.Series:
    """Downcast float64 columns to float32 for parquet efficiency."""
    if pd.api.types.is_float_dtype(s):
        return s.astype(np.float32)
    if pd.api.types.is_integer_dtype(s):
        if s.dtype == np.int64:
            min_val, max_val = s.min(), s.max()
            if min_val >= 0:
                if max_val <= 255:
                    return s.astype(np.uint8)
                elif max_val <= 65535:
                    return s.astype(np.uint16)
                elif max_val <= 4294967295:
                    return s.astype(np.uint32)
            else:
                if min_val >= -128 and max_val <= 127:
                    return s.astype(np.int8)
                elif min_val >= -32768 and max_val <= 32767:
                    return s.astype(np.int16)
                elif min_val >= -2147483648 and max_val <= 2147483647:
                    return s.astype(np.int32)
    return s


def process_obs(adata, dataset_id: str):
    """
    Build umap_metadata DataFrame:
      - UMAP_1, UMAP_2 (and UMAP_3 if 3D)
      - All obs columns (categorical-encoded, downcasted)
      - dataset_id column
    """
    umap_coords = adata.obsm["X_umap"]
    n_umap_dims = umap_coords.shape[1]

    df = pd.DataFrame(index=adata.obs_names)
    for i in range(n_umap_dims):
        df[f"UMAP_{i + 1}"] = umap_coords[:, i].astype(np.float32)

    for col in adata.obs.columns:
        s = adata.obs[col]
        s = _categorize_series(s)
        s = _downcast_numeric(s)
        df[col] = s

    df["dataset_id"] = dataset_id

    return df


# ---------------------------------------------------------------------------
# Gene expression processing
# ---------------------------------------------------------------------------

def extract_top_hvg_expression(adata, n_top=1000, dataset_id: str = ""):
    """
    Extract expression of top n_top highly variable genes.
    Returns (DataFrame, list_of_gene_names).
    """
    if "highly_variable" not in adata.var.columns:
        print(f"    [WARN] No 'highly_variable' column found — skipping gene expression")
        return None, []

    hvg = adata.var["highly_variable"].values
    hvg_indices = np.where(hvg)[0]

    if len(hvg_indices) == 0:
        print(f"    [WARN] No highly variable genes marked — skipping gene expression")
        return None, []

    n_available = min(len(hvg_indices), n_top)
    top_indices = hvg_indices[:n_available]
    gene_names = adata.var_names[top_indices].tolist()

    X_subset = adata[:, top_indices].X
    if hasattr(X_subset, "toarray"):
        X_subset = X_subset.toarray()

    X_float16 = X_subset.astype(np.float16)

    df = pd.DataFrame(X_float16, index=adata.obs_names, columns=gene_names)
    df.index.name = "cell"

    if dataset_id:
        df["dataset_id"] = dataset_id

    return df, gene_names


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def process_dataset(h5ad_path: str, output_dir: str, dataset_id: str):
    """Process a single dataset: extract parquet files and info JSON."""
    print(f"\n{'=' * 60}")
    print(f"[PROCESS] {dataset_id}")
    print(f"  Source: {h5ad_path}")
    print(f"  Output: {output_dir}")
    print(f"{'=' * 60}")

    t0 = time.time()

    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes ({time.time() - t0:.1f}s)")

    os.makedirs(output_dir, exist_ok=True)

    # ---- 1. UMAP + metadata ----
    t1 = time.time()
    print(f"  Building umap_metadata...")
    meta_df = process_obs(adata, dataset_id)

    meta_table = pa.Table.from_pandas(meta_df, preserve_index=False)
    meta_path = os.path.join(output_dir, "umap_metadata.parquet")
    pq.write_table(meta_table, meta_path)
    print(f"    → {meta_path}  ({meta_df.shape[0]} rows × {meta_df.shape[1]} cols, {time.time() - t1:.1f}s)")

    # ---- 2. Gene expression ----
    t2 = time.time()
    print(f"  Extracting top 1000 HVG expression...")
    expr_df, hvg_genes = extract_top_hvg_expression(adata, n_top=1000, dataset_id=dataset_id)

    if expr_df is not None:
        expr_table = pa.Table.from_pandas(expr_df, preserve_index=False)
        expr_path = os.path.join(output_dir, "gene_expression.parquet")
        pq.write_table(expr_table, expr_path)
        print(f"    → {expr_path}  ({expr_df.shape[0]} rows × {expr_df.shape[1]} cols, float16, {time.time() - t2:.1f}s)")
    else:
        expr_path = None
        print(f"    [SKIP] No gene expression parquet generated")

    # ---- 3. Dataset info ----
    meta_info = DATASET_METADATA.get(dataset_id, {})
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    available_columns = sorted(adata.obs.columns.tolist())
    obs_column_types = {col: str(adata.obs[col].dtype) for col in adata.obs.columns}

    cell_type_col = None
    for candidate in ("cell_type", "predicted_celltype", "Cell_type", "celltype"):
        if candidate in adata.obs.columns:
            cell_type_col = candidate
            break

    info = {
        "dataset_id": dataset_id,
        "n_cells": int(n_cells),
        "n_genes": int(n_genes),
        "n_hvgs_available": int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var.columns else 0,
        "n_hvgs_extracted": len(hvg_genes) if hvg_genes else 0,
        "hvg_genes": hvg_genes if hvg_genes else [],
        "available_columns": available_columns,
        "column_types": obs_column_types,
        "cell_type_column": cell_type_col,
        "species": meta_info.get("species", "Unknown"),
        "tissue": meta_info.get("tissue", "Unknown"),
        "citation": meta_info.get("citation", ""),
        "umap_dims": int(adata.obsm["X_umap"].shape[1]),
        "files": {
            "umap_metadata": "umap_metadata.parquet",
        },
    }
    if expr_path is not None:
        info["files"]["gene_expression"] = "gene_expression.parquet"

    info_path = os.path.join(output_dir, "dataset_info.json")
    with open(info_path, "w") as f:
        json.dump(info, f, indent=2, default=str)
    print(f"    → {info_path}")
    print(f"  [DONE] {dataset_id} in {time.time() - t0:.1f}s")

    return info


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build data cache: extract 04_clustered.h5ad → parquet + JSON."
    )
    parser.add_argument(
        "--data-dir",
        default=".",
        help="Root directory containing GSE* subdirectories (default: current dir)",
    )
    parser.add_argument(
        "--output-dir",
        default="vis_website/data",
        help="Output directory for parquet files and registry (default: vis_website/data)",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip datasets that already have umap_metadata.parquet in output",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    data_dir = os.path.abspath(args.data_dir)
    output_dir = os.path.abspath(args.output_dir)

    print(f"Data dir:   {data_dir}")
    print(f"Output dir: {output_dir}")
    print(f"Skip existing: {args.skip_existing}")
    print()

    # Find all dataset paths
    dataset_paths = find_dataset_paths(data_dir)

    if not dataset_paths:
        print("[ERROR] No 04_clustered.h5ad files found.")
        sys.exit(1)

    print(f"Found {len(dataset_paths)} dataset(s) with 04_clustered.h5ad:\n")
    for gse, path in dataset_paths.items():
        fsize = os.path.getsize(path) / (1024**3)
        print(f"  {gse}: {path} ({fsize:.2f} GB)")

    print()

    # Process each dataset
    registry_entries = []
    for gse_id, h5ad_path in dataset_paths.items():
        ds_output_dir = os.path.join(output_dir, gse_id)

        # Check if already processed
        if args.skip_existing:
            existing_meta = os.path.join(ds_output_dir, "umap_metadata.parquet")
            existing_info = os.path.join(ds_output_dir, "dataset_info.json")
            if os.path.isfile(existing_meta) and os.path.isfile(existing_info):
                print(f"[SKIP]  {gse_id}: already processed (use --skip-existing to re-run)")
                # Still load existing info for registry
                try:
                    with open(existing_info) as f:
                        registry_entries.append(json.load(f))
                except Exception:
                    pass
                continue

        info = process_dataset(h5ad_path, ds_output_dir, gse_id)
        registry_entries.append(info)

    # Write aggregated registry
    os.makedirs(output_dir, exist_ok=True)
    registry_path = os.path.join(output_dir, "dataset_registry.json")
    with open(registry_path, "w") as f:
        json.dump(registry_entries, f, indent=2, default=str)
    print(f"\n{'=' * 60}")
    print(f"[REGISTRY] {registry_path} ({len(registry_entries)} datasets)")
    print(f"{'=' * 60}")

    # Summary
    total_cells = sum(e["n_cells"] for e in registry_entries)
    total_genes = sum(e["n_genes"] for e in registry_entries)
    print(f"\nSummary: {len(registry_entries)} datasets, {total_cells:,} total cells, {total_genes:,} total genes")
    print("Done.")


if __name__ == "__main__":
    main()

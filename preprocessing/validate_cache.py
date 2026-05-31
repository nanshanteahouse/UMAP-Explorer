#!/usr/bin/env python3
"""
validate_cache.py — Validate parquet data cache produced by build_data_cache.py.

For each dataset directory, validates:
  - umap_metadata.parquet exists, non-empty, has dataset_id column
  - gene_expression.parquet exists with matching row count
  - UMAP coordinates are reasonable: no NaN, within [-50, 50], not all identical
  - dataset_registry.json entry matches actual row count
  - Metadata column coverage (non-null percentage)

Usage:
    python validate_cache.py --cache-dir vis_website/data/
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

UMAP_COL_PREFIXES = ("UMAP", "X_umap", "X_UMAP")


def _find_umap_columns(df: pd.DataFrame) -> list[str]:
    """Return column names that look like UMAP coordinates."""
    return [
        c for c in df.columns
        if c.upper().startswith("UMAP_") or c.upper().startswith("X_UMAP_")
    ]


def _print_coverage(coverage: dict) -> None:
    """Print a formatted coverage report sorted by ascending coverage."""
    if not coverage:
        return
    print("  Metadata coverage (non-null %):")
    sorted_cols = sorted(coverage.items(), key=lambda x: x[1]["pct"])
    for col, cov in sorted_cols:
        pct = cov["pct"]
        bar_len = max(1, int(pct / 5))
        bar = "█" * bar_len + "░" * (20 - bar_len)
        print(
            f"    {col:30s} {pct:6.1f}% {bar} "
            f"({cov['non_null']}/{cov['total']})"
        )


# ---------------------------------------------------------------------------
# Per-dataset validation
# ---------------------------------------------------------------------------

def validate_dataset(
    dataset_id: str,
    dataset_dir: Path,
    registry_entry: dict | None,
) -> tuple[bool, list[str], dict]:
    """
    Validate a single dataset directory.

    Returns
    -------
    (passed, issues, coverage)
        coverage maps column name → {non_null, total, pct}
    """
    issues: list[str] = []
    coverage: dict = {}

    umap_path = dataset_dir / "umap_metadata.parquet"
    expr_path = dataset_dir / "gene_expression.parquet"

    n_cells_meta = 0
    meta_df: pd.DataFrame | None = None

    # ---- 1. umap_metadata.parquet ----
    if not umap_path.exists():
        issues.append("umap_metadata.parquet missing")
    else:
        try:
            meta_df = pd.read_parquet(umap_path)
        except Exception as exc:
            issues.append(f"umap_metadata.parquet read error: {exc}")

        if meta_df is not None:
            n_cells_meta = len(meta_df)

            if n_cells_meta == 0:
                issues.append("umap_metadata.parquet has 0 rows")
            else:
                # ---- 2. dataset_id column ----
                if "dataset_id" not in meta_df.columns:
                    issues.append("umap_metadata missing 'dataset_id' column")

                # ---- 3. UMAP coordinate checks ----
                umap_cols = _find_umap_columns(meta_df)
                if not umap_cols:
                    issues.append(
                        "No UMAP coordinate columns found "
                        "(expected columns starting with 'UMAP' or 'X_umap')"
                    )
                else:
                    umap_data = meta_df[umap_cols]

                    nan_count = int(umap_data.isnull().sum().sum())
                    if nan_count > 0:
                        issues.append(
                            f"UMAP coordinates contain {nan_count} NaN value(s)"
                        )

                    valid_umap = umap_data.dropna()
                    if len(valid_umap) > 0:
                        for col in umap_cols:
                            col_data = valid_umap[col]
                            if not pd.api.types.is_numeric_dtype(col_data):
                                issues.append(
                                    f"UMAP column '{col}' has non-numeric "
                                    f"dtype ({col_data.dtype})"
                                )
                                continue
                            lo = float(col_data.min())
                            hi = float(col_data.max())
                            if lo < -50.0 or hi > 50.0:
                                issues.append(
                                    f"UMAP column '{col}' range "
                                    f"[{lo:.2f}, {hi:.2f}] exceeds [-50, 50]"
                                )

                        for col in umap_cols:
                            col_data = valid_umap[col]
                            if not pd.api.types.is_numeric_dtype(col_data):
                                continue
                            nunique = int(col_data.nunique())
                            if nunique <= 1:
                                val_str = f"{float(col_data.iloc[0]):.4f}"
                                issues.append(
                                    f"UMAP column '{col}' has all identical "
                                    f"values ({val_str})"
                                )

                # ---- 4. Metadata column coverage ----
                for col in meta_df.columns:
                    # skip UMAP coordinate columns
                    if col.upper().startswith("UMAP_") or col.upper().startswith("X_UMAP_"):
                        continue
                    n_non_null = int(meta_df[col].notna().sum())
                    pct = float((n_non_null / n_cells_meta) * 100)
                    coverage[col] = {
                        "non_null": n_non_null,
                        "total": n_cells_meta,
                        "pct": round(pct, 1),
                    }

    # ---- 5. gene_expression.parquet ----
    if not expr_path.exists():
        issues.append("gene_expression.parquet missing")
    else:
        try:
            expr_df = pd.read_parquet(expr_path)
        except Exception as exc:
            issues.append(f"gene_expression.parquet read error: {exc}")
            expr_df = None

        if expr_df is not None:
            n_cells_expr = len(expr_df)
            if n_cells_expr == 0:
                issues.append("gene_expression.parquet has 0 rows")
            elif meta_df is not None and n_cells_expr != n_cells_meta:
                issues.append(
                    f"Row count mismatch: umap_metadata has {n_cells_meta} rows, "
                    f"gene_expression has {n_cells_expr} rows"
                )

    # ---- 6. Validate against registry ----
    if registry_entry is not None and meta_df is not None:
        expected_n = registry_entry.get("n_cells")
        if expected_n is not None and expected_n != n_cells_meta:
            issues.append(
                f"registry n_cells={expected_n} does not match "
                f"actual row count {n_cells_meta}"
            )

    passed = len(issues) == 0
    return passed, issues, coverage


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def validate_cache(cache_dir: str) -> int:
    """Run full validation.  Returns 0 on success, 1 on failure."""
    cache_path = Path(cache_dir)

    if not cache_path.exists():
        print(f"ERROR: Cache directory not found: {cache_dir}")
        print("Run build_data_cache.py first to generate the data cache.")
        return 1

    # ---- Load registry ----
    registry: dict = {}
    registry_path = cache_path / "dataset_registry.json"
    if registry_path.exists():
        with open(registry_path) as fh:
            registry = json.load(fh)
        if isinstance(registry, list):
            registry = {
                entry.get("dataset_id", ""): entry
                for entry in registry
                if entry.get("dataset_id")
            }
        print(f"Loaded dataset_registry.json — {len(registry)} entries\n")
    else:
        print("WARNING: dataset_registry.json not found "
              "(registry cross-checks skipped)\n")

    # ---- Discover dataset directories ----
    dataset_dirs = sorted(
        d for d in cache_path.iterdir()
        if d.is_dir() and d.name.upper().startswith("GSE")
    )

    if not dataset_dirs:
        print("No GSE dataset directories found in cache directory.")
        print("Run build_data_cache.py first to generate the data cache.")
        return 1

    print(f"Found {len(dataset_dirs)} dataset directories\n")

    # ---- Validate each dataset ----
    total = len(dataset_dirs)
    passed_count = 0
    results: list[tuple[str, bool, list[str]]] = []

    for d in dataset_dirs:
        dataset_id = d.name
        registry_entry = registry.get(dataset_id)

        print("=" * 60)
        print(f"  Dataset: {dataset_id}")
        print("=" * 60)

        passed, issues, coverage = validate_dataset(
            dataset_id, d, registry_entry,
        )

        if passed:
            passed_count += 1
            print("  ✅ PASS")
        else:
            print("  ❌ FAIL")
            for issue in issues:
                print(f"     • {issue}")

        _print_coverage(coverage)
        print()
        results.append((dataset_id, passed, issues))

    # ---- Overall summary ----
    print("=" * 60)
    print("  OVERALL SUMMARY")
    print("=" * 60)
    print(f"  {passed_count} / {total} datasets passed")

    for dataset_id, passed, issues in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        if issues:
            detail = " — ".join(issues[:2])
            if len(issues) > 2:
                detail += f" (+{len(issues)-2} more)"
            print(f"    {status}  {dataset_id}  ({detail})")
        else:
            print(f"    {status}  {dataset_id}")

    print()

    if passed_count == total:
        print("✓ All datasets validated successfully.")
        return 0
    else:
        print(f"✗ {total - passed_count} dataset(s) failed validation.")
        return 1


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate parquet data cache produced by build_data_cache.py",
    )
    parser.add_argument(
        "--cache-dir",
        default="vis_website/data/",
        help="Path to the data cache directory (default: vis_website/data/)",
    )
    args = parser.parse_args()

    sys.exit(validate_cache(args.cache_dir))


if __name__ == "__main__":
    main()

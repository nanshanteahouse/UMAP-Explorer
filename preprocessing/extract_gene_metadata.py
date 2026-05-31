"""
extract_gene_metadata.py
=======================
Builds a unified gene index across all scRNA-seq datasets for cross-dataset
gene searching and comparison.

Reads all 04_clustered.h5ad files, extracts gene names (var_names) and HVG
status, annotates known retinal marker genes, and computes cross-dataset
overlap statistics.

Output: vis_website/data/gene_index.json
"""

import argparse
import json
import sys
from pathlib import Path

import scanpy as sc


# ---------------------------------------------------------------------------
# Known retinal marker genes (grouped by cell type)
# ---------------------------------------------------------------------------
MARKER_GENES = {
    "RPC":       ["VSX2", "PAX6", "SOX2", "HES1", "NOTCH1"],
    "Neurogenic": ["ASCL1", "NEUROG2", "TUBB3", "DCX", "STMN2"],
    "RGC":       ["POU4F2", "POU4F1", "NEFM", "NEFL", "ELAVL4"],
    "Cone":      ["ARR3", "OPN1SW", "GNAT2", "PDE6C", "RCVRN"],
    "Rod":       ["NRL", "RHO", "GNAT1", "PDE6B", "SAG"],
    "Bipolar":   ["VSX1", "PRKCA", "TRPM1", "GRM6", "CABP5"],
    "Amacrine":  ["GAD1", "GAD2", "SLC6A9", "TFAP2A", "TFAP2B"],
    "Muller":    ["RLBP1", "GFAP", "VIM", "SLC1A3", "GLUL"],
    "Horizontal": ["ONECUT1", "ONECUT2", "PROX1", "CALB1", "LHX1"],
}

# Build reverse mapping: gene_symbol -> list of cell types
_GENE_TO_CELL_TYPES: dict[str, list[str]] = {}
for cell_type, genes in MARKER_GENES.items():
    for g in genes:
        _GENE_TO_CELL_TYPES.setdefault(g, []).append(cell_type)


# ---------------------------------------------------------------------------
# Dataset discovery
# ---------------------------------------------------------------------------

def discover_datasets(data_dir: Path) -> list[tuple[str, Path]]:
    """
    Find all 04_clustered.h5ad files under *data_dir*.

    Returns a list of ``(dataset_id, filepath)`` tuples sorted by dataset ID.

    Handles three known directory layouts:

    * ``GSE<id>/results/h5ad/04_clustered.h5ad``  (standard)
    * ``GSE<id>/GSE<id>_results/h5ad/04_clustered.h5ad``  (nested)
    * ``GSE268630/pipeline/results/h5ad/04_clustered.h5ad``  (pipeline layout)
    """
    datasets: list[tuple[str, Path]] = []

    # Pattern 1 & 2: standard + nested
    for h5ad_path in data_dir.rglob("04_clustered.h5ad"):
        ds_id = _infer_dataset_id(h5ad_path)
        if ds_id:
            datasets.append((ds_id, h5ad_path))

    # Sort for reproducibility
    datasets.sort(key=lambda x: x[0])
    return datasets


def _infer_dataset_id(h5ad_path: Path) -> str | None:
    """
    Extract the GEO accession (e.g. ``GSE107618``) from an h5ad path.

    Returns ``None`` if no GEO-like component can be identified.
    """
    for part in h5ad_path.parts:
        if part.startswith("GSE") and len(part) > 3 and part[3:].isdigit():
            return part
    return None


# ---------------------------------------------------------------------------
# Gene index construction
# ---------------------------------------------------------------------------

def build_gene_index(
    datasets: list[tuple[str, Path]],
) -> dict:
    """
    Read every dataset's ``04_clustered.h5ad`` and assemble a unified gene
    index together with high-level statistics.
    """
    # Per-dataset gene sets
    dataset_genes: dict[str, set[str]] = {}
    dataset_hvg: dict[str, set[str]] = {}
    dataset_info: dict[str, dict] = {}

    for ds_id, h5ad_path in datasets:
        print(f"  [{ds_id}] Reading {h5ad_path} ...")
        adata = sc.read_h5ad(str(h5ad_path), backed="r")

        all_genes: set[str] = set(adata.var_names)
        dataset_genes[ds_id] = all_genes

        # HVG status
        if "highly_variable" in adata.var.columns:
            hvg_genes = set(
                adata.var_names[adata.var["highly_variable"].values]
            )
        else:
            hvg_genes = set()
        dataset_hvg[ds_id] = hvg_genes

        dataset_info[ds_id] = {
            "n_genes": len(all_genes),
            "n_hvg": len(hvg_genes),
        }
        adata.file.close()
        print(
            f"    → {len(all_genes)} genes, "
            f"{len(hvg_genes)} HVGs"
        )

    # ── Unified gene index ──────────────────────────────────────────────
    all_unique_genes: set[str] = set()
    for genes in dataset_genes.values():
        all_unique_genes |= genes

    # Genes shared across ALL datasets
    shared_genes: set[str] = set.intersection(*dataset_genes.values()) if dataset_genes else set()

    # Marker gene sets
    all_marker_symbols: set[str] = set(_GENE_TO_CELL_TYPES.keys())

    # Build per-gene records
    genes_index: dict[str, dict] = {}
    for gene in sorted(all_unique_genes):
        present_in: list[str] = sorted(
            ds for ds, genes in dataset_genes.items() if gene in genes
        )
        hvg_status: dict[str, bool] = {
            ds: (gene in dataset_hvg[ds])
            for ds in present_in
        }

        entry: dict = {
            "datasets": present_in,
            "is_hvg": hvg_status,
            "is_marker": gene in all_marker_symbols,
        }
        if gene in _GENE_TO_CELL_TYPES:
            entry["cell_types"] = _GENE_TO_CELL_TYPES[gene]

        genes_index[gene] = entry

    # ── Per-dataset marker gene counts ──────────────────────────────────
    marker_found: dict[str, int] = {}
    for ds_id, genes in dataset_genes.items():
        marker_found[ds_id] = len(genes & all_marker_symbols)

    # ── Marker genes found in each dataset ──────────────────────────────
    marker_gene_datasets: dict[str, list[str]] = {}
    for marker in sorted(all_marker_symbols):
        if marker in genes_index:
            marker_gene_datasets[marker] = genes_index[marker]["datasets"]

    # ── Pairwise overlap matrix (Jaccard) ───────────────────────────────
    ds_ids = sorted(dataset_genes.keys())
    jaccard_matrix: dict[str, dict[str, float]] = {}
    for ds_a in ds_ids:
        jaccard_matrix[ds_a] = {}
        genes_a = dataset_genes[ds_a]
        for ds_b in ds_ids:
            genes_b = dataset_genes[ds_b]
            intersection = genes_a & genes_b
            union = genes_a | genes_b
            jaccard = round(len(intersection) / len(union), 4) if union else 0.0
            jaccard_matrix[ds_a][ds_b] = jaccard

    # ── Shared gene intersection counts ─────────────────────────────────
    shared_marker_genes = sorted(
        g for g in shared_genes if g in all_marker_symbols
    )

    # ── Statistics ──────────────────────────────────────────────────────
    statistics = {
        "total_unique_genes": len(all_unique_genes),
        "shared_across_all": sorted(shared_genes),
        "n_shared_across_all": len(shared_genes),
        "per_dataset_gene_counts": {
            ds: {
                "total": dataset_info[ds]["n_genes"],
                "hvg": dataset_info[ds]["n_hvg"],
                "marker_genes_found": marker_found[ds],
            }
            for ds in ds_ids
        },
        "marker_genes": {
            marker: {
                "datasets": ds_list,
                "cell_types": _GENE_TO_CELL_TYPES.get(marker, []),
            }
            for marker, ds_list in marker_gene_datasets.items()
        },
        "shared_marker_genes": shared_marker_genes,
        "n_shared_marker_genes": len(shared_marker_genes),
        "jaccard_overlap": jaccard_matrix,
    }

    return {
        "genes": genes_index,
        "statistics": statistics,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build unified gene index across all scRNA-seq datasets."
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default=".",
        help="Root directory containing GSE* dataset folders (default: .)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="vis_website/data",
        help="Directory to write gene_index.json (default: vis_website/data)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    if not data_dir.is_dir():
        print(f"ERROR: data-dir not found: {data_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")

    # ── Discover datasets ───────────────────────────────────────────────
    print("\nDiscovering datasets ...")
    datasets = discover_datasets(data_dir)
    if not datasets:
        print("ERROR: No 04_clustered.h5ad files found!", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(datasets)} dataset(s):")
    for ds_id, h5ad_path in datasets:
        print(f"  {ds_id:<15s} {h5ad_path.relative_to(data_dir)}")

    # ── Build gene index ────────────────────────────────────────────────
    print("\nBuilding gene index ...")
    gene_index = build_gene_index(datasets)
    stats = gene_index["statistics"]

    # ── Write output ────────────────────────────────────────────────────
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "gene_index.json"
    print(f"\nWriting {output_path} ...")
    with open(output_path, "w") as f:
        json.dump(gene_index, f, indent=2)

    # ── Summary ─────────────────────────────────────────────────────────
    n_markers_found = sum(
        1 for g in gene_index["genes"].values() if g["is_marker"]
    )
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Total unique genes:       {stats['total_unique_genes']}")
    print(f"  Genes shared across all:  {stats['n_shared_across_all']}")
    print(f"  Marker genes indexed:     {n_markers_found}")
    print(f"  Marker genes shared:      {stats['n_shared_marker_genes']}")
    print(f"  Datasets processed:       {len(datasets)}")
    print()
    print("Per-dataset gene counts:")
    for ds_id in sorted(datasets, key=lambda x: x[0]):
        info = stats["per_dataset_gene_counts"][ds_id[0]]
        print(
            f"    {ds_id[0]:<15s}  "
            f"total={info['total']:<6d}  "
            f"hvg={info['hvg']:<6d}  "
            f"markers_found={info['marker_genes_found']}"
        )
    if stats["shared_marker_genes"]:
        print(f"\n  Marker genes in ALL datasets ({stats['n_shared_marker_genes']}):")
        print(f"    {', '.join(stats['shared_marker_genes'])}")
    print()
    print(f"Output: {output_path}")
    print("Done.")


if __name__ == "__main__":
    main()

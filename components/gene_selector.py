"""
gene_selector.py — Gene search & autocomplete for the UMAP explorer.

Provides three utility functions that generate ``dcc.Dropdown``-style option
lists from the cross-dataset gene index (``gene_index.json``).  Used by
``controls.py`` and page-level callbacks to populate the gene-search
dropdown with grouped, labelled options.

Typical usage::

    from components.gene_selector import get_gene_options

    options = get_gene_options(
        dataset_id="GSE107618",
        gene_index=DataLoader().load_gene_index(),
    )
    # → [{"label": "PAX6 ★ (RPC)", "value": "PAX6"}, …]
"""

from __future__ import annotations

from typing import Any, Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# ── Maximum number of dropdown options returned by any function ────────────
# Keeps the UI responsive and avoids overwhelming dcc.Dropdown rendering.
_MAX_OPTIONS = 8000

# ── Display symbols ───────────────────────────────────────────────────────
_MARKER_SYMBOL = "\u2605"  # ★


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def _format_marker_label(gene: str, gene_info: dict[str, Any]) -> str:
    """Build a dropdown label for a marker gene.

    Appends the cell type(s) in parentheses when available, e.g.
    ``"PAX6 ★ (RPC)"``.
    """
    cell_types: list[str] = gene_info.get("cell_types", [])
    if cell_types:
        return f"{gene} {_MARKER_SYMBOL} ({', '.join(cell_types)})"
    return f"{gene} {_MARKER_SYMBOL}"


def _format_shared_label(gene: str, gene_info: dict[str, Any], n_datasets: int) -> str:
    """Build a dropdown label for cross-dataset search results.

    Marker genes get a star suffix; all results show the dataset count, e.g.
    ``"PAX6 ★ (7 datasets)"`` or ``"A2M (9 datasets)"``.
    """
    suffix = f" ({n_datasets} dataset{'s' if n_datasets != 1 else ''})"
    if gene_info.get("is_marker"):
        return f"{gene} {_MARKER_SYMBOL}{suffix}"
    return f"{gene}{suffix}"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_gene_options(
    dataset_id: str,
    gene_index: dict[str, Any],
    available_genes: Optional[set[str]] = None,
    data_loader: object = None,
) -> list[dict[str, str]]:
    """Return ``dcc.Dropdown`` options for genes in *dataset_id*.

    Options are grouped in this order:
      1. Marker genes present in the dataset (sorted alphabetically).
      2. HVG (highly variable) genes that are **not** markers (sorted).
      3. All remaining genes (sorted).

    Each option is a dict with ``label`` and ``value`` keys.  Marker gene
    labels include the star symbol and associated cell type(s), e.g.
    ``"PAX6 ★ (RPC)"``.  Non-marker genes show only the gene name.

    Parameters
    ----------
    dataset_id : str
        GSE identifier, e.g. ``'GSE107618'``.
    gene_index : dict
        The ``gene_index.json`` contents — a dict with a ``"genes"`` key
        mapping gene symbols to their metadata.
    available_genes : set[str] or None
        When provided, only genes present in this set are included in
        the returned options.  Use this to restrict the dropdown to
        genes that actually have expression data in the cached parquet.
    data_loader : object, optional
        Reserved for future use.  Not currently used.

    Returns
    -------
    list[dict[str, str]]
        Up to 200 options, each ``{"label": str, "value": str}``.
    """
    genes: dict[str, dict[str, Any]] = gene_index.get("genes", {})

    markers: list[dict[str, str]] = []
    hvgs: list[dict[str, str]] = []
    others: list[dict[str, str]] = []

    for gene_name, info in genes.items():
        # Skip genes not present in the requested dataset.
        if dataset_id not in info.get("datasets", []):
            continue

        # Filter to only genes with expression data in the parquet cache.
        if available_genes is not None and gene_name not in available_genes:
            continue

        is_marker: bool = info.get("is_marker", False)
        is_hvg: bool = info.get("is_hvg", {}).get(dataset_id, False)

        if is_marker:
            label = _format_marker_label(gene_name, info)
            markers.append({"label": label, "value": gene_name})
        elif is_hvg:
            hvgs.append({"label": gene_name, "value": gene_name})
        else:
            others.append({"label": gene_name, "value": gene_name})

    # Sort each group alphabetically by gene name.
    markers.sort(key=lambda o: o["value"])
    hvgs.sort(key=lambda o: o["value"])
    others.sort(key=lambda o: o["value"])

    # Concatenate groups and cap at MAX_OPTIONS.
    result = markers + hvgs + others
    return result[:_MAX_OPTIONS]


def get_marker_gene_options(gene_index: dict[str, Any]) -> list[dict[str, str]]:
    """Return dropdown options for all known marker genes.

    Only genes whose ``is_marker`` flag is ``True`` in *gene_index* are
    included.  Labels include the star symbol and the associated cell
    type(s), e.g. ``"PAX6 ★ (RPC)"``.

    Parameters
    ----------
    gene_index : dict
        The ``gene_index.json`` contents.

    Returns
    -------
    list[dict[str, str]]
        Up to 200 options sorted alphabetically.
    """
    genes: dict[str, dict[str, Any]] = gene_index.get("genes", {})
    options: list[dict[str, str]] = []

    for gene_name, info in genes.items():
        if not info.get("is_marker", False):
            continue
        label = _format_marker_label(gene_name, info)
        options.append({"label": label, "value": gene_name})

    options.sort(key=lambda o: o["value"])
    return options[:_MAX_OPTIONS]


def search_genes(
    query: str,
    gene_index: dict[str, Any],
    limit: int = 50,
) -> list[dict[str, str]]:
    """Search across *all* datasets for genes matching *query*.

    Performs a **case-insensitive substring** match on gene symbols.
    Results are sorted with marker genes first, then alphabetically.
    Each label includes the number of datasets the gene appears in, e.g.
    ``"PAX6 ★ (7 datasets)"`` or ``"A2M (9 datasets)"``.

    Parameters
    ----------
    query : str
        Substring to search for.  Case-insensitive.
    gene_index : dict
        The ``gene_index.json`` contents.
    limit : int, optional
        Maximum number of results (default 50).  Hard-capped at 200.

    Returns
    -------
    list[dict[str, str]]
        Each ``{"label": str, "value": str}``, sorted markers-first,
        then alphabetically.
    """
    if not query:
        return []

    genes: dict[str, dict[str, Any]] = gene_index.get("genes", {})
    query_lower: str = query.lower()

    markers: list[dict[str, str]] = []
    others: list[dict[str, str]] = []

    for gene_name, info in genes.items():
        if query_lower not in gene_name.lower():
            continue

        n_datasets: int = len(info.get("datasets", []))
        label = _format_shared_label(gene_name, info, n_datasets)
        entry: dict[str, str] = {"label": label, "value": gene_name}

        if info.get("is_marker", False):
            markers.append(entry)
        else:
            others.append(entry)

    markers.sort(key=lambda o: o["value"])
    others.sort(key=lambda o: o["value"])

    result = markers + others
    # Apply the caller's limit but never exceed the hard cap.
    effective_limit = min(limit, _MAX_OPTIONS)
    return result[:effective_limit]

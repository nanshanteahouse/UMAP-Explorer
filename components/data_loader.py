"""
data_loader.py — Shared parquet data loading utility for all Dash components.

This module is the single point of contact with the preprocessed data.
All pages and components use this loader to access dataset metadata,
UMAP coordinates, observation annotations, and gene expression values.

Typical usage:
    from components.data_loader import DataLoader

    loader = DataLoader()
    registry = loader.load_registry()
    meta_df = loader.load_metadata("GSE107618")
    expr_series = loader.load_gene_expression("GSE107618", "A2M")

Uses pyarrow memory mapping via pandas.read_parquet() and caches loaded
DataFrames in a module-level dict for efficient reuse within a session.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Module-level cache
# ---------------------------------------------------------------------------
# Key: (dataset_id, file_type) -> DataFrame
# Cleared when the Python process exits (no explicit flush needed).
_cache: dict[tuple[str, str], pd.DataFrame] = {}

# ---------------------------------------------------------------------------
# Module-level constant: absolute path to the data directory
# ---------------------------------------------------------------------------
# Resolved relative to this file's location so the loader works regardless
# of the process working directory.
DATA_DIR = str(
    Path(__file__).resolve().parent.parent / "data"
)


# ---------------------------------------------------------------------------
# DataLoader
# ---------------------------------------------------------------------------

class DataLoader:
    """Loads parquet + JSON data files from the preprocessed cache.

    All loaded DataFrames are cached in the module-level ``_cache`` dict so
    that repeated calls for the same dataset + file type return quickly.

    Parameters
    ----------
    data_dir : str, optional
        Path to the data directory.  Defaults to ``DATA_DIR``
        (``vis_website/data/`` resolved from this file's location).
    """

    def __init__(self, data_dir: str | None = None) -> None:
        self.data_dir = data_dir if data_dir is not None else DATA_DIR

    # -- Registry (JSON) ---------------------------------------------------

    def load_registry(self) -> list[dict]:
        """Load the dataset registry.

        Returns
        -------
        list[dict]
            List of dataset info dictionaries from
            ``dataset_registry.json``.
        """
        path = os.path.join(self.data_dir, "dataset_registry.json")
        with open(path, "r") as f:
            return json.load(f)

    def load_gene_index(self) -> dict:
        """Load the cross-dataset gene index.

        Returns
        -------
        dict
            The ``gene_index.json`` contents — a mapping of gene symbols to
            their dataset membership and HVG/marker flags.
        """
        path = os.path.join(self.data_dir, "gene_index.json")
        with open(path, "r") as f:
            return json.load(f)

    # -- Per-dataset parquet loaders ---------------------------------------

    def load_metadata(self, dataset_id: str) -> pd.DataFrame:
        """Load UMAP metadata parquet for *dataset_id*.

        Returns a DataFrame with UMAP coordinates (UMAP_1, UMAP_2, …),
        all observation annotation columns, and the ``dataset_id`` column.

        UMAP coordinate columns are cast to ``float32`` to save memory.

        Parameters
        ----------
        dataset_id : str
            GSE identifier, e.g. ``'GSE107618'``.

        Returns
        -------
        pd.DataFrame
        """
        key = (dataset_id, "umap_metadata")
        if key in _cache:
            return _cache[key]

        path = os.path.join(self.data_dir, dataset_id, "umap_metadata.parquet")
        df: pd.DataFrame = pd.read_parquet(path)

        # Downcast UMAP coordinates to float32 for memory efficiency.
        umap_cols = [c for c in df.columns if c.startswith("UMAP_")]
        for col in umap_cols:
            df[col] = df[col].astype("float32")

        _cache[key] = df
        return df

    def load_gene_expression(
        self,
        dataset_id: str,
        gene_name: str | None = None,
    ) -> pd.DataFrame | pd.Series:
        """Load gene expression data for *dataset_id*.

        Parameters
        ----------
        dataset_id : str
            GSE identifier.
        gene_name : str, optional
            If provided, returns a **Series** for that single gene.
            If ``None`` (default), returns the full expression DataFrame
            (all gene columns plus ``dataset_id``).

        Returns
        -------
        pd.DataFrame or pd.Series
            Series when *gene_name* is given; DataFrame otherwise.
        """
        key = (dataset_id, "gene_expression")
        if key in _cache:
            df = _cache[key]
        else:
            path = os.path.join(
                self.data_dir, dataset_id, "gene_expression.parquet"
            )
            df = pd.read_parquet(path)
            _cache[key] = df

        if gene_name is not None:
            return df[gene_name]
        return df

    # -- Metadata helpers --------------------------------------------------

    def get_available_columns(self, dataset_id: str) -> list[str]:
        """List metadata column names for *dataset_id*.

        Excludes UMAP coordinate columns (``UMAP_*``) and the
        ``dataset_id`` column.

        Parameters
        ----------
        dataset_id : str
            GSE identifier.

        Returns
        -------
        list[str]
            Column names suitable for colouring / selection in the UI.
        """
        df = self.load_metadata(dataset_id)
        return [
            c
            for c in df.columns
            if not c.startswith("UMAP_") and c != "dataset_id"
        ]

    def get_available_genes(self, dataset_id: str) -> list[str]:
        """List gene names available in *dataset_id*'s expression data.

        The ``dataset_id`` column is excluded from the result.

        Parameters
        ----------
        dataset_id : str
            GSE identifier.

        Returns
        -------
        list[str]
            Sorted gene names.
        """
        df = self.load_gene_expression(dataset_id)
        return sorted([c for c in df.columns if c != "dataset_id"])

    def get_column_type(
        self, dataset_id: str, column_name: str
    ) -> str:
        """Return the display type of *column_name* in the metadata.

        Used by the UI to decide whether to colour cells with a discrete
        palette (categorical), a continuous colour scale (numeric), or a
        binary two-colour scheme (bool).

        Parameters
        ----------
        dataset_id : str
            GSE identifier.
        column_name : str
            Name of the column to inspect.

        Returns
        -------
        str
            ``'categorical'``, ``'numeric'``, or ``'bool'``.
        """
        df = self.load_metadata(dataset_id)
        dtype = df[column_name].dtype

        if pd.api.types.is_bool_dtype(dtype):
            return "bool"
        if pd.api.types.is_categorical_dtype(dtype) or pd.api.types.is_object_dtype(dtype):
            return "categorical"
        if pd.api.types.is_numeric_dtype(dtype):
            return "numeric"

        # Fallback — should not normally be reached.
        return "categorical"


# ---------------------------------------------------------------------------
# Module-level convenience functions
# ---------------------------------------------------------------------------
# These let callers avoid instantiating DataLoader when they just need
# one-off access.  A single default loader is lazily created and reused.

_default_loader: DataLoader | None = None


def _get_default() -> DataLoader:
    global _default_loader
    if _default_loader is None:
        _default_loader = DataLoader()
    return _default_loader


def load_registry() -> list[dict]:
    """Convenience — call ``DataLoader().load_registry()``."""
    return _get_default().load_registry()


def load_gene_index() -> dict:
    """Convenience — call ``DataLoader().load_gene_index()``."""
    return _get_default().load_gene_index()

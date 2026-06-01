"""
umap_figure.py — Core UMAP rendering function for WebGL-accelerated visualization.

Provides a single entry point ``make_umap_figure`` that produces a
``plotly.graph_objects.Figure`` with a ``go.Scattergl`` (2D) or
``go.Scatter3d`` (3D) trace.  The function is designed to handle up to
500k cells via adaptive random sampling and WebGL rendering.

Typical usage (from a Dash callback)::

    from components.umap_figure import make_umap_figure

    fig = make_umap_figure(
        df=metadata_df,
        color_col="cell_type",
        dim="2d",
        color_type="categorical",
    )
    return fig
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import qualitative

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Default marker styling applied to all points.
_MARKER_SIZE = 3
_MARKER_OPACITY = 0.7
_MARKER_LINE_WIDTH = 0

# Hover fields in priority order — the function checks which columns actually
# exist in the DataFrame and includes only those.
_HOVER_FIELDS = ["cell_type", "stage", "sample"]

# User-facing label for each known hover field (title-cased fallback for
# unknown columns).
_HOVER_LABELS: dict[str, str] = {
    "cell_type": "Cell type",
    "stage": "Stage",
    "sample": "Sample",
    "tissue": "Tissue",
    "cell_type_sub": "Cell type (sub)",
    "predicted_celltype": "Predicted cell type",
}

# Default sequential colours for categorical colouring.
_CATEGORICAL_PALETTE = qualitative.Plotly

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def make_umap_figure(
    df: pd.DataFrame,
    color_col: str | None = None,
    dim: str = "2d",
    sample_size: int = 50000,
    color_type: str = "categorical",
    gene_series: pd.Series | None = None,
) -> go.Figure:
    """Build a WebGL-accelerated UMAP scatter plot.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with at least ``UMAP_1`` and ``UMAP_2`` columns (and
        optionally ``UMAP_3`` for 3D views), plus metadata columns for
        hover and colouring.
    color_col : str or None
        Column in *df* to colour points by.  When ``None``, all points
        are rendered in the default colour.
    dim : {'2d', '3d'}
        Dimensionality of the UMAP projection.
    sample_size : int
        Maximum number of cells to render.  When the DataFrame has fewer
        rows than *sample_size*, all cells are plotted (no subsampling).
    color_type : {'categorical', 'numeric', 'gene'}
        - ``'categorical'`` — discrete colouring using the
          :data:`_CATEGORICAL_PALETTE` (``qualitative.Plotly``).
        - ``'numeric'`` — continuous colouring using the *color_col*
          column as numeric values with a Viridis colour scale.
        - ``'gene'`` — continuous colouring using *gene_series*
          expression values with a Viridis colour scale.
    gene_series : pd.Series or None
        Gene expression values indexed by cell barcode.  Required when
        *color_type* is ``'gene'``; ignored otherwise.

    Returns
    -------
    go.Figure
        Figure with a single trace (``go.Scattergl`` for 2D,
        ``go.Scatter3d`` for 3D).
    """
    # ------------------------------------------------------------------
    # 1. Adaptive sampling
    # ------------------------------------------------------------------
    n_plot = min(len(df), sample_size)
    df_plot = df.sample(n=n_plot, random_state=42) if n_plot < len(df) else df

    # ------------------------------------------------------------------
    # 2. UMAP coordinate arrays
    # ------------------------------------------------------------------
    x: np.ndarray = df_plot["UMAP_1"].values
    y: np.ndarray = df_plot["UMAP_2"].values

    is_3d = dim == "3d" and "UMAP_3" in df_plot.columns
    z: np.ndarray | None = df_plot["UMAP_3"].values if is_3d else None

    # ------------------------------------------------------------------
    # 3. Hover customdata and template
    # ------------------------------------------------------------------
    hover_cols: list[str] = [c for c in _HOVER_FIELDS if c in df_plot.columns]

    # Build a 2-D numpy array: one column per hover field, plus the
    # cell barcode (DataFrame index) as the last column.
    arrays: list[np.ndarray] = [df_plot[c].values for c in hover_cols]
    arrays.append(df_plot.index.values)
    customdata: np.ndarray = np.column_stack(arrays)

    # Dynamic template built from the actual columns.
    parts: list[str] = []
    for i, col in enumerate(hover_cols):
        label = _HOVER_LABELS.get(col, col.replace("_", " ").title())
        parts.append(f"{label}: %{{customdata[{i}]}}")
    # Last column is always the cell barcode.
    parts.append(f"Barcode: %{{customdata[{len(hover_cols)}]}}")
    hover_template = "<br>".join(parts) + "<extra></extra>"

    # ------------------------------------------------------------------
    # 4. Marker colouring
    # ------------------------------------------------------------------
    marker: dict[str, Any] = _build_marker(
        df_plot=df_plot,
        color_col=color_col,
        color_type=color_type,
        gene_series=gene_series,
    )

    # ------------------------------------------------------------------
    # 5. Trace creation
    # ------------------------------------------------------------------
    if is_3d:
        trace = go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            marker=marker,
            customdata=customdata,
            hovertemplate=hover_template,
        )
    else:
        trace = go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            marker=marker,
            customdata=customdata,
            hovertemplate=hover_template,
        )

    fig = go.Figure(data=[trace])

    # ------------------------------------------------------------------
    # 6. Layout configuration
    # ------------------------------------------------------------------
    _apply_layout(fig, is_3d=is_3d)

    return fig


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _build_marker(
    df_plot: pd.DataFrame,
    color_col: str | None,
    color_type: str,
    gene_series: pd.Series | None,
) -> dict[str, Any]:
    """Construct the *marker* dict for the scatter trace.

    Three colouring modes are supported:

    * **No colouring** (*color_col* is ``None``) — single default colour.
    * **Categorical** — integer-coded categories mapped through a
      discrete colour scale with a labelled colour-bar.
    * **Numeric / gene** — continuous Viridis colour scale with a
      colour-bar title.
    """
    base: dict[str, Any] = {
        "size": _MARKER_SIZE,
        "opacity": _MARKER_OPACITY,
        "line": {"width": _MARKER_LINE_WIDTH},
    }

    # -- Gene expression colouring (BEFORE column-existence check) ----------
    # Gene names don't appear in the metadata DataFrame columns — the
    # expression values come from *gene_series* instead.  We must handle
    # this colouring mode **before** the "No colouring" guard so that a
    # valid gene name (e.g. "PAX6") isn't silently rejected just because
    # it's absent from the metadata parquet.
    if color_type == "gene" and gene_series is not None:
        values = gene_series.loc[df_plot.index].values
        base["color"] = values
        base["colorscale"] = "Viridis"
        base["colorbar"] = {
            "title": color_col.replace("_", " ").title() if color_col else "Expression"
        }
        base["showscale"] = True
        return base

    # -- No colouring -------------------------------------------------------
    if color_col is None or color_col not in df_plot.columns:
        base["color"] = "#636efa"
        return base

    # -- Numeric colouring --------------------------------------------------
    if color_type == "numeric":
        values = df_plot[color_col].values
        base["color"] = values
        base["colorscale"] = "Viridis"
        base["colorbar"] = {"title": color_col.replace("_", " ").title()}
        base["showscale"] = True
        return base

    # -- Categorical colouring (default) ------------------------------------
    categories: list[str] = sorted(df_plot[color_col].astype(str).unique())
    n_cats = len(categories)

    # Map each category to an integer id.
    cat_to_id: dict[str, int] = {cat: i for i, cat in enumerate(categories)}
    color_ids: np.ndarray = np.array(
        [cat_to_id[v] for v in df_plot[color_col].astype(str).values],
        dtype=np.int32,
    )

    # Pick colours, cycling the palette if needed.
    palette = _CATEGORICAL_PALETTE * (1 + n_cats // len(_CATEGORICAL_PALETTE))

    # Build a discrete colour scale — each category gets a contiguous
    # colour-block so the colour-bar renders solid blocks rather than
    # a smooth gradient (misleading for categorical data).
    cmax = float(n_cats)
    colorscale: list[list[float | str]] = []
    for i, cat_color in enumerate(palette[:n_cats]):
        lo = i / cmax
        hi = (i + 1) / cmax
        colorscale.append([lo, cat_color])
        colorscale.append([hi, cat_color])

    base["color"] = color_ids
    base["colorscale"] = colorscale
    base["cmin"] = 0
    base["cmax"] = cmax
    base["colorbar"] = {
        "tickvals": [i + 0.5 for i in range(n_cats)],
        "ticktext": categories,
        "title": color_col.replace("_", " ").title(),
    }
    base["showscale"] = True
    return base


def _apply_layout(fig: go.Figure, is_3d: bool) -> None:
    """Apply the standard UMAP figure layout."""
    axis_style: dict[str, Any] = {
        "showgrid": False,
        "zeroline": False,
        "showticklabels": False,
        "title": {"text": ""},
    }

    if is_3d:
        fig.update_layout(
            scene={
                "xaxis": {**axis_style, "title": {"text": "UMAP 1"}},
                "yaxis": {**axis_style, "title": {"text": "UMAP 2"}},
                "zaxis": {**axis_style, "title": {"text": "UMAP 3"}},
                "aspectmode": "cube",
            },
            dragmode="lasso",
            uirevision="constant",
            margin={"l": 0, "r": 0, "t": 0, "b": 0},
            hovermode="closest",
        )
    else:
        fig.update_layout(
            xaxis={**axis_style, "title": {"text": "UMAP 1"}},
            yaxis={**axis_style, "title": {"text": "UMAP 2"}},
            dragmode="lasso",
            uirevision="constant",
            margin={"l": 0, "r": 0, "t": 0, "b": 0},
            hovermode="closest",
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
        )

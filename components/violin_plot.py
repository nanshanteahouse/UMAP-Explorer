"""
violin_plot.py — Violin plot rendering for scRNA-seq gene expression data.

Provides a single entry point ``make_violin_plot`` that produces a
``plotly.graph_objects.Figure`` with one ``go.Violin`` trace per cell
group (cluster, cell type, etc.), allowing users to compare expression
distributions across populations.

Typical usage (from a Dash callback)::

    from components.violin_plot import make_violin_plot

    fig = make_violin_plot(
        df=metadata_df,
        gene_series=gene_series,
        group_col="cell_type",
        gene_name="PAX6",
    )
    return fig
"""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

from components.umap_figure import _CATEGORICAL_PALETTE, _MARKER_OPACITY

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def make_violin_plot(
    df: pd.DataFrame,
    gene_series: pd.Series | None,
    group_col: str | None,
    gene_name: str | None,
    max_groups: int = 30,
) -> go.Figure:
    """Build a violin plot of gene expression across cell groups.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata DataFrame whose index contains cell barcodes and which
        includes the column specified by *group_col*.
    gene_series : pd.Series or None
        Gene expression values indexed by cell barcode.  When ``None``
        (or when *gene_name* is ``None``), an empty figure with a prompt
        annotation is returned.
    group_col : str or None
        Column in *df* defining the cell groups (clusters, cell types,
        etc.) to split the violin traces by.
    gene_name : str or None
        Display name of the gene, used as the figure title and y-axis
        label.  When ``None`` (or when *gene_series* is ``None``), an
        empty figure with a prompt annotation is returned.
    max_groups : int
        Maximum number of groups to display.  Groups are sorted by cell
        count descending and truncated to this limit.

    Returns
    -------
    go.Figure
        Figure with one ``go.Violin`` trace per group.  Returns an
        empty figure with a centred annotation when the input data is
        missing or insufficient.
    """
    # ------------------------------------------------------------------
    # 1. Edge case — missing gene input
    # ------------------------------------------------------------------
    if gene_series is None or gene_name is None:
        return _empty_figure("Select a gene to view violin plot")

    # ------------------------------------------------------------------
    # 2. Edge case — missing group column
    # ------------------------------------------------------------------
    if group_col is None or group_col not in df.columns:
        return _empty_figure(f"Group column '{group_col}' not found in data")

    # ------------------------------------------------------------------
    # 3. Merge expression values with group annotations
    # ------------------------------------------------------------------
    merge_df = df[[group_col]].join(gene_series.rename("expression"), how="inner")

    # Drop cells with no expression measurement.
    merge_df = merge_df.dropna(subset=["expression"])

    if merge_df.empty:
        return _empty_figure("No expression data available for this selection")

    # ------------------------------------------------------------------
    # 4. Sort groups by cell count, limit to max_groups
    # ------------------------------------------------------------------
    group_counts = merge_df[group_col].value_counts()
    top_groups = group_counts.head(max_groups).index.tolist()

    # Filter to top groups only.
    merge_df = merge_df[merge_df[group_col].isin(top_groups)]

    # ------------------------------------------------------------------
    # 5. Build one violin trace per group
    # ------------------------------------------------------------------
    palette = _CATEGORICAL_PALETTE * (1 + len(top_groups) // len(_CATEGORICAL_PALETTE))

    traces: list[go.Violin] = []
    for i, group_name in enumerate(top_groups):
        group_mask = merge_df[group_col].values == group_name
        expression_values = merge_df.loc[group_mask, "expression"].values

        traces.append(
            go.Violin(
                x=[group_name] * len(expression_values),
                y=expression_values,
                name=str(group_name),
                box_visible=True,
                points="outliers",
                line_color=palette[i],
                fillcolor=palette[i],
                opacity=_MARKER_OPACITY,
                meanline_visible=True,
            )
        )

    # ------------------------------------------------------------------
    # 6. Build figure with layout
    # ------------------------------------------------------------------
    fig = go.Figure(data=traces)

    fig.update_layout(
        title=f"{gene_name} Expression",
        xaxis={"tickangle": -45, "title": {"text": ""}},
        yaxis={"title": {"text": gene_name}},
        template="plotly_white",
        margin={"l": 80, "r": 20, "t": 50, "b": 100},
        hovermode="closest",
        legend={"title": {"text": group_col.replace("_", " ").title()}},
    )

    return fig


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _empty_figure(message: str) -> go.Figure:
    """Return an empty figure with a centred annotation."""
    fig = go.Figure()
    fig.update_layout(
        template="plotly_white",
        xaxis={"visible": False},
        yaxis={"visible": False},
        annotations=[
            {
                "text": message,
                "xref": "paper",
                "yref": "paper",
                "x": 0.5,
                "y": 0.5,
                "showarrow": False,
                "font": {"size": 14, "color": "#888"},
            }
        ],
        margin={"l": 20, "r": 20, "t": 20, "b": 20},
    )
    return fig

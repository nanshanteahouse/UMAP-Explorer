"""
info_panel.py — Bottom information panel displaying statistics about selected/lassoed cells.

Provides two public functions:
    create_info_panel()
        Returns a placeholder ``dash.html.Div`` with the info-panel container.
    build_info_panel_content(selected_df, color_col, total_cells)
        Returns updated children (stats grid + compact charts) reflecting the
        current selection.

These functions are pure layout builders — they contain **no callback logic**.
Callbacks are defined in ``pages/single_dataset.py`` (or equivalent page module)
and call ``build_info_panel_content()`` to produce new children when the user
lasso- or box-selects cells on the UMAP plot.
"""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
from dash import dcc, html


def create_info_panel() -> html.Div:
    """Create an empty / placeholder info panel.

    Returns
    -------
    html.Div
        A ``Div`` with ``id='info-panel-content'`` and CSS class
        ``info-panel``, containing a short instruction message.  The contents
        are replaced by ``build_info_panel_content()`` when the user makes a
        selection on the UMAP plot.
    """
    return html.Div(
        id="info-panel-content",
        className="info-panel",
        children=[
            html.P(
                "Click and drag to select cells on the UMAP",
                style={"color": "var(--text-secondary)", "text-align": "center"},
            )
        ],
    )


def build_info_panel_content(
    selected_df: pd.DataFrame,
    color_col: str,
    total_cells: int,
) -> list:
    """Build updated content for the info panel when cells are selected.

    Parameters
    ----------
    selected_df : pd.DataFrame
        Subset of the metadata DataFrame corresponding to the currently
        selected cells.  May be empty (0 rows) — the function returns a
        placeholder message in that case.
    color_col : str
        Name of the column currently used for colouring the UMAP points
        (e.g. ``'cell_type'``, ``'stage'``, or a gene symbol).
    total_cells : int
        Total number of cells in the full dataset (used to compute the
        selection percentage).

    Returns
    -------
    list
        List of Dash children to place inside the ``info-panel-content`` Div.
        Contains a stats grid, an optional cell-type pie chart, and a
        distribution chart for *color_col*.
    """
    if selected_df is None or len(selected_df) == 0:
        return [
            html.P(
                "Click and drag to select cells on the UMAP",
                style={"color": "var(--text-secondary)", "text-align": "center"},
            )
        ]

    n_selected = len(selected_df)
    pct = n_selected / total_cells * 100

    stat_cards = [
        html.Div(
            className="info-stat",
            children=[
                html.Div("Selected", className="info-stat-label"),
                html.Div(f"{n_selected:,}", className="info-stat-value"),
            ],
        ),
        html.Div(
            className="info-stat",
            children=[
                html.Div("of total", className="info-stat-label"),
                html.Div(f"{pct:.1f}%", className="info-stat-value"),
            ],
        ),
    ]

    if "cell_type" in selected_df.columns:
        n_types = int(selected_df["cell_type"].nunique())
        stat_cards.append(
            html.Div(
                className="info-stat",
                children=[
                    html.Div("Cell types", className="info-stat-label"),
                    html.Div(str(n_types), className="info-stat-value"),
                ],
            ),
        )

    children: list = [
        html.Div(className="info-grid", children=stat_cards),
    ]

    if "cell_type" in selected_df.columns:
        ct_counts = selected_df["cell_type"].value_counts()
        pie = go.Figure(
            data=[
                go.Pie(
                    labels=ct_counts.index.tolist(),
                    values=ct_counts.values,
                    textinfo="label+percent",
                    textposition="inside",
                    insidetextfont={"size": 10},
                    hovertemplate=(
                        "%{label}<br>%{value} cells (%{percent})<extra></extra>"
                    ),
                )
            ],
            layout=go.Layout(
                title={"text": "Cell Type Distribution", "font": {"size": 11}},
                height=140,
                margin={"t": 20, "b": 5, "l": 5, "r": 5},
                template="plotly_white",
                showlegend=True,
                legend={
                    "font": {"size": 10},
                    "orientation": "h",
                    "y": -0.2,
                },
            ),
        )
        children.append(
            dcc.Graph(
                figure=pie,
                style={"height": "140px"},
                config={"displayModeBar": False},
            )
        )

    if color_col in selected_df.columns:
        col_data = selected_df[color_col]
        is_numeric = pd.api.types.is_numeric_dtype(col_data)

        if is_numeric:
            dist = go.Figure(
                data=[
                    go.Histogram(
                        x=col_data.dropna(),
                        marker={
                            "color": "#1a73e8",
                            "line": {"color": "white", "width": 0.5},
                        },
                        hovertemplate=(
                            "Range: %{x}<br>Count: %{y}<extra></extra>"
                        ),
                    )
                ],
                layout=go.Layout(
                    title={
                        "text": f"Distribution of {color_col}",
                        "font": {"size": 11},
                    },
                    height=140,
                    margin={"t": 20, "b": 20, "l": 35, "r": 5},
                    template="plotly_white",
                    xaxis={
                        "title": {"text": color_col, "font": {"size": 9}}
                    },
                    yaxis={"title": {"text": "Count", "font": {"size": 9}}},
                ),
            )
        else:
            vc = col_data.value_counts().head(20)
            dist = go.Figure(
                data=[
                    go.Bar(
                        x=vc.index.tolist(),
                        y=vc.values,
                        marker={"color": "#1a73e8"},
                        hovertemplate="%{x}<br>%{y} cells<extra></extra>",
                    )
                ],
                layout=go.Layout(
                    title={
                        "text": f"Distribution of {color_col}",
                        "font": {"size": 11},
                    },
                    height=140,
                    margin={"t": 20, "b": 40, "l": 35, "r": 5},
                    template="plotly_white",
                    xaxis={
                        "tickangle": -45,
                        "title": {
                            "text": color_col,
                            "font": {"size": 9},
                        },
                    },
                    yaxis={
                        "title": {"text": "Count", "font": {"size": 9}}
                    },
                ),
            )

        children.append(
            dcc.Graph(
                figure=dist,
                style={"height": "140px"},
                config={"displayModeBar": False},
            )
        )

    return children

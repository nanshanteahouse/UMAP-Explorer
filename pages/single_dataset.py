"""Single-dataset UMAP exploration page with all Dash callbacks."""
from __future__ import annotations

import dash
import pandas as pd
from dash import (
    Input,
    Output,
    State,
    callback,
    dcc,
    html,
    no_update,
)
from typing import Any

from components.controls import create_controls_panel
from components.data_loader import DataLoader
from components.gene_selector import get_gene_options
from components.info_panel import build_info_panel_content, create_info_panel
from components.umap_figure import make_umap_figure
from components.violin_plot import make_violin_plot

loader = DataLoader()
_metadata_cache: dict[str, pd.DataFrame] = {}

dash.register_page(
    __name__,
    path="/",
    name="Single Dataset",
    title="Single Dataset Explorer",
)


def _get_metadata(dataset_id: str) -> pd.DataFrame:
    if dataset_id not in _metadata_cache:
        _metadata_cache[dataset_id] = loader.load_metadata(dataset_id)
    return _metadata_cache[dataset_id]


def _error_figure(gene: str | None, group_col: str | None) -> Any:
    """Return a figure showing a gene-not-found message."""
    import plotly.graph_objects as go
    message = f"Gene '{gene}' expression not available in this dataset"
    if group_col:
        message += f"\n(Try selecting a different gene or check Group by: '{group_col}')"
    fig = go.Figure()
    fig.update_layout(
        template="plotly_white",
        xaxis={"visible": False},
        yaxis={"visible": False},
        annotations=[{
            "text": message,
            "xref": "paper", "yref": "paper",
            "x": 0.5, "y": 0.5,
            "showarrow": False,
            "font": {"size": 13, "color": "#c5221f"},
        }],
        margin={"l": 20, "r": 20, "t": 20, "b": 20},
    )
    return fig


def get_single_dataset_layout() -> html.Div:
    registry = loader.load_registry()

    return html.Div(
        [
            html.Div(
                create_controls_panel(registry),
                className="sidebar",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            dcc.Graph(
                                id="umap-graph",
                                style={"height": "100%"},
                                config={
                                    "displayModeBar": True,
                                    "scrollZoom": True,
                                },
                            ),
                            dcc.Graph(
                                id="violin-graph",
                                style={"display": "none"},
                                config={
                                    "displayModeBar": True,
                                    "scrollZoom": False,
                                },
                            ),
                        ],
                        className="umap-container",
                    ),
                    create_info_panel(),
                ],
                className="content-area",
            ),
        ],
        className="main-container",
    )


layout = get_single_dataset_layout


@callback(
    Output("color-selector", "options"),
    Output("color-selector", "value"),
    Output("gene-selector", "options"),
    Output("gene-selector", "value"),
    Output("dataset-store", "data"),
    Output("sample-info", "children"),
    Output("violin-groupby-selector", "options"),
    Output("violin-groupby-selector", "value"),
    Input("dataset-selector", "value"),
)
def on_dataset_select(dataset_id: str) -> tuple:
    if not dataset_id:
        return [], None, [], None, None, "No dataset selected", [], None

    df = _get_metadata(dataset_id)
    n_cells = len(df)

    columns = loader.get_available_columns(dataset_id)
    cat_opts: list[dict] = []
    num_opts: list[dict] = []

    for col in columns:
        col_type = loader.get_column_type(dataset_id, col)
        base_label = col.replace("_", " ").title()
        entry = {"label": f"[C] {base_label}", "value": col}
        if col_type in ("categorical", "bool"):
            cat_opts.append(entry)
        else:
            entry["label"] = f"[N] {base_label}"
            num_opts.append(entry)

    cat_opts.sort(key=lambda x: x["label"])
    num_opts.sort(key=lambda x: x["label"])
    color_options = cat_opts + num_opts

    gene_index = loader.load_gene_index()
    available_genes = set(loader.get_available_genes(dataset_id))
    gene_options = get_gene_options(dataset_id, gene_index, available_genes=available_genes)

    sample_size = min(n_cells, 50000)
    sample_info = f"{n_cells:,} cells available, rendering {sample_size:,}"

    cat_cols = [c for c in columns if loader.get_column_type(dataset_id, c) in ("categorical", "bool")]
    groupby_options = [
        {"label": c.replace("_", " ").title(), "value": c} for c in sorted(cat_cols)
    ]
    default_groupby = None
    for candidate in ("cell_type", "cell_type_sub", "predicted_celltype", "predicted_cell_type"):
        if candidate in cat_cols:
            default_groupby = candidate
            break
    if default_groupby is None and cat_cols:
        default_groupby = cat_cols[0]

    return (
        color_options, no_update, gene_options, None,
        dataset_id, sample_info,
        groupby_options, default_groupby,
    )


@callback(
    Output("umap-graph", "figure"),
    Input("dataset-store", "data"),
    Input("color-selector", "value"),
    Input("gene-selector", "value"),
    Input("dim-toggle", "value"),
    Input("size-slider", "value"),
    Input("opacity-slider", "value"),
    prevent_initial_call=True,
)
def _update_umap_figure(
    dataset_id: str | None,
    color_col: str | None,
    gene: str | None,
    dim: str,
    size: float,
    opacity: float,
) -> Any:
    if not dataset_id:
        return no_update

    df = _get_metadata(dataset_id)

    if gene:
        try:
            gene_series = loader.load_gene_expression(dataset_id, gene)
        except (KeyError, Exception):
            return no_update

        fig = make_umap_figure(
            df,
            color_col=gene,
            dim=dim.lower(),
            color_type="gene",
            gene_series=gene_series,
        )
        fig.update_traces(marker={"size": size, "opacity": opacity})
        return fig

    if color_col and color_col in df.columns:
        col_type = loader.get_column_type(dataset_id, color_col)
        fig = make_umap_figure(
            df,
            color_col=color_col,
            dim=dim.lower(),
            color_type=col_type,
        )
    else:
        fig = make_umap_figure(df, color_col=None, dim=dim.lower())

    fig.update_traces(marker={"size": size, "opacity": opacity})
    return fig


@callback(
    Output("info-panel-content", "children"),
    Input("umap-graph", "selectedData"),
    State("dataset-store", "data"),
    State("color-selector", "value"),
)
def _on_lasso_select(
    selected_data: dict | None,
    dataset_id: str | None,
    color_col: str | None,
) -> list:
    if not dataset_id:
        return build_info_panel_content(None, "", 0)

    df = _get_metadata(dataset_id)
    total_n = len(df)

    if not selected_data or not selected_data.get("points"):
        return build_info_panel_content(None, color_col or "", total_n)

    barcodes = [p["customdata"][-1] for p in selected_data["points"]]
    selected_df = df.loc[df.index.isin(barcodes)]

    return build_info_panel_content(selected_df, color_col or "", total_n)


@callback(
    Output("violin-graph", "figure"),
    Input("dataset-store", "data"),
    Input("gene-selector", "value"),
    Input("violin-groupby-selector", "value"),
    prevent_initial_call=True,
)
def _update_violin_figure(
    dataset_id: str | None,
    gene: str | None,
    group_col: str | None,
) -> Any:
    if not dataset_id:
        return no_update

    df = _get_metadata(dataset_id)

    if not gene:
        return make_violin_plot(df, None, group_col, None)

    try:
        gene_series = loader.load_gene_expression(dataset_id, gene)
    except (KeyError, Exception):
        return _error_figure(gene, group_col)

    return make_violin_plot(df, gene_series, group_col, gene)


@callback(
    Output("umap-graph", "style"),
    Output("violin-graph", "style"),
    Input("plot-type-toggle", "value"),
)
def _toggle_plot_view(view_mode: str) -> tuple:
    if view_mode == "violin":
        return {"display": "none"}, {"display": "block", "height": "100%"}
    return {"display": "block", "height": "100%"}, {"display": "none"}

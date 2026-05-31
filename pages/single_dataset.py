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

from components.controls import create_controls_panel
from components.data_loader import DataLoader
from components.gene_selector import get_gene_options
from components.info_panel import build_info_panel_content, create_info_panel
from components.umap_figure import make_umap_figure

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
    Input("dataset-selector", "value"),
)
def on_dataset_select(dataset_id: str) -> tuple:
    if not dataset_id:
        return [], None, [], None, None, "No dataset selected"

    df = _get_metadata(dataset_id)
    n_cells = len(df)

    columns = loader.get_available_columns(dataset_id)
    cat_options: list[dict] = []
    num_options: list[dict] = []

    for col in columns:
        col_type = loader.get_column_type(dataset_id, col)
        label = col.replace("_", " ").title()
        entry: dict[str, str] = {"label": label, "value": col}
        if col_type in ("categorical", "bool"):
            cat_options.append(entry)
        else:
            num_options.append(entry)

    color_options: list[dict] = []
    if cat_options:
        color_options.append({"label": "Categorical", "options": cat_options})
    if num_options:
        color_options.append({"label": "Numeric", "options": num_options})

    gene_index = loader.load_gene_index()
    gene_options = get_gene_options(dataset_id, gene_index)

    sample_size = min(n_cells, 50000)
    sample_info = f"{n_cells:,} cells available, rendering {sample_size:,}"

    return color_options, None, gene_options, None, dataset_id, sample_info


@callback(
    Output("umap-graph", "figure"),
    Output("color-selector", "value", allow_duplicate=True),
    Input("dataset-store", "data"),
    Input("color-selector", "value"),
    Input("gene-selector", "value"),
    Input("dim-toggle", "value"),
    Input("size-slider", "value"),
    Input("opacity-slider", "value"),
)
def _update_umap_figure(
    dataset_id: str | None,
    color_col: str | None,
    gene: str | None,
    dim: str,
    size: float,
    opacity: float,
) -> tuple:
    if not dataset_id:
        return no_update, no_update

    df = _get_metadata(dataset_id)

    if gene:
        try:
            gene_series = loader.load_gene_expression(dataset_id, gene)
        except (KeyError, Exception):
            return no_update, no_update

        fig = make_umap_figure(
            df,
            color_col=gene,
            dim=dim.lower(),
            color_type="gene",
            gene_series=gene_series,
        )
        fig.update_traces(marker={"size": size, "opacity": opacity})

        # Only clear color-selector when transitioning into gene mode;
        # skip if already None to avoid cascading re-triggers.
        if color_col is not None:
            return fig, None
        return fig, no_update

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
    return fig, no_update


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

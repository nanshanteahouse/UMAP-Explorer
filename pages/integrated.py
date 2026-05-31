"""
integrated.py — Integrated cross-dataset UMAP view page.

Displays all datasets in a shared Harmony-corrected UMAP embedding.
Requires outputs from build_integrated.py to exist in data/integrated/.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import dash
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, callback, dcc, html, no_update

from components.umap_figure import make_umap_figure

# ---------------------------------------------------------------------------
# Path constants
# ---------------------------------------------------------------------------
INTEGRATED_DIR = str(Path(__file__).resolve().parent.parent / "data" / "integrated")
UMAP_PARQUET = os.path.join(INTEGRATED_DIR, "integrated_umap.parquet")
INFO_PATH = os.path.join(INTEGRATED_DIR, "integrated_info.json")

dash.register_page(__name__, path="/integrated", name="Integrated", title="Integrated View")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_COLORBY_BASE = [
    {"label": "Dataset", "value": "dataset_id"},
    {"label": "Unified Cell Type", "value": "unified_cell_type"},
]

# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def _check_data_exists() -> bool:
    return os.path.isfile(UMAP_PARQUET) and os.path.isfile(INFO_PATH)


def _load_info() -> dict | None:
    if not _check_data_exists():
        return None
    with open(INFO_PATH) as f:
        return json.load(f)


def _load_umap() -> pd.DataFrame | None:
    if not _check_data_exists():
        return None
    df = pd.read_parquet(UMAP_PARQUET)
    for col in [c for c in df.columns if c.startswith("UMAP_")]:
        df[col] = df[col].astype("float32")
    return df


# ---------------------------------------------------------------------------
# Placeholder layout (shown when data is missing)
# ---------------------------------------------------------------------------

_PLACEHOLDER = html.Div(
    [
        html.H3("Integrated View Not Available"),
        html.P("Please run the following command to generate the integrated dataset:"),
        html.Pre(
            "python vis_website/preprocessing/build_integrated.py \\\n"
            "    --data-dir /path/to/neurobiology \\\n"
            "    --output-dir vis_website/data/integrated",
            style={
                "background": "#f0f0f0",
                "padding": "16px",
                "borderRadius": "8px",
                "margin": "16px 0",
                "fontSize": "13px",
                "maxWidth": "700px",
            },
        ),
    ],
    style={
        "display": "flex",
        "flexDirection": "column",
        "alignItems": "center",
        "justifyContent": "center",
        "height": "calc(100vh - 56px)",
        "color": "#5f6368",
        "textAlign": "center",
        "padding": "24px",
    },
)


# ---------------------------------------------------------------------------
# Page builder
# ---------------------------------------------------------------------------

def _build_page(umap_df: pd.DataFrame) -> html.Div:
    """Construct the integrated page layout once data is confirmed to exist."""
    n_cells = len(umap_df)
    n_datasets = umap_df["dataset_id"].nunique() if "dataset_id" in umap_df.columns else 0

    available = [c for c in umap_df.columns if not c.startswith("UMAP_") and c != "dataset_id"]
    extra_cats = sorted(c for c in available if umap_df[c].dtype.name == "category" or umap_df[c].nunique() < 50)
    colorby_options = list(_COLORBY_BASE)
    for col in extra_cats:
        colorby_options.append({"label": col.replace("_", " ").title(), "value": col})

    dataset_ids = sorted(umap_df["dataset_id"].unique()) if "dataset_id" in umap_df.columns else []
    cell_types = sorted(umap_df["unified_cell_type"].dropna().unique()) if "unified_cell_type" in umap_df.columns else []

    size_info = f"{n_cells:,} cells across {n_datasets} datasets"
    if n_cells > 50000:
        size_info += " (rendering 50,000)"

    controls = html.Div(
        className="sidebar",
        style={"width": "300px", "minWidth": "300px"},
        children=[
            html.Div(
                className="sidebar-section",
                children=[
                    html.Div("Color By", className="sidebar-section-title"),
                    dcc.Dropdown(
                        id="integrated-color-selector",
                        options=colorby_options,
                        value="dataset_id",
                        clearable=False,
                        searchable=True,
                    ),
                ],
            ),
            html.Div(
                className="sidebar-section",
                children=[
                    html.Div("Datasets", className="sidebar-section-title"),
                    dcc.Checklist(
                        id="integrated-dataset-filter",
                        options=[{"label": d, "value": d} for d in dataset_ids],
                        value=list(dataset_ids),
                        labelStyle={"display": "block", "fontSize": "13px", "marginBottom": "4px", "cursor": "pointer"},
                        inputStyle={"marginRight": "6px", "cursor": "pointer"},
                    ),
                    html.Button(
                        "Select All", id="integrated-dataset-select-all",
                        style={"fontSize": "12px", "padding": "2px 8px", "marginRight": "4px", "cursor": "pointer"},
                    ),
                    html.Button(
                        "Clear", id="integrated-dataset-clear",
                        style={"fontSize": "12px", "padding": "2px 8px", "cursor": "pointer"},
                    ),
                ],
            ),
            html.Div(
                className="sidebar-section",
                children=[
                    html.Div("Cell Types", className="sidebar-section-title"),
                    dcc.Dropdown(
                        id="integrated-celltype-filter",
                        options=[{"label": ct, "value": ct} for ct in cell_types],
                        value=None,
                        clearable=True,
                        searchable=True,
                        placeholder="All cell types…",
                        multi=True,
                    ),
                ],
            ),
            html.Div(
                className="sidebar-section",
                children=[
                    html.Div("Overview", className="sidebar-section-title"),
                    html.Div(
                        id="integrated-overview",
                        children=size_info,
                        style={"fontSize": "13px", "color": "var(--text-secondary)"},
                    ),
                ],
            ),
        ],
    )

    content = html.Div(
        className="content-area",
        style={"flexDirection": "row"},
        children=[
            html.Div(
                className="umap-container",
                style={"flex": "3", "position": "relative"},
                children=[
                    dcc.Graph(
                        id="integrated-umap-graph",
                        style={"height": "100%"},
                        config={"displayModeBar": True, "scrollZoom": True},
                    ),
                ],
            ),
            html.Div(
                style={
                    "width": "300px", "minWidth": "300px", "padding": "12px",
                    "overflowY": "auto", "display": "flex", "flexDirection": "column", "gap": "12px",
                },
                children=[
                    html.Div(
                        className="sidebar-section",
                        children=[dcc.Graph(id="integrated-dataset-pie", style={"height": "250px"}, config={"displayModeBar": False})],
                    ),
                    html.Div(
                        className="sidebar-section",
                        children=[dcc.Graph(id="integrated-ct-bar", style={"height": "300px"}, config={"displayModeBar": False})],
                    ),
                ],
            ),
        ],
    )

    return html.Div([controls, content], className="main-container")


# ---------------------------------------------------------------------------
# Module-level layout assignment
# ---------------------------------------------------------------------------
if _check_data_exists():
    _df = _load_umap()
    layout = _build_page(_df) if _df is not None and len(_df) > 0 else _PLACEHOLDER
else:
    layout = _PLACEHOLDER


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

def _filter_df(df, datasets, celltypes):
    if datasets is not None:
        df = df[df["dataset_id"].isin(datasets)]
    if celltypes and "unified_cell_type" in df.columns:
        df = df[df["unified_cell_type"].isin(celltypes)]
    return df


@callback(
    Output("integrated-umap-graph", "figure"),
    Input("integrated-color-selector", "value"),
    Input("integrated-dataset-filter", "value"),
    Input("integrated-celltype-filter", "value"),
    prevent_initial_call=False,
)
def _update_umap(color_col, datasets, celltypes):
    df = _load_umap()
    if df is None or len(df) == 0:
        return no_update
    df = _filter_df(df, datasets, celltypes)
    if len(df) == 0:
        return go.Figure({"data": [], "layout": {"annotations": [{"text": "No cells match filters", "showarrow": False}]}})
    is_num = color_col and color_col in df.columns and pd.api.types.is_numeric_dtype(df[color_col].dtype) and not pd.api.types.is_bool_dtype(df[color_col].dtype) and df[color_col].nunique() > 20
    return make_umap_figure(df, color_col=color_col, color_type="numeric" if is_num else "categorical", sample_size=50000)


@callback(
    Output("integrated-dataset-pie", "figure"),
    Input("integrated-dataset-filter", "value"),
    Input("integrated-celltype-filter", "value"),
    prevent_initial_call=False,
)
def _update_pie(datasets, celltypes):
    df = _load_umap()
    if df is None or "dataset_id" not in df.columns:
        return go.Figure()
    df = _filter_df(df, None, celltypes)
    counts = df["dataset_id"].value_counts()
    return go.Figure(
        data=[go.Pie(labels=counts.index.tolist(), values=counts.values, textinfo="label+percent", textposition="inside", insidetextfont={"size": 9}, hovertemplate="%{label}<br>%{value} cells (%{percent})<extra></extra>")],
        layout={"title": {"text": "Dataset Distribution", "font": {"size": 12}}, "height": 250, "margin": {"t": 30, "b": 10, "l": 10, "r": 10}, "template": "plotly_white", "showlegend": False},
    )


@callback(
    Output("integrated-ct-bar", "figure"),
    Input("integrated-dataset-filter", "value"),
    Input("integrated-celltype-filter", "value"),
    prevent_initial_call=False,
)
def _update_ct_bar(datasets, celltypes):
    df = _load_umap()
    if df is None or "unified_cell_type" not in df.columns or "dataset_id" not in df.columns:
        return go.Figure()
    df = _filter_df(df, datasets, celltypes)
    ct = df.groupby(["dataset_id", "unified_cell_type"], observed=True).size().reset_index(name="count")
    fig = go.Figure()
    for ct_name in sorted(ct["unified_cell_type"].unique()):
        sub = ct[ct["unified_cell_type"] == ct_name]
        fig.add_trace(go.Bar(name=ct_name, x=sub["dataset_id"].tolist(), y=sub["count"].values, hovertemplate="%{x}<br>%{y} cells<extra></extra>"))
    fig.update_layout(barmode="stack", title={"text": "Cell Types by Dataset", "font": {"size": 12}}, height=300, margin={"t": 30, "b": 60, "l": 40, "r": 10}, template="plotly_white", legend={"font": {"size": 9}, "orientation": "h", "y": -0.25}, xaxis={"tickangle": -45, "title": {"text": ""}}, yaxis={"title": {"text": "Cells", "font": {"size": 10}}})
    return fig


@callback(
    Output("integrated-dataset-filter", "value"),
    Input("integrated-dataset-select-all", "n_clicks"),
    Input("integrated-dataset-clear", "n_clicks"),
    State("integrated-dataset-filter", "options"),
    prevent_initial_call=True,
)
def _handle_dataset_buttons(select_clicks, clear_clicks, options):
    ctx = dash.callback_context
    if not ctx.triggered:
        return no_update
    bid = ctx.triggered[0]["prop_id"].split(".")[0]
    return [opt["value"] for opt in options] if bid == "integrated-dataset-select-all" else []

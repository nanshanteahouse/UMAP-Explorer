"""comparison.py — Cross-dataset comparison page with side-by-side UMAPs."""
from __future__ import annotations

import dash
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, callback, dcc, html, no_update
from dash.exceptions import PreventUpdate
from components.data_loader import DataLoader
from components.gene_selector import get_gene_options
from components.umap_figure import make_umap_figure

loader = DataLoader()
_cache: dict[str, pd.DataFrame] = {}

dash.register_page(__name__, path="/comparison", name="Compare", title="Cross-Dataset Comparison")


def _meta(dsid: str) -> pd.DataFrame:
    if dsid not in _cache:
        _cache[dsid] = loader.load_metadata(dsid)
    return _cache[dsid]


def _color_opts(dsid: str) -> list[dict]:
    cols = loader.get_available_columns(dsid)
    opts: list[dict] = []
    for c in cols:
        t = loader.get_column_type(dsid, c)
        base_label = c.replace("_", " ").title()
        prefix = "[C] " if t in ("categorical", "bool") else "[N] "
        opts.append({"label": f"{prefix}{base_label}", "value": c})
    return opts


def _shared_color_opts(lid: str, rid: str) -> list[dict]:
    common = sorted(set(loader.get_available_columns(lid)) & set(loader.get_available_columns(rid)))
    opts: list[dict] = []
    for c in common:
        t = loader.get_column_type(lid, c)
        base_label = c.replace("_", " ").title()
        prefix = "[C] " if t in ("categorical", "bool") else "[N] "
        opts.append({"label": f"{prefix}{base_label}", "value": c})
    return opts


def _shared_gene_opts(lid: str, rid: str) -> list[dict]:
    gi = loader.load_gene_index().get("genes", {})
    lset = {g for g, v in gi.items() if lid in v.get("datasets", [])}
    rset = {g for g, v in gi.items() if rid in v.get("datasets", [])}
    return [{"label": g, "value": g} for g in sorted(lset & rset)]


def _fig(dsid: str | None, cc: str | None, gene: str | None, side: str = "left") -> go.Figure | type[no_update]:
    if not dsid:
        return no_update
    df = _meta(dsid)
    if gene:
        try:
            gs = loader.load_gene_expression(dsid, gene)
        except (KeyError, Exception):
            return no_update
        fig = make_umap_figure(df, color_col=gene, color_type="gene", gene_series=gs)
    elif cc and cc in df.columns:
        fig = make_umap_figure(df, color_col=cc, color_type=loader.get_column_type(dsid, cc))
    else:
        fig = make_umap_figure(df, color_col=None)
    fig.update_layout(uirevision=f"{side}-{dsid}")
    return fig


def _apply_zoom(fig: dict | go.Figure, rd: dict) -> go.Figure:
    fig = go.Figure(fig)
    xr: list = [None, None]
    yr: list = [None, None]
    for k, v in rd.items():
        if k == "xaxis.range[0]":
            xr[0] = v
        elif k == "xaxis.range[1]":
            xr[1] = v
        elif k == "yaxis.range[0]":
            yr[0] = v
        elif k == "yaxis.range[1]":
            yr[1] = v
        elif k == "xaxis.range" and isinstance(v, list):
            xr = v
        elif k == "yaxis.range" and isinstance(v, list):
            yr = v
        elif k == "xaxis.autorange":
            fig.update_xaxes(autorange=v)
        elif k == "yaxis.autorange":
            fig.update_yaxes(autorange=v)
    if xr[0] is not None and xr[1] is not None:
        fig.update_xaxes(range=xr)
    if yr[0] is not None and yr[1] is not None:
        fig.update_yaxes(range=yr)
    return fig


def _panel(side: str, opts: list[dict]) -> html.Div:
    return html.Div([
        html.Div(dcc.Dropdown(id=f"{side}-dataset", options=opts, placeholder=f"Select {side} dataset..."),
                 className="comparison-panel-header"),
        html.Div([
            dcc.Dropdown(id=f"{side}-color", placeholder="Color by...", style={"marginBottom": "4px"}),
            dcc.Dropdown(id=f"{side}-gene", placeholder="Gene expression...", searchable=True,
                         style={"marginBottom": "4px"}),
            dcc.Graph(id=f"{side}-umap", style={"height": "100%", "minHeight": "400px"},
                      config={"displayModeBar": True, "scrollZoom": True}),
        ], className="comparison-panel-body", style={"padding": "8px", "display": "flex", "flexDirection": "column"}),
    ], className="comparison-panel")


S = {"fontSize": "13px", "fontWeight": 600, "marginRight": "8px"}


def layout(**kwargs) -> html.Div:
    registry = loader.load_registry()
    dopts = [{"label": f"{ds['dataset_id']} ({ds.get('species', '?')}, {ds['n_cells']:,} cells)", "value": ds["dataset_id"]} for ds in registry]
    return html.Div([
        html.Div([
            html.Span("Synchronize:", style={"fontSize": "13px", "fontWeight": 600}),
            dcc.Checklist(id="sync-mode", options=[{"label": " Zoom", "value": "zoom"},
                                                    {"label": " Color", "value": "color"},
                                                    {"label": " Gene", "value": "gene"}],
                          inline=True, style={"marginRight": "16px"}),
            html.Span("Color by:", style=S),
            html.Div(dcc.Dropdown(id="shared-color", placeholder="Shared color..."),
                     style={"width": "200px", "marginRight": "16px"}),
            html.Span("Gene:", style=S),
            html.Div(dcc.Dropdown(id="shared-gene", placeholder="Shared gene...", searchable=True),
                     style={"width": "250px"}),
        ], style={"padding": "12px 16px", "background": "white", "borderBottom": "1px solid var(--border)",
                  "display": "flex", "alignItems": "center", "gap": "8px", "flexWrap": "wrap"}),
        html.Div([_panel("left", dopts), _panel("right", dopts)],
                 className="comparison-container", style={"padding": "8px"}),
        dcc.Store(id="left-dataset-store"),
        dcc.Store(id="right-dataset-store"),
    ], className="main-container", style={"flexDirection": "column"})


@callback(
    Output("left-color", "options"), Output("left-color", "value"),
    Output("left-gene", "options"), Output("left-gene", "value"),
    Output("left-dataset-store", "data"),
    Input("left-dataset", "value"),
)
def _on_left_dataset(dsid: str | None) -> tuple:
    if not dsid:
        return [], None, [], None, None
    return _color_opts(dsid), None, get_gene_options(dsid, loader.load_gene_index()), None, dsid


@callback(
    Output("right-color", "options"), Output("right-color", "value"),
    Output("right-gene", "options"), Output("right-gene", "value"),
    Output("right-dataset-store", "data"),
    Input("right-dataset", "value"),
)
def _on_right_dataset(dsid: str | None) -> tuple:
    if not dsid:
        return [], None, [], None, None
    return _color_opts(dsid), None, get_gene_options(dsid, loader.load_gene_index()), None, dsid


@callback(
    Output("shared-color", "options"), Output("shared-gene", "options"),
    Input("left-dataset-store", "data"), Input("right-dataset-store", "data"),
)
def _on_shared_opts(lid: str | None, rid: str | None) -> tuple:
    if not lid or not rid:
        return [], []
    return _shared_color_opts(lid, rid), _shared_gene_opts(lid, rid)


@callback(
    Output("left-umap", "figure"),
    Input("left-dataset-store", "data"), Input("left-color", "value"), Input("left-gene", "value"),
    Input("shared-color", "value"), Input("shared-gene", "value"), Input("sync-mode", "value"),
)
def _left_umap(dsid: str | None, lc: str | None, lg: str | None,
               sc: str | None, sg: str | None, sync: list[str] | None) -> go.Figure | type[no_update]:
    sync = sync or []
    return _fig(dsid, sc if "color" in sync and sc else lc, sg if "gene" in sync and sg else lg, "left")


@callback(
    Output("right-umap", "figure"),
    Input("right-dataset-store", "data"), Input("right-color", "value"), Input("right-gene", "value"),
    Input("shared-color", "value"), Input("shared-gene", "value"), Input("sync-mode", "value"),
)
def _right_umap(dsid: str | None, lc: str | None, lg: str | None,
                sc: str | None, sg: str | None, sync: list[str] | None) -> go.Figure | type[no_update]:
    sync = sync or []
    return _fig(dsid, sc if "color" in sync and sc else lc, sg if "gene" in sync and sg else lg, "right")


@callback(
    Output("right-umap", "figure", allow_duplicate=True),
    Input("left-umap", "relayoutData"),
    State("right-umap", "figure"), State("sync-mode", "value"),
    prevent_initial_call=True,
)
def _zoom_l2r(rd: dict | None, fig: dict | None, sync: list[str] | None) -> go.Figure | type[no_update]:
    if not sync or "zoom" not in sync or not rd or fig is None:
        raise PreventUpdate
    return _apply_zoom(fig, rd)


@callback(
    Output("left-umap", "figure", allow_duplicate=True),
    Input("right-umap", "relayoutData"),
    State("left-umap", "figure"), State("sync-mode", "value"),
    prevent_initial_call=True,
)
def _zoom_r2l(rd: dict | None, fig: dict | None, sync: list[str] | None) -> go.Figure | type[no_update]:
    if not sync or "zoom" not in sync or not rd or fig is None:
        raise PreventUpdate
    return _apply_zoom(fig, rd)

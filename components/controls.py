"""Sidebar controls panel for the UMAP explorer."""

from __future__ import annotations

from dash import dcc, html


def create_controls_panel(registry: list[dict]) -> html.Div:
    """Build the sidebar controls panel.

    Parameters
    ----------
    registry : list[dict]
        Dataset registry entries; each must contain at least
        ``dataset_id``, ``n_cells``, and ``species``.

    Returns
    -------
    html.Div
        Sidebar controls panel ready to be placed inside
        a ``className='sidebar'`` container.
    """
    dataset_options = [
        {
            "label": f'{d["dataset_id"]} ({d["n_cells"]} cells, {d["species"]})',
            "value": d["dataset_id"],
        }
        for d in registry
    ]
    default_dataset = dataset_options[0]["value"] if dataset_options else None

    dataset_section = html.Div(
        className="sidebar-section",
        children=[
            html.Div("Dataset", className="sidebar-section-title"),
            dcc.Dropdown(
                id="dataset-selector",
                options=dataset_options,
                value=default_dataset,
                clearable=False,
                searchable=True,
                placeholder="Select dataset…",
            ),
            dcc.Store(id="dataset-store"),
        ],
    )

    coloring_section = html.Div(
        className="sidebar-section",
        children=[
            html.Div("Coloring", className="sidebar-section-title"),
            html.Label("Color by", style={"fontSize": "13px", "marginBottom": "4px", "display": "block"}),
            dcc.Dropdown(
                id="color-selector",
                options=[],
                value=None,
                clearable=True,
                searchable=True,
                placeholder="Select column…",
            ),
            html.Div(style={"height": "8px"}),
            html.Label("Gene expression", style={"fontSize": "13px", "marginBottom": "4px", "display": "block"}),
            dcc.Dropdown(
                id="gene-selector",
                options=[],
                value=None,
                clearable=False,
                searchable=True,
                placeholder="Search gene…",
            ),
            html.Div(style={"height": "8px"}),
            html.Label("Group by (Violin)", style={"fontSize": "13px", "marginBottom": "4px", "display": "block"}),
            dcc.Dropdown(
                id="violin-groupby-selector",
                options=[],
                value=None,
                clearable=False,
                searchable=True,
                placeholder="cell_type (auto)...",
            ),
        ],
    )

    display_section = html.Div(
        className="sidebar-section",
        children=[
            html.Div("Display", className="sidebar-section-title"),
            html.Label("Plot mode", style={"fontSize": "13px", "marginBottom": "4px", "display": "block"}),
            dcc.RadioItems(
                id="dim-toggle",
                options=[
                    {"label": " 2D", "value": "2D"},
                    {"label": " 3D", "value": "3D"},
                ],
                value="2D",
                labelStyle={
                    "display": "inline-block",
                    "marginRight": "16px",
                    "fontSize": "13px",
                    "cursor": "pointer",
                },
                inputStyle={"marginRight": "4px"},
            ),
            html.Div(style={"height": "8px"}),
            html.Label("View", style={"fontSize": "13px", "marginBottom": "4px", "display": "block"}),
            dcc.RadioItems(
                id="plot-type-toggle",
                options=[
                    {"label": " UMAP", "value": "umap"},
                    {"label": " Violin", "value": "violin"},
                ],
                value="umap",
                labelStyle={
                    "display": "inline-block",
                    "marginRight": "16px",
                    "fontSize": "13px",
                    "cursor": "pointer",
                },
                inputStyle={"marginRight": "4px"},
                inline=True,
            ),
            html.Div(style={"height": "8px"}),
            html.Label("Point size", style={"fontSize": "13px", "marginBottom": "4px", "display": "block"}),
            dcc.Slider(
                id="size-slider",
                min=1,
                max=10,
                step=0.5,
                value=3,
                marks={1: "1", 3: "3", 5: "5", 7: "7", 10: "10"},
                tooltip={"placement": "bottom", "always_visible": False},
            ),
            html.Div(style={"height": "12px"}),
            html.Label("Opacity", style={"fontSize": "13px", "marginBottom": "4px", "display": "block"}),
            dcc.Slider(
                id="opacity-slider",
                min=0.1,
                max=1.0,
                step=0.1,
                value=0.7,
                marks={0.1: "0.1", 0.3: "0.3", 0.5: "0.5", 0.7: "0.7", 1.0: "1.0"},
                tooltip={"placement": "bottom", "always_visible": False},
            ),
        ],
    )

    selection_section = html.Div(
        className="sidebar-section",
        children=[
            html.Div("Selection Info", className="sidebar-section-title"),
            html.Div(
                id="sample-info",
                children="Rendering 0 of 0 cells",
                style={"fontSize": "13px", "color": "var(--text-secondary)"},
            ),
        ],
    )

    return html.Div(
        children=[dataset_section, coloring_section, display_section, selection_section],
    )

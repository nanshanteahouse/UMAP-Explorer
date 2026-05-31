"""
app.py — Interactive UMAP Explorer Dash application entry point.

This module initialises the Dash app, defines the global layout
(header + navigation), registers all page modules, and exports the
``server`` WSGI object for gunicorn.

Usage
-----
Development:
    python vis_website/app.py

Production:
    gunicorn vis_website.app:server --workers 4 --timeout 120
"""

from __future__ import annotations

import dash
from dash import dcc, html

# ---------------------------------------------------------------------------
# Dash app initialisation
# ---------------------------------------------------------------------------

app = dash.Dash(
    __name__,
    use_pages=True,
    pages_folder="pages",
    suppress_callback_exceptions=True,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"},
    ],
    title="Interactive UMAP Explorer",
    update_title="Loading...",
)

# Expose the underlying Flask server for WSGI runners (gunicorn).
server = app.server

# ---------------------------------------------------------------------------
# Global layout — rendered on every page
# ---------------------------------------------------------------------------
# Each individual page defines its own content inside a .main-container div
# (or equivalent) so that the global layout only needs the fixed header and
# the page container.

app.layout = html.Div(
    [
        # -- Fixed header ---------------------------------------------------
        html.Header(
            [
                html.Div(
                    [
                        html.Span("\U0001F9EC", className="icon"),
                        html.Span(
                            "Interactive UMAP Explorer",
                            className="header-title",
                        ),
                    ]
                ),
                html.Nav(
                    [
                        dcc.Link("Single Dataset", href="/"),
                        dcc.Link("Compare", href="/comparison"),
                        dcc.Link("Integrated", href="/integrated"),
                    ],
                    className="header-nav",
                ),
            ],
            className="header",
        ),
        # -- Page content (switched by Dash multipage routing) --------------
        dash.page_container,
    ]
)

# ---------------------------------------------------------------------------
# Development server
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)

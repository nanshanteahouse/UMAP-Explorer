# =============================================================================
# Dockerfile — Interactive UMAP Explorer
#
# Base: python:3.11-slim (scanpy/numba compatibility)
#
# Production entrypoint: gunicorn with 4 workers
# Data directory (/app/data/) should be mounted as a volume.
# =============================================================================

FROM python:3.11-slim

# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
LABEL org.opencontainers.image.title="Interactive UMAP Explorer"
LABEL org.opencontainers.image.description="Dash-based UMAP visualization for scRNAseq datasets"
LABEL org.opencontainers.image.source="https://github.com/example/neurobiology"

# ---------------------------------------------------------------------------
# System dependencies
#   libgomp1 — required by scanpy / numba for OpenMP threading
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Non-root user
# ---------------------------------------------------------------------------
RUN adduser --disabled-password --gecos '' appuser

# ---------------------------------------------------------------------------
# Python dependencies
# ---------------------------------------------------------------------------
WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ---------------------------------------------------------------------------
# Application code
# ---------------------------------------------------------------------------
COPY app.py .
COPY assets/ assets/
COPY components/ components/
COPY pages/ pages/
COPY preprocessing/ preprocessing/

# Data should be mounted as a volume in production
# COPY data/ data/

# ---------------------------------------------------------------------------
# Switch to non-root user
# ---------------------------------------------------------------------------
RUN chown -R appuser:appuser /app
USER appuser

# ---------------------------------------------------------------------------
# Runtime configuration
# ---------------------------------------------------------------------------
EXPOSE 8050

HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8050/')"

# Run with gunicorn (4 workers, 120 s timeout for large UMAP rendering)
CMD ["gunicorn", "--workers=4", "--timeout=120", "--bind=0.0.0.0:8050", "app:server"]

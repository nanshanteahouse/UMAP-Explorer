#!/usr/bin/env bash
# ===========================================================================
# start.sh — Start Interactive UMAP Explorer with optional custom port.
#
# Usage:
#   ./start.sh                  # default port 8050
#   ./start.sh -p 8051          # custom port
#   ./start.sh --port 8051      # custom port
#   PORT=8051 ./start.sh        # via environment variable
# ===========================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEFAULT_PORT=8050
PORT="${PORT:-$DEFAULT_PORT}"

# ---- Parse CLI arguments ------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--port)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --port requires a value." >&2
                exit 1
            fi
            PORT="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [-p|--port PORT]"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 [-p|--port PORT]"
            exit 1
            ;;
    esac
done

# ---- Validate port ------------------------------------------------------
if ! [[ "$PORT" =~ ^[0-9]+$ ]] || [ "$PORT" -lt 1 ] || [ "$PORT" -gt 65535 ]; then
    echo "Error: Port must be a number between 1 and 65535, got: $PORT" >&2
    exit 1
fi

# ---- Activate virtual environment ---------------------------------------
VENV="$SCRIPT_DIR/../.venv/bin/activate"
if [ ! -f "$VENV" ]; then
    echo "Error: Virtual environment not found at $VENV" >&2
    exit 1
fi
# shellcheck disable=SC1091
source "$VENV"

# ---- Launch app ---------------------------------------------------------
echo "Starting Interactive UMAP Explorer on http://0.0.0.0:$PORT"
export PORT
python "$SCRIPT_DIR/app.py"

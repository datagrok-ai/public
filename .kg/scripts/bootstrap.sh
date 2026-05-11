#!/usr/bin/env bash
# Bootstrap the .kg/ Python venv and verify the graph works.
# Idempotent: safe to re-run.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KG_DIR="$(dirname "$SCRIPT_DIR")"
VENV="$KG_DIR/.venv"
REQS="$KG_DIR/requirements.txt"
DB="$KG_DIR/kg.kuzu"

# Pick the right python launcher
if command -v py >/dev/null 2>&1; then
    PY=py
elif command -v python3 >/dev/null 2>&1; then
    PY=python3
elif command -v python >/dev/null 2>&1; then
    PY=python
else
    echo "[bootstrap] no python found on PATH (need py / python3 / python)" >&2
    exit 1
fi

# Pick the right venv python path (Windows vs POSIX)
if [[ -f "$VENV/Scripts/python.exe" ]]; then
    VENV_PY="$VENV/Scripts/python.exe"
elif [[ -f "$VENV/bin/python" ]]; then
    VENV_PY="$VENV/bin/python"
else
    VENV_PY=""
fi

echo "[bootstrap] kg dir: $KG_DIR"

# 1. Create venv if missing
if [[ -z "$VENV_PY" ]]; then
    echo "[bootstrap] creating venv at $VENV ..."
    "$PY" -m venv "$VENV"
    if [[ -f "$VENV/Scripts/python.exe" ]]; then
        VENV_PY="$VENV/Scripts/python.exe"
    else
        VENV_PY="$VENV/bin/python"
    fi
else
    echo "[bootstrap] venv already exists, skipping creation"
fi

# 2. Install requirements
echo "[bootstrap] installing requirements ..."
"$VENV_PY" -m pip install --quiet --upgrade pip
"$VENV_PY" -m pip install --quiet -r "$REQS"

# 3. Build DB if missing
if [[ ! -e "$DB" ]]; then
    echo "[bootstrap] kg.kuzu/ missing — building from JSONL (~3 min) ..."
    "$VENV_PY" "$KG_DIR/build.py"
else
    echo "[bootstrap] kg.kuzu/ already exists, skipping build"
fi

# 4. Smoke test
echo "[bootstrap] smoke test ..."
"$VENV_PY" "$KG_DIR/qq.py" "MATCH (n) RETURN count(n) AS nodes LIMIT 1"
"$VENV_PY" "$KG_DIR/qq.py" --stop >/dev/null || true

cat <<EOF

[bootstrap] OK. Use:
  $VENV_PY $KG_DIR/qq.py "<cypher>"
EOF

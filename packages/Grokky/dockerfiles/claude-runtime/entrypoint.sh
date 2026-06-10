#!/usr/bin/env bash
# Container entrypoint: refresh the materialized knowledge-graph DB from the
# (possibly newer, git-synced) .kg-dist artifact, then start the runtime server.
set -u

KG_DIR="${CLAUDE_WORKSPACE:-/workspace}/.kg"
VENV_PY="$KG_DIR/.venv/bin/python"

if [[ -x "$VENV_PY" && -f "$KG_DIR/scripts/unpack.py" ]]; then
  # --force so a synced/updated kg.kuzu.xz is picked up on restart. Non-fatal.
  "$VENV_PY" "$KG_DIR/scripts/unpack.py" --force || echo "[entrypoint] kg unpack skipped"
else
  echo "[entrypoint] kg venv not found — knowledge graph unavailable"
fi

exec node dist/server.js

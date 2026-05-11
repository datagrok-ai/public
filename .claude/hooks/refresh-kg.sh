#!/usr/bin/env bash
# KG refresh script. Called explicitly from the /dg-task orchestrator at
# step 7, or by hand if you want to drain queues without /dg-task.
#
# Reads two queue files and acts on them:
#
#   .kg/.dirty-packages   — one package folder name per line
#   .kg/.learned/*.md     — newly discovered facts
#
# Idempotent: empty queues = no-op. Best-effort: missing venv = bail
# silently (so calling this from a fresh checkout doesn't error).
#
# NOT wired as a session-Stop hook (intentionally) — see
# .claude/settings.json. The /dg-task skill calls this directly so
# the graph only refreshes when the user asked for work.

set -u

# Anchor on script location so the hook works whether the project root
# is the monorepo root (with public/ as a subdir) or public/ itself.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KG="$(cd "$SCRIPT_DIR/../.." && pwd)/.kg"

# Locate the venv python (Windows or POSIX layout)
if [[ -x "$KG/.venv/Scripts/python.exe" ]]; then
    PY="$KG/.venv/Scripts/python.exe"
elif [[ -x "$KG/.venv/bin/python" ]]; then
    PY="$KG/.venv/bin/python"
else
    # No venv yet — silently skip
    exit 0
fi

DIRTY="$KG/.dirty-packages"
LEARNED_DIR="$KG/.learned"

# 1. Drain learned facts into the next enrichment slice
if [[ -d "$LEARNED_DIR" ]] && compgen -G "$LEARNED_DIR/*.md" >/dev/null; then
    "$PY" "$KG/scripts/learned_to_enrichment.py" 2>/dev/null || true
fi

# 2. Refresh dirty packages
if [[ -s "$DIRTY" ]]; then
    # Dedupe and trim
    PKGS=$(sort -u "$DIRTY" | grep -v '^$' | paste -sd, -)
    if [[ -n "$PKGS" ]]; then
        # Stop server first to release the DB lock
        "$PY" "$KG/qq.py" --stop >/dev/null 2>&1 || true
        # Re-extract just the dirty packages, then rebuild Kuzu
        "$PY" "$KG/build.py" --packages "$PKGS" >/dev/null 2>&1 || true
    fi
    # Truncate the queue
    : > "$DIRTY"
fi

exit 0

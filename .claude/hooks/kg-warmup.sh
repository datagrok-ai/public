#!/usr/bin/env bash
# KG server warmup. Called explicitly from the /dg-task orchestrator at
# step 0, or by hand if you want a warm server for ad-hoc qq.py use.
#
# NOT wired as a SessionStart hook (intentionally) — see
# .claude/settings.json. The /dg-task skill calls this directly so
# the kg_server only spawns when the user asked for /dg-task work.
#
# Best-effort: silent skip if the venv or DB isn't there.

set -u

# Anchor on script location so the hook works whether the project root
# is the monorepo root (with public/ as a subdir) or public/ itself.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KG="$(cd "$SCRIPT_DIR/../.." && pwd)/.kg"

if [[ -x "$KG/.venv/Scripts/python.exe" ]]; then
    PY="$KG/.venv/Scripts/python.exe"
elif [[ -x "$KG/.venv/bin/python" ]]; then
    PY="$KG/.venv/bin/python"
else
    exit 0
fi

# Skip if DB doesn't exist yet (fresh clone) — bootstrap will handle it
[[ -e "$KG/kg.kuzu" ]] || exit 0

# Already running?
if "$PY" "$KG/qq.py" --no-start >/dev/null 2>&1; then
    exit 0
fi

# Start it in the background, swallow output
nohup "$PY" "$KG/kg_server.py" >/dev/null 2>&1 &
disown 2>/dev/null || true

exit 0

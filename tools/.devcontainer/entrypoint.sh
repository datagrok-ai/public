#!/bin/bash
set -e

# Map DG_VERSION to a git branch name for cloning the public repo.
# Docker tags like "latest" or "bleeding-edge" aren't git branches — use master.
# Explicit version numbers (e.g. "0.151.9") are tried as-is with a master fallback.
resolve_branch() {
  local ver="${1:-latest}"
  case "$ver" in
    latest|bleeding-edge) echo "master" ;;
    *) echo "$ver" ;;
  esac
}

PUBLIC_DIR="${DG_PUBLIC_DIR:-/workspace/datagrok}"
# Detect if /workspace/repo IS the public repo.
# Require .git to avoid false positives from packages that have public/js-api via npm link.
if [ -e "/workspace/repo/.git" ] && [ -d "/workspace/repo/js-api" ]; then
  echo "[tools-dev] Public repo detected at /workspace/repo — skipping clone."
elif [ -e "/workspace/repo/.git" ] && [ -d "/workspace/repo/public/js-api" ]; then
  echo "[tools-dev] Monorepo detected — public repo at /workspace/repo/public."
elif [ -d "$PUBLIC_DIR/js-api" ]; then
  echo "[tools-dev] Public repo already at $PUBLIC_DIR."
else
  REPO="${DG_PUBLIC_REPO:-https://github.com/datagrok-ai/public.git}"
  # Resolve branch: explicit DG_PUBLIC_BRANCH > resolved DG_VERSION > master
  if [ -n "$DG_PUBLIC_BRANCH" ]; then
    BRANCH="$DG_PUBLIC_BRANCH"
  else
    BRANCH=$(resolve_branch "${DG_VERSION:-latest}")
  fi
  # Sparse checkout: only fetch dirs needed for package context (js-api, libraries, help, packages)
  # plus root files (CLAUDE.md, .claude/, etc.). Much faster than a full clone.
  sparse_clone() {
    local branch="$1"
    git clone --depth 1 --branch "$branch" --filter=blob:none --no-checkout "$REPO" "$PUBLIC_DIR" \
      && git -C "$PUBLIC_DIR" sparse-checkout set --no-cone \
           '/*' '!connectors/' '!docker/' '!docusaurus/' '!docusaurus-static/' \
           '!environments/' '!hooks/' '!misc/' '!python-api/' '!datagrok-celery-task/' \
           '/js-api/**' '/libraries/**' '/help/**' '/packages/**' '/tools/**' \
      && git -C "$PUBLIC_DIR" checkout
  }
  echo "[tools-dev] Cloning public repo ($BRANCH, sparse) into $PUBLIC_DIR..."
  sparse_clone "$BRANCH" \
    || { echo "[tools-dev] Branch '$BRANCH' not found, falling back to master."
         rm -rf "$PUBLIC_DIR"
         sparse_clone "master"; }
  echo "[tools-dev] Public repo ready at $PUBLIC_DIR (branch: $(git -C "$PUBLIC_DIR" branch --show-current))."
fi

# ── Mount workspace inside public repo for non-public workspaces ──
# When TASK_KEY is set and workspace is not the public repo, link it into packages/
# so Claude Code walks up to find all public repo context (CLAUDE.md, .claude/, js-api/, etc.)
if [ -e "/workspace/repo/.git" ] && [ -d "/workspace/repo/js-api" ]; then
  echo "[tools-dev] Public repo IS the workspace — no linking needed."
elif [ -e "/workspace/repo/.git" ] && [ -d "/workspace/repo/public/js-api" ]; then
  echo "[tools-dev] Monorepo detected — public context at /workspace/repo/public/."
elif [ -d "$PUBLIC_DIR/js-api" ]; then
  # Cloned public repo — link workspace into packages/ and expose context at /workspace/
  LINK_NAME="${FOLDER_NAME:-$TASK_KEY}"
  if [ -n "$LINK_NAME" ]; then
    mkdir -p "$PUBLIC_DIR/packages" 2>/dev/null || true
    [ ! -e "$PUBLIC_DIR/packages/$LINK_NAME" ] && ln -s /workspace/repo "$PUBLIC_DIR/packages/$LINK_NAME"
    echo "[tools-dev] Workspace linked at $PUBLIC_DIR/packages/$LINK_NAME"
  fi
  PUBLIC_BASENAME=$(basename "$PUBLIC_DIR")
  [ ! -e /workspace/.claude ] && ln -s "$PUBLIC_BASENAME/.claude" /workspace/.claude
  [ ! -e /workspace/CLAUDE.md ] && ln -s "$PUBLIC_BASENAME/CLAUDE.md" /workspace/CLAUDE.md
  echo "[tools-dev] Linked public repo context at /workspace/"
fi

# Auto-configure grok CLI to point to the local Datagrok instance.
# Only creates the config if it doesn't already exist (preserves host config).
GROK_DIR="/home/node/.grok"
GROK_CFG="$GROK_DIR/config.yaml"
if [ ! -f "$GROK_CFG" ] || [ ! -s "$GROK_CFG" ]; then
  mkdir -p "$GROK_DIR" 2>/dev/null || true
  cat > "$GROK_CFG" <<'YAML' 2>/dev/null || true
default: local
servers:
  local:
    url: http://datagrok:8080/api
    key: admin
YAML
  [ -s "$GROK_CFG" ] && echo "[tools-dev] Created grok config at $GROK_CFG"
elif ! grep -q "datagrok:8080" "$GROK_CFG" 2>/dev/null; then
  echo "[tools-dev] Note: existing grok config found. Add 'local' server with:"
  echo "  grok config add --alias local --server http://datagrok:8080/api --key admin --default"
fi

exec "$@"

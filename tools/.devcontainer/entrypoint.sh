#!/bin/bash
set -e

# Resolve the actual version when DG_VERSION is "latest" — needed for git branch matching.
resolve_version() {
  local ver="${1:-latest}"
  if [ "$ver" = "latest" ]; then
    echo "[tools-dev] Resolving latest Datagrok version from Docker Hub..." >&2
    local resolved
    resolved=$(curl -sf "https://hub.docker.com/v2/repositories/datagrok/datagrok/tags?page_size=100&ordering=-name" \
      | node -e "
        const data = JSON.parse(require('fs').readFileSync('/dev/stdin','utf8'));
        const tags = (data.results||[]).map(t=>t.name).filter(n=>/^\d+\.\d+\.\d+$/.test(n));
        tags.sort((a,b)=>{
          const pa=a.split('.').map(Number), pb=b.split('.').map(Number);
          return pb[0]-pa[0]||pb[1]-pa[1]||pb[2]-pa[2];
        });
        if(tags[0]) process.stdout.write(tags[0]);
      " 2>/dev/null) || true
    if [ -n "$resolved" ]; then
      echo "[tools-dev] Resolved latest version: $resolved" >&2
      echo "$resolved"
    else
      echo "[tools-dev] Could not resolve latest version, falling back to master" >&2
      echo "master"
    fi
  else
    echo "$ver"
  fi
}

PUBLIC_DIR="/workspace/public"
# Detect if /workspace/repo IS the public repo (has js-api/ at root)
if [ -d "/workspace/repo/js-api" ]; then
  echo "[tools-dev] Public repo detected at /workspace/repo — skipping clone."
elif [ -d "/workspace/repo/public/js-api" ]; then
  echo "[tools-dev] Monorepo detected — public repo at /workspace/repo/public."
elif [ -d "$PUBLIC_DIR/js-api" ]; then
  echo "[tools-dev] Public repo already at $PUBLIC_DIR."
else
  REPO="${DG_PUBLIC_REPO:-https://github.com/datagrok-ai/public.git}"
  # Resolve branch: explicit DG_PUBLIC_BRANCH > resolved DG_VERSION > master
  if [ -n "$DG_PUBLIC_BRANCH" ]; then
    BRANCH="$DG_PUBLIC_BRANCH"
  else
    BRANCH=$(resolve_version "${DG_VERSION:-latest}")
  fi
  echo "[tools-dev] Cloning public repo ($BRANCH) into $PUBLIC_DIR..."
  git clone --depth 1 --branch "$BRANCH" "$REPO" "$PUBLIC_DIR" 2>/dev/null \
    || { echo "[tools-dev] Branch '$BRANCH' not found, falling back to master."
         git clone --depth 1 --branch master "$REPO" "$PUBLIC_DIR"; }
  echo "[tools-dev] Public repo ready at $PUBLIC_DIR (branch: $(git -C "$PUBLIC_DIR" branch --show-current))."
fi

# ── Mount workspace inside public repo for non-public workspaces ──
# When TASK_KEY is set and workspace is not the public repo, link it into packages/
# so Claude Code walks up to find all public repo context (CLAUDE.md, .claude/, js-api/, etc.)
if [ -d "/workspace/repo/js-api" ]; then
  echo "[tools-dev] Public repo IS the workspace — no linking needed."
elif [ -d "/workspace/repo/public/js-api" ]; then
  echo "[tools-dev] Monorepo detected — public context at /workspace/repo/public/."
elif [ -d "$PUBLIC_DIR/js-api" ]; then
  # Cloned public repo — link workspace into packages/ and expose context at /workspace/
  if [ -n "$TASK_KEY" ]; then
    mkdir -p "$PUBLIC_DIR/packages" 2>/dev/null || true
    [ ! -e "$PUBLIC_DIR/packages/$TASK_KEY" ] && ln -s /workspace/repo "$PUBLIC_DIR/packages/$TASK_KEY"
    echo "[tools-dev] Workspace linked at $PUBLIC_DIR/packages/$TASK_KEY"
  fi
  [ ! -e /workspace/.claude ] && ln -s public/.claude /workspace/.claude
  [ ! -e /workspace/CLAUDE.md ] && ln -s public/CLAUDE.md /workspace/CLAUDE.md
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

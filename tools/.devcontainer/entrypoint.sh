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
# Ensure workspace dirs are writable by node user (bind-mounts create parents as root)
chown node:node /workspace "$PUBLIC_DIR" 2>/dev/null || true
# Allow git to operate on dirs owned by different users (root runs git, dir owned by node)
git config --global --add safe.directory "$PUBLIC_DIR" 2>/dev/null || true
git config --global init.defaultBranch master 2>/dev/null || true
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
  # Sparse checkout (cone mode) with partial clone: only fetch js-api, libraries, and
  # ApiSamples. Cone mode integrates with --filter=blob:none so the server only sends
  # blobs for the included directories (~3 MB vs 1.67 GB for the full tree).
  sparse_clone() {
    local branch="$1"
    # Use init+fetch instead of clone to handle pre-existing directories (e.g. mount points)
    git init -q "$PUBLIC_DIR" \
      && (git -C "$PUBLIC_DIR" remote add origin "$REPO" 2>/dev/null || git -C "$PUBLIC_DIR" remote set-url origin "$REPO") \
      && git -C "$PUBLIC_DIR" config remote.origin.promisor true \
      && git -C "$PUBLIC_DIR" config remote.origin.partialclonefilter blob:none \
      && git -C "$PUBLIC_DIR" sparse-checkout set --cone .claude js-api libraries packages/ApiSamples \
      && git -C "$PUBLIC_DIR" fetch -q --depth 1 --filter=blob:none origin "$branch" \
      && git -C "$PUBLIC_DIR" checkout -q -B "$branch" FETCH_HEAD
  }
  # Clear git state without touching bind-mounted subdirs (e.g. packages/awesome)
  clear_git() {
    rm -rf "$1/.git" 2>/dev/null || true
  }
  echo "[tools-dev] Cloning public repo ($BRANCH, sparse) into $PUBLIC_DIR..."
  if sparse_clone "$BRANCH"; then
    echo "[tools-dev] Public repo ready."
  elif [ "$BRANCH" != "master" ]; then
    echo "[tools-dev] Branch '$BRANCH' clone failed, falling back to master."
    clear_git "$PUBLIC_DIR"
    if sparse_clone "master"; then
      echo "[tools-dev] Public repo ready (branch: master)."
    else
      echo "[tools-dev] WARNING: Failed to clone public repo."
    fi
  else
    echo "[tools-dev] WARNING: Failed to clone public repo."
  fi
  # Ensure node user can read cloned files (entrypoint runs as root)
  chown -R node:node "$PUBLIC_DIR" 2>/dev/null || true
fi

# ── Mount workspace inside public repo for non-public workspaces ──
# When TASK_KEY is set and workspace is not the public repo, link it into packages/
# so Claude Code walks up to find all public repo context (CLAUDE.md, .claude/, js-api/, etc.)
if [ -e "/workspace/repo/.git" ] && [ -d "/workspace/repo/js-api" ]; then
  echo "[tools-dev] Public repo IS the workspace — no linking needed."
elif [ -e "/workspace/repo/.git" ] && [ -d "/workspace/repo/public/js-api" ]; then
  echo "[tools-dev] Monorepo detected — public context at /workspace/repo/public/."
elif [ -d "$PUBLIC_DIR/js-api" ]; then
  LINK_NAME="${FOLDER_NAME:-$TASK_KEY}"
  if [ -n "$LINK_NAME" ]; then
    mkdir -p "$PUBLIC_DIR/packages" 2>/dev/null || true
    if [ ! -d "$PUBLIC_DIR/packages/$LINK_NAME" ]; then
      ln -sfn /workspace/repo "$PUBLIC_DIR/packages/$LINK_NAME"
    fi
    echo "[tools-dev] Workspace at $PUBLIC_DIR/packages/$LINK_NAME"
  fi
  PUBLIC_BASENAME=$(basename "$PUBLIC_DIR")
  ln -sfn "$PUBLIC_BASENAME/.claude" /workspace/.claude 2>/dev/null || true
  ln -sfn "$PUBLIC_BASENAME/CLAUDE.md" /workspace/CLAUDE.md 2>/dev/null || true
fi

# Auto-configure grok CLI to point to the local Datagrok instance.
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
fi

# Drop to node user for the main process (entrypoint runs as root for permission fixes)
if [ "$(id -u)" = "0" ]; then
  exec setpriv --reuid=node --regid=node --init-groups "$@"
fi
exec "$@"

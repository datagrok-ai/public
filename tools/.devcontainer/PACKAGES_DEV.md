# Packages Dev Container

Isolated environment for JS/TS package development against a running Datagrok instance.
All images are pre-built on Docker Hub — no local build step needed. Install `datagrok-tools`
globally and run `grok claude` from any directory.

## Quick start

```bash
npm i -g datagrok-tools        # one-time install

# From any package directory:
grok claude GROK-12345                        # start working on a task
grok claude GROK-12345 --version 1.22.0       # pin Datagrok version
grok claude GROK-12345 --profile full --keep  # all services, keep running
grok claude GROK-12345 --prompt "fix the bug" # one-shot command
grok claude GROK-12345 --in-place             # use current directory (no worktree)
grok claude destroy GROK-12345                # tear down a task
grok claude destroy-all                       # tear down everything
```

Authentication: Claude Code credentials are copied from the host `~/.claude/` into the
container on startup. Alternatively, set `ANTHROPIC_API_KEY` in the environment.

**Project name restrictions:** `master` and `main` are rejected.

### What it does

1. Creates a git worktree at `~/pkg-worktrees/<project>` (unless `--in-place` or not in a git repo)
2. Generates `docker-compose.yaml`, `.env`, and optional `docker-compose.override.yaml` in `$TMPDIR/dg-pkg-<project>`
3. Runs `docker compose up -d --wait` (pulls pre-built images from Docker Hub)
4. Copies Claude credentials into the container and fixes ownership
5. Detects working directory based on repo type (see below)
6. Launches `claude --dangerously-skip-permissions` inside the `tools-dev` container
7. On exit, stops containers (unless `--keep`)

### Manual setup

For manual setup or customization, use the `docker-compose.yaml` in this directory
directly. The `Dockerfile.pkg_dev` and `entrypoint.sh` are the build recipe for the
`datagrok/tools-dev` Docker Hub image — end users don't need them.

## Architecture

```
Host machine
┌──────────────────────────────────────────────────────────────────┐
│  Docker network: dg-pkg-<project>-net                            │
│  ┌────────────────────────────────────────────────────────────┐  │
│  │  datagrok     (datagrok/datagrok:${VERSION})   → :${PORT} │  │
│  │  postgres     (pgvector/pgvector:pg17)                     │  │
│  │  rabbitmq     (rabbitmq:4.0.5-management)                  │  │
│  │  grok_pipe    (datagrok/grok_pipe:latest)                  │  │
│  │  grok_connect (datagrok/grok_connect:latest)               │  │
│  │  grok_spawner, jkg          (optional profiles)            │  │
│  │  world, test_db, northwind  (optional demo DBs)            │  │
│  │                                                            │  │
│  │  tools-dev  (datagrok/tools-dev:latest, runs as root)      │  │
│  │    ├── /workspace/repo       ← bind mount of worktree     │  │
│  │    ├── /workspace/datagrok/  ← sparse clone of public repo │  │
│  │    │   ├── .claude/, CLAUDE.md                             │  │
│  │    │   ├── js-api/                                         │  │
│  │    │   ├── libraries/                                      │  │
│  │    │   └── packages/                                       │  │
│  │    │       ├── ApiSamples/   ← checked out                 │  │
│  │    │       └── <folder>/     ← bind mount of user's repo   │  │
│  │    ├── /var/run/docker.sock  ← host Docker for stand mgmt  │  │
│  │    ├── entrypoint.sh         ← bind-mounted from tools/    │  │
│  │    ├── Claude Code + MCP plugins (Jira, GitHub)            │  │
│  │    ├── Playwright + Chrome for testing                     │  │
│  │    └── grok CLI for build/publish/test                     │  │
│  └────────────────────────────────────────────────────────────┘  │
│                                                                   │
│  Host browser → http://localhost:${PORT} (Datagrok UI)            │
│               → chrome://inspect → localhost:9222 (Debug)         │
└───────────────────────────────────────────────────────────────────┘
```

## Working directory detection

Claude Code launches in different directories depending on the repo type:

| Repo type | Working directory | How detected |
|-----------|------------------|--------------|
| Public repo root | `/workspace/repo` | Host has `js-api/` at git root |
| Monorepo | `/workspace/repo/public` | Host has `public/js-api/` at git root |
| External repo | `/workspace/datagrok/packages/<folder>` | Everything else |

For **external repos**, the user's directory is bind-mounted at both `/workspace/repo`
and `/workspace/datagrok/packages/<folder>`. Claude Code starts in the latter so it can
walk up the directory tree to find `CLAUDE.md`, `.claude/`, and `js-api/` from the
sparse-cloned public repo.

## Sparse clone (external repos)

When the workspace is not the public repo, the entrypoint sparse-clones it into
`/workspace/datagrok/` using **cone mode** with `--filter=blob:none`:

```
git sparse-checkout set --cone .claude js-api libraries packages/ApiSamples
```

This fetches only ~3 MB (vs 1.67 GB for the full repo) and completes in ~10 seconds.
The checked-out directories provide:
- `.claude/` — Claude Code skills and settings
- `CLAUDE.md` — project instructions for the agent
- `js-api/` — JS API source for reference
- `libraries/` — shared library source
- `packages/ApiSamples/` — runnable code examples

### Branch resolution

1. `DG_PUBLIC_BRANCH` (explicit override)
2. Current branch of the host repo (if inside the public repo)
3. `DG_VERSION` mapped to a branch (version tags like `latest`/`bleeding-edge` → `master`)
4. Falls back to `master`

To use a private fork: `DG_PUBLIC_REPO=https://github.com/myorg/public-fork.git`

## Entrypoint

The entrypoint (`entrypoint.sh`) runs as **root** (via `user: root` in compose) to
handle bind-mount permission issues, then drops to the `node` user (UID 1000) via
`setpriv` for the main process (`sleep infinity`).

Steps:
1. Fix ownership of `/workspace/datagrok` (bind-mounts create parents as root)
2. Detect repo type (public, monorepo, or external)
3. Sparse-clone public repo if needed (cone mode, ~10s)
4. `chown -R node:node` on cloned files
5. Link workspace into `packages/` (skipped if bind-mount already exists)
6. Create symlinks at `/workspace/` for CLAUDE.md and .claude/
7. Auto-configure `~/.grok/config.yaml` pointing to `http://datagrok:8080/api`
8. Drop privileges: `exec setpriv --reuid=node --regid=node --init-groups "$@"`

The entrypoint is **bind-mounted from the host repo** (`tools/.devcontainer/entrypoint.sh`)
via the compose override, so changes take effect without rebuilding the Docker image.

## Compose template

The compose configuration is embedded in `bin/commands/claude.ts` — no external files
needed. It generates three files in `$TMPDIR/dg-pkg-<project>/`:

- `docker-compose.yaml` — main compose template
- `.env` — resolved variables (paths, versions, ports, tokens)
- `docker-compose.override.yaml` — host-specific volumes (entrypoint, credentials mount for external repos)

### Key differences from `.devcontainer/docker-compose.yaml`

| | Embedded (grok claude) | .devcontainer/ |
|---|---|---|
| Version vars | Separate per service (`DATAGROK_VERSION`, `GROK_CONNECT_VERSION`, etc.) | Single `DG_VERSION` |
| Workspace mount | `/workspace/repo` | `/workspace` |
| grok_connect | Always on | Profile `full` only |
| grok_pipe | `grok_pipe:latest` | `grok_pipe:${DG_VERSION}` |
| Host ~/.grok | Not mounted (auto-created by entrypoint) | Mounted from host |
| Entrypoint | Bind-mounted from repo via override | Baked into image |
| User | `root` (drops to `node` via setpriv) | `node` (from Dockerfile) |

## Credential handling

Claude Code credentials are **copied** into the container after startup (not
bind-mounted) to avoid permission issues — host files may have different UID/GID
than the container's `node` user.

Files copied from `~/.claude/`:
- `.credentials.json` — OAuth credentials
- `settings.json` — Claude Code settings
- `settings.local.json` — local settings

After copying, ownership is fixed: `chown -R node:node /home/node/.claude`

The host `~/.claude` directory is found by searching (in order):
1. `CLAUDE_HOME` environment variable
2. `$HOME/.claude`
3. `$USERPROFILE/.claude` (Windows)
4. `$APPDATA/../.claude` (Windows fallback)

## Profiles

| Profile | Services added | Use case |
|---------|---------------|----------|
| (none) | postgres, rabbitmq, grok_pipe, datagrok, grok_connect, tools-dev | Basic package dev |
| `scripting` | + jupyter_kernel_gateway | Python/R/Julia scripts |
| `demo` | + world, test_db, northwind | Demo databases |
| `full` | + grok_spawner, JKG, demo DBs | Everything |

## tools-dev container

Based on `node:22-bookworm-slim`. Pre-installed:
- Google Chrome stable (for Puppeteer)
- Playwright + Chromium
- `datagrok-tools` (grok CLI) — global npm
- `@anthropic-ai/claude-code` — global npm
- git, curl, jq, docker CLI
- Node.js 22

Runs as **root** in compose (entrypoint drops to `node` user for the main process).

## MCP plugins: Jira and GitHub

The agent uses MCP plugins to search for similar issues, read context, and update
status. Set tokens in `.env` or the agent will ask interactively on first use.

### Jira (Atlassian)

```bash
# Inside tools-dev, register the MCP server:
claude mcp add mcp-atlassian -s user -- \
  npx -y mcp-atlassian \
  --jira-url "$JIRA_URL" \
  --jira-username "$JIRA_USERNAME" \
  --jira-token "$JIRA_TOKEN"
```

### GitHub

```bash
claude mcp add github -s user -- \
  npx -y @modelcontextprotocol/server-github
```

Requires `GITHUB_TOKEN` in the environment (already passed from `.env`).

## Testing

### grok test (Puppeteer-based)

```bash
cd /workspace/datagrok/packages/MyPkg
grok test --host local              # headless
grok test --host local --gui        # visible browser
grok test --host local --verbose    # detailed output
```

### Playwright

```bash
cd /workspace/datagrok/packages/MyPkg
npx playwright test --project chromium
```

### Chrome remote debugging

`grok test --gui` runs Chrome with `--remote-debugging-port=9222`.
On host: `chrome://inspect` → Configure → add `localhost:9222`

## Publishing packages

```bash
cd /workspace/datagrok/packages/MyPkg
grok publish                    # debug mode (visible only to dev)
grok publish --release          # public release
grok publish --build            # build webpack first
```

### To an external server

```bash
grok config add --alias prod \
  --server https://example.datagrok.ai/api \
  --key <dev-key>
grok publish prod --release --build
```

## Exposing the client

- Datagrok UI: `http://localhost:${DG_PORT}` (port shown at startup, or use `--port` to fix it)
- Default credentials: **admin / admin**
- grok CLI inside container: auto-configured to `http://datagrok:8080/api` with key `admin`

## Troubleshooting

### Datagrok not starting
```bash
docker compose logs datagrok
# Wait for DB migrations — first start takes 1-2 minutes
```

### grok publish fails
Ensure `~/.grok/config.yaml` inside the container points to `http://datagrok:8080/api`.
The container reaches Datagrok via Docker DNS (`datagrok`), not `localhost`.

### Permission denied on /workspace/datagrok
The entrypoint runs as root and fixes ownership. If using the image directly (not via
`grok claude`), add `user: root` to the compose service definition.

### Agent can't manage Docker
The Docker socket must be mounted. Check: `docker exec <tools-dev> docker ps`

### Container can't resolve `datagrok` hostname
All services must be on the same Docker network:
```bash
docker network inspect dg-pkg-<project>-net
```

# Packages Dev Container

Isolated environment for JS/TS package development against a running Datagrok instance.
All images are pre-built on Docker Hub — no local build step needed. Install `datagrok-tools`
globally and run `grok claude` from any git repo.

## Quick start

```bash
npm i -g datagrok-tools        # one-time install
export ANTHROPIC_API_KEY=...   # or set in shell profile

# From any package directory:
grok claude

# With options:
grok claude --version 1.22.0 --profile full --task GROK-123

# Stop containers:
grok claude --stop --task GROK-123
```

The `grok claude` command is fully self-contained — the compose configuration is
embedded in the CLI, all container images are pulled from Docker Hub (`datagrok/tools-dev`,
`datagrok/datagrok`, etc.). No files to copy, no Dockerfile to build.
See `grok claude --help` for all options.

### What it does

1. Generates a `docker-compose.yaml` and `.env` in a temp directory
2. Runs `docker compose up -d --wait` (pulls pre-built images)
3. Waits for Datagrok to be healthy
4. Launches `claude --dangerously-skip-permissions` inside the `tools-dev` container
5. On exit, stops containers (unless `--keep`)

### Manual setup

For manual setup or customization, use the `docker-compose.yaml` in this directory
directly. The `Dockerfile.pkg_dev` and `entrypoint.sh` are the build recipe for the
`datagrok/tools-dev` Docker Hub image — end users don't need them.

## What you can do

| Workflow | How it works |
|----------|-------------|
| **Develop packages (public repo)** | Mount public repo worktree as `/workspace`. Agent reads js-api source, help docs, ApiSamples, and CLAUDE.md skills directly. |
| **Develop packages (separate repo)** | Mount your repo as `/workspace`. The public repo is auto-cloned to `/workspace/public` on first start, branch matching `DG_VERSION`. |
| **Learn from examples** | Agent reads `packages/ApiSamples/scripts/` for runnable code samples and `help/develop/` for guides. |
| **Search Jira and GitHub** | MCP plugins for Atlassian Jira and GitHub — agent searches for similar issues, reads context, updates status. |
| **Write and run tests** | Playwright for browser automation, `grok test` for Puppeteer-based package tests, Chrome remote debugging. |
| **Publish packages** | `grok publish` to the local Datagrok instance or any external server. Agent handles build + publish. |
| **Manage the Datagrok stand** | Agent changes `DG_VERSION`, pulls new images, resets DB, redeploys — all from inside the container via Docker socket. |
| **Interactive setup** | Agent asks the user for tokens, auth keys, server URLs via Claude Code prompts. No pre-configuration required. |

## Architecture

```
Host machine
┌──────────────────────────────────────────────────────────────────┐
│  Worktree: ~/pkg-worktrees/TASK-123/                            │
│    └── public/ (auto-cloned if workspace is not public repo)    │
│                                                                  │
│  Docker network: dg-pkg-TASK-123-net                             │
│  ┌────────────────────────────────────────────────────────────┐  │
│  │  datagrok    (datagrok/datagrok:${VERSION})    → :8080     │  │
│  │  postgres    (pgvector/pgvector:pg17)                      │  │
│  │  rabbitmq    (rabbitmq:4.0.5-management)                   │  │
│  │  grok_pipe   (datagrok/grok_pipe:${VERSION})               │  │
│  │  grok_connect, grok_spawner, jkg (optional profiles)       │  │
│  │  world, test_db, northwind (optional demo DBs)             │  │
│  │                                                             │  │
│  │  tools-dev     (datagrok/tools-dev:latest)                     │  │
│  │    ├── /workspace          ← bind mount of worktree        │  │
│  │    ├── /var/run/docker.sock ← host Docker for stand mgmt   │  │
│  │    ├── ~/.claude/          ← bind mount from host          │  │
│  │    ├── Claude Code + MCP plugins (Jira, GitHub)            │  │
│  │    ├── Playwright + Chrome for testing                     │  │
│  │    └── grok CLI for build/publish/test                     │  │
│  └────────────────────────────────────────────────────────────┘  │
│                                                                  │
│  Host browser → http://localhost:8080 (Datagrok UI)              │
│               → chrome://inspect → localhost:9222 (Debug)        │
└──────────────────────────────────────────────────────────────────┘
```

## docker-compose.yaml

The `docker-compose.yaml` in this directory uses pre-built images from Docker Hub.
The `grok claude` command embeds this same configuration and generates it automatically
— you don't need this file unless you want manual control.

For manual use, set variables in `.env` or export them, then:

```bash
docker compose up -d
```

### Building the tools-dev image

The `Dockerfile.pkg_dev` and `entrypoint.sh` are the build recipe for the
`datagrok/tools-dev` image published to Docker Hub. To build locally:

```bash
cd public/tools/.devcontainer
docker build -f Dockerfile.pkg_dev -t datagrok/tools-dev:latest .
```

## Working with repos

### Public repo (default)

```bash
cd /path/to/public
git worktree add ~/pkg-worktrees/TASK-123 -b TASK-123

# Set WORKTREE_PATH in .env or export
echo "WORKTREE_PATH=$HOME/pkg-worktrees/TASK-123" >> .env
docker compose up -d
```

Inside the container, the agent sees:
- `/workspace/js-api/` — JS API source (read directly, no build needed for reference)
- `/workspace/packages/ApiSamples/scripts/` — runnable code examples
- `/workspace/help/develop/` — development guides
- `/workspace/packages/` — all existing packages as reference

### Separate repo (e.g., private packages)

```bash
git -C /path/to/my-repo worktree add ~/pkg-worktrees/TASK-123 -b TASK-123

echo "WORKTREE_PATH=$HOME/pkg-worktrees/TASK-123" >> .env
docker compose up -d
```

On first start, the entrypoint detects that `/workspace` is not the public repo
(no `js-api/` at root) and automatically clones it to `/workspace/public`:

```
[tools-dev] Cloning public repo (bleeding-edge) into /workspace/public...
[tools-dev] Public repo ready at /workspace/public (branch: bleeding-edge).
```

The branch is resolved in this order:
1. `DG_PUBLIC_BRANCH` (explicit override in `.env`)
2. `DG_VERSION` (matches the Datagrok image tag — keeps API and server in sync)
3. Falls back to `master` if the branch/tag doesn't exist

Inside the container:
- `/workspace/` — your repo
- `/workspace/public/js-api/` — JS API source (matching Datagrok version)
- `/workspace/public/packages/ApiSamples/scripts/` — code examples
- `/workspace/public/help/` — docs
- `/workspace/public/packages/` — reference packages

To use a private fork of public, set `DG_PUBLIC_REPO` in `.env`:
```bash
DG_PUBLIC_REPO=https://github.com/myorg/public-fork.git
```

Worktrees keep full git connectivity — commit, push, fetch, PR all work.

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

The agent can then:
- Search for issues: "find similar bugs to GROK-12345"
- Read issue details, comments, attachments
- Update issue status, add comments
- Create new issues

### GitHub

```bash
# Inside tools-dev, register the MCP server:
claude mcp add github -s user -- \
  npx -y @modelcontextprotocol/server-github
```

Requires `GITHUB_TOKEN` in the environment (already passed from `.env`). The agent can:
- Search issues and PRs across repos
- Read issue comments and PR diffs
- Create issues, PRs, and comments
- Check CI status

### Interactive token setup

If tokens are not pre-configured, the agent asks the user:

```
Agent: I need a Jira API token to search for similar issues.
       Go to https://id.atlassian.com/manage-profile/security/api-tokens
       and create a token. Paste it here.
User:  <pastes token>
Agent: <configures MCP server and proceeds>
```

The same flow works for GitHub tokens and grok dev keys.

## Testing

### Playwright (browser automation)

Playwright is pre-installed with Chromium. Use it for custom browser automation,
E2E tests, and UI interaction scenarios.

```bash
# Inside tools-dev:
cd /workspace/packages/MyPkg

# Run Playwright tests (if the package has them)
npx playwright test --project chromium

# Or write ad-hoc automation
npx playwright codegen http://datagrok:8080
```

### grok test (Puppeteer-based package tests)

Standard Datagrok package testing via the `grok` CLI:

```bash
# Inside tools-dev:
cd /workspace/packages/MyPkg
grok test --host local              # headless
grok test --host local --gui        # visible browser
grok test --host local --verbose    # detailed output
grok test --host local --category "MyCategory"  # filter tests
```

### Chrome remote debugging

`grok test --gui` runs Chrome with `--remote-debugging-port=9222`. Port 9222 is
exposed to the host.

On host: open `chrome://inspect` → Configure → add `localhost:9222` → click "inspect"
on the test session. Set breakpoints in package source during test execution.

### Developing test scenarios

The agent can:
1. Read existing tests in `packages/ApiTests/` and other packages for patterns
2. Create new test files following the `package-test.ts` template
3. Run tests and analyze failures
4. Generate test cases from Jira issue descriptions or GitHub issues

```bash
# Scaffold a test file in a package
cd /workspace/packages/MyPkg
grok add test
```

## Publishing packages

### To the local Datagrok instance

```bash
# Inside tools-dev:
cd /workspace/packages/MyPkg
grok publish                    # debug mode (visible only to dev)
grok publish --release          # public release
grok publish --build            # build webpack first
```

### To an external server

```bash
# Add a server config
grok config add --alias prod \
  --server https://example.datagrok.ai/api \
  --key <dev-key>

# Publish to it
grok publish prod --release --build
```

The agent handles the full cycle: build, check, publish, verify in the UI.

## Managing the Datagrok stand

The tools-dev container has the Docker socket mounted, so the agent can manage the
entire compose stack from inside.

### Version switching

```bash
# Inside tools-dev (or from host):
DG_VERSION=1.22.0 docker compose up -d --pull always
```

When the workspace is a separate repo (not public), the auto-cloned public repo
should also be updated to match the new version:

```bash
# Inside tools-dev — update the public clone to match new DG_VERSION:
cd /workspace/public && git fetch && git checkout 1.22.0
```

The agent does this autonomously when asked:
```
User:  Switch to Datagrok 1.22.0
Agent: <updates DG_VERSION, pulls new images, updates public branch>
       Datagrok 1.22.0 is running. Public repo updated to 1.22.0.
```

### DB reset / redeploy

```bash
docker compose down -v && docker compose up -d
```

### Adding profiles on the fly

```bash
# Need demo databases now
docker compose --profile demo up -d

# Need everything
docker compose --profile full up -d
```

## Profiles

| Profile | Services added | Use case |
|---------|---------------|----------|
| (none) | postgres, rabbitmq, grok_pipe, datagrok, tools-dev | Basic package dev |
| `demo` | + world, test_db, northwind | Need demo databases |
| `scripting` | + jupyter_kernel_gateway | Need Python/R/Julia scripts |
| `full` | + grok_connect, grok_spawner, demo DBs, JKG | Everything |

## JS API and code reference

The agent reads JS API source and documentation directly from the workspace — no
generated docs needed.

### Key paths (public repo at `/workspace`)

| What | Path |
|------|------|
| JS API source (types, classes, methods) | `js-api/src/` |
| JS API entry points (grok, ui, dg) | `js-api/grok.ts`, `js-api/ui.ts`, `js-api/dg.ts` |
| JS API CLAUDE.md (module map, patterns) | `js-api/CLAUDE.md` |
| Runnable code samples | `packages/ApiSamples/scripts/` |
| Package development guide | `help/develop/packages/` |
| All platform help docs | `help/` |
| Existing packages (patterns) | `packages/` |
| CLI tool source | `tools/` |

### Key paths (separate repo — public auto-cloned)

Same as above, prefixed with `public/`:
- `public/js-api/src/`, `public/packages/ApiSamples/scripts/`, etc.

The public branch matches `DG_VERSION` by default, so js-api types stay in sync
with the running Datagrok server.

### CLAUDE.md for the agent

Add to your project's CLAUDE.md so the agent knows where to look:

```markdown
## Reference

- JS API: read source in `js-api/src/` — see `js-api/CLAUDE.md` for module map
- Code samples: `packages/ApiSamples/scripts/` — runnable examples for eval
- Help docs: `help/develop/` — guides for packages, viewers, functions
- Existing packages: `packages/` — real-world patterns
```

## Exposing the client

- Datagrok UI: `http://localhost:${DG_PORT:-8080}`
- Default credentials: **admin / admin** (created on first deploy)

### grok CLI config inside the container

The entrypoint auto-creates `~/.grok/config.yaml` with the local Datagrok instance
(dev key `admin`). The grok CLI is ready to use immediately — no manual config needed.

If you need to reconfigure:

```bash
# Inside tools-dev:
grok config add --alias local \
  --server http://datagrok:8080/api \
  --key admin --default
```

## Claude Code inside the container

Mount `~/.claude/` from host (already in the compose file) for credentials and
config. Pass `ANTHROPIC_API_KEY` via `.env` or shell export.

```bash
# Single entry point — launch Claude Code interactively
docker exec -it <tools-dev_container> claude --dangerously-skip-permissions

# One-shot command
docker exec <tools-dev_container> claude -p "publish the Chem package" \
  --dangerously-skip-permissions
```

The agent can:
- Build, test, and publish packages
- Search Jira/GitHub for context (via MCP plugins)
- Manage the Datagrok stand (version switch, redeploy, DB reset)
- Ask the user for tokens and auth when needed
- Read js-api source and ApiSamples for reference

## Quick reference

```bash
# Start (basic)
docker compose up -d

# Start (everything)
docker compose --profile full up -d

# Launch Claude agent (single entry point)
docker exec -it <tools-dev> claude --dangerously-skip-permissions

# Shell into tools-dev
docker exec -it <tools-dev> bash

# Publish a package
docker exec <tools-dev> bash -c "cd /workspace/packages/MyPkg && grok publish"

# Run grok tests
docker exec <tools-dev> bash -c "cd /workspace/packages/MyPkg && grok test --host local"

# Run Playwright tests
docker exec <tools-dev> bash -c "cd /workspace/packages/MyPkg && npx playwright test"

# Version switch
DG_VERSION=1.22.0 docker compose up -d --pull always

# DB reset
docker compose down -v && docker compose up -d

# Logs
docker compose logs -f datagrok

# Stop
docker compose down
```

## Troubleshooting

### Datagrok not starting
```bash
docker compose logs datagrok
# Wait for DB migrations — first start takes 1-2 minutes
# Look for: "Server started on port 8080"
```

### grok publish fails
Ensure `~/.grok/config.yaml` inside the container points to `http://datagrok:8080/api`
with a valid dev key. The container reaches Datagrok via Docker DNS, not `localhost`.

### Chrome / Playwright not working
```bash
docker exec <tools-dev> google-chrome-stable --version
docker exec <tools-dev> npx playwright --version
```
If missing, pull the latest image: `docker pull datagrok/tools-dev:latest`

### MCP plugins not connecting
```bash
# Verify env vars are set
docker exec <tools-dev> env | grep -E 'JIRA|GITHUB'
# Re-register if needed
docker exec -it <tools-dev> claude mcp add mcp-atlassian -s user -- \
  npx -y mcp-atlassian --jira-url "$JIRA_URL" \
  --jira-username "$JIRA_USERNAME" --jira-token "$JIRA_TOKEN"
```

### Agent can't manage Docker (version switch, redeploy)
The Docker socket must be mounted and the `dev` user must be in the `docker` group:
```bash
docker exec <tools-dev> docker ps
```
If permission denied, check that `/var/run/docker.sock` is mounted and accessible.

### Container can't resolve `datagrok` hostname
All services must be on the same Docker network:
```bash
docker network inspect dg-pkg-${TASK_KEY:-default}-net
```

---
name: deploy-grokky
description: Deploy Grokky — start Docker, ensure grok_spawner, build images, publish package
disable-model-invocation: true
context: fork
---

# Deploy Grokky

This skill performs a full deployment of the Grokky package:
- Ensures Docker is running
- Ensures grok_spawner is running
- Builds required Docker images
- Starts runtime containers
- Builds and publishes the package

## Arguments

`$ARGUMENTS` may include:
- `key: <ANTHROPIC_API_KEY>`
- `host: <publish host>`

### Resolve arguments

1. If `host` is missing:
   Ask: "Which host should I publish to? (e.g., `localhost`, `dev`, `public`)"
   Wait for user input before continuing.

2. If `key` is missing:
   Run:
   ```bash
   echo $ANTHROPIC_API_KEY
   ```
   If empty, ask the user to provide the key and wait.

## Resolve paths

Set variables explicitly:

```bash
GROKKY_ROOT=$(cd "${CLAUDE_SKILL_DIR}/../../.." && pwd)
DOCKER_DIR=$(cd "$GROKKY_ROOT/../../docker" && pwd)
COMPOSE_FILE="$DOCKER_DIR/localhost.bleeding-edge.docker-compose.yaml"
```

---

## Step 1 — Ensure Docker is running

Run:

```bash
docker info > /dev/null 2>&1
```

If it fails:

```bash
open -a Docker
```

Then retry every 3 seconds (max 10 attempts):

```bash
for i in {1..10}; do
  sleep 3
  docker info > /dev/null 2>&1 && break
done
```

If still failing after 10 attempts:
Stop execution and tell user: "Docker is not running. Please start Docker manually."

---

## Step 2 — Ensure grok_spawner is running and fresh

Check if running:

```bash
docker ps --filter "name=grok_spawner" --format "{{.Names}}"
```

If no output, start it:

```bash
cd "$DOCKER_DIR"
GROK_SPAWNER_CORE_MODE=true docker compose -f "$COMPOSE_FILE" up -d grok_spawner
```

If running, check its age:

```bash
docker inspect --format '{{.Created}}' grok_spawner
```

Parse the timestamp and compare to now. If the container is older than 2 days, recreate it:

```bash
cd "$DOCKER_DIR"
GROK_SPAWNER_CORE_MODE=true docker compose -f "$COMPOSE_FILE" up -d --force-recreate grok_spawner
```

Do NOT start any other services.

---

## Step 3 — Detect the Docker network

The network name varies by setup. Detect it from the running grok_spawner container:

```bash
DOCKER_NETWORK=$(docker ps --filter "name=grok_spawner" --format "{{.Networks}}" | head -1)
```

If `DOCKER_NETWORK` is empty, fall back:

```bash
DOCKER_NETWORK=$(docker network ls --format "{{.Name}}" | grep -E "datagrok" | head -1)
```

If still empty, stop and tell user: "No datagrok Docker network found. Is the platform running?"

---

## Step 4 — Detect required image builds

```bash
cd "$GROKKY_ROOT"
CHANGED_FILES=$(git diff --name-only HEAD -- dockerfiles/; git status --short -- dockerfiles/)
```

Check existing images:

```bash
EXISTING_IMAGES=$(docker images --format "{{.Repository}}:{{.Tag}}" | grep grokky || true)
```

Determine build flags:

- Build `claude-runtime` if:
  - `dockerfiles/claude-runtime/` appears in `$CHANGED_FILES` OR
  - `grokky-claude-runtime:admin` not found in `$EXISTING_IMAGES`

- Build `mcp-server` if:
  - `dockerfiles/mcp-server/` appears in `$CHANGED_FILES` OR
  - `grokky-mcp-server:admin` not found in `$EXISTING_IMAGES`

If neither needs building, skip to Step 6.

---

## Step 5 — Build images (only those flagged in Step 4)

### claude-runtime

```bash
cd "$GROKKY_ROOT/dockerfiles/claude-runtime"
docker build --build-arg ANTHROPIC_API_KEY=$KEY -t grokky-claude-runtime:admin .
docker rm -f grokky-claude-runtime 2>/dev/null || true
```

### mcp-server

```bash
cd "$GROKKY_ROOT/dockerfiles/mcp-server"
docker build -t grokky-mcp-server:admin .
docker rm -f grokky-mcp-server 2>/dev/null || true
```

---

## Step 6 — Build and publish

```bash
cd "$GROKKY_ROOT"
npm run build
```

If build fails, stop and show the error output.

If build succeeds:

```bash
grok publish $HOST
```

If publish fails, stop and show the error output.

---

## Step 7 — Report results

Print:
- Whether Docker was started
- Whether grok_spawner was started
- Which images were built (or "all cached")
- Container statuses
- Publish result

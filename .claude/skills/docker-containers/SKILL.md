---
name: docker-containers
description: Ship a Docker container with a Datagrok package and call it over HTTP/WebSocket via the platform proxy
---

# docker-containers

## When to use

Your package needs to run code that doesn't fit in the browser or the
package's JS — a Python/CUDA model, a long-running service, a side-car
HTTP API. Triggers: "ship a Dockerfile in my package", "call my Python
container from the plugin", "stream a WebSocket from a container",
"make Datagrok build and run my image."

## Prerequisites

- A package scaffold (`grok create <Name>`). Paths below are relative to
  the package root.
- `datagrok-api` available — entry point is
  `grok.dapi.docker.dockerContainers` (knowledge `DG-FACT-129`). The
  sibling `grok.dapi.dockerfiles.*` form named in the article prose
  does NOT exist (drift `DG-FACT-DRIFT-049`).
- A target Datagrok instance where `grok-spawner` is running (required
  for image build + container lifecycle).
- Familiarity with `Dockerfile` basics (one `EXPOSE <port>` per image is
  the hard limit — knowledge `DG-FACT-136`).

## Steps

1. **Add the Dockerfile under `dockerfiles/`.**
   Single Dockerfile → place it directly in `dockerfiles/Dockerfile`.
   Multiple containers → one subfolder per image, each with its own
   `Dockerfile` (knowledge `DG-FACT-134`). Expose exactly one port — the
   platform proxy routes to that single port (knowledge `DG-FACT-136`).
   ```dockerfile
   FROM python:3.11-slim
   WORKDIR /app
   COPY . /app
   RUN pip install --no-cache-dir -r requirements.txt
   EXPOSE 5000
   CMD ["python", "app.py"]
   ```
   Expected: `<package>/dockerfiles/Dockerfile` exists. Real reference:
   `packages/Admetica/dockerfiles/Dockerfile`.

2. **Add `container.json` next to the Dockerfile (optional).**
   Same directory as the `Dockerfile`. All fields optional — defaults
   apply per knowledge `DG-FACT-132` (`cpu` 0.25, `gpu` 0, `memory` 512
   MB, `on_demand` false, `shutdown_timeout` null, `storage` 21 GB,
   `shm_size` 64 MB). Use `environmentVars` (production-canonical), NOT
   `env` (drift `DG-FACT-DRIFT-050`). Use `shutdown_timeout`, NOT
   `timeout_minutes` (drift `DG-FACT-DRIFT-051`).
   ```json
   {
     "cpu": 1.5,
     "memory": 2048,
     "on_demand": true,
     "shutdown_timeout": 60,
     "environmentVars": {
       "DATA_CONNECTION": "#{x.MyPackage:MyConnection}",
       "LOG_LEVEL": "info"
     }
   }
   ```
   Expected: file at `dockerfiles/container.json`. The
   `#{x.<Package>:<Entity>}` token is JSON-serialized at start time;
   only same-package entities resolve (knowledge `DG-FACT-133`).

3. **Publish the package.**
   The first publish queues an image build; subsequent publishes
   re-queue only if the build context changed.
   ```bash
   webpack
   grok publish dev
   ```
   Expected: exit code `0`. Image and container appear in
   **Platform → Dockers** (top: containers, bottom: images). Status
   indicators: green dot = running/ready, red dot = error, blinking
   grey dot = pending (starting/stopping/rebuilding), no dot on a
   container card = stopped (knowledge `DG-FACT-135`).

4. **Resolve the container by friendly name.**
   Friendly name is `<PackageName>` for a single Dockerfile, or
   `<PackageName>-<FolderName>` for multi-container packages
   (knowledge `DG-FACT-134`).
   ```typescript
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';

   const container: DG.DockerContainer = await grok.dapi.docker
     .dockerContainers.filter('mypackage').first();
   ```
   Expected: `container.id` is a non-empty GUID. `.list()`, `.run(id)`,
   `.stop(id)`, `.getContainerLogs(id, limit?)` are siblings on the
   same DataSource (knowledge `DG-FACT-129`).

5. **Start the container if it isn't already (`on_demand: true`).**
   Production check: only call `run` when status is neither `started`
   nor `checking` (knowledge `DG-FACT-135`). Pass `awaitStart=true` so
   the call doesn't return until the container is up.
   ```typescript
   if (container.status !== 'started' && container.status !== 'checking')
     await grok.dapi.docker.dockerContainers.run(container.id, true);
   ```
   Expected: subsequent `fetchProxy` returns `200`, not `400`.

6. **Send an HTTP request through the proxy.**
   Method shape mirrors the Fetch API; `params` is `RequestInit`
   (knowledge `DG-FACT-130`). Leading `/` on `path` is added if missing;
   `params.method` defaults to `GET`. Errors come back with HTTP 400
   (bad container status), 404 (container/path not found), or 500
   (server-side workflow error) and a `datagrok-error` JSON field.
   ```typescript
   const resp: Response = await grok.dapi.docker.dockerContainers
     .fetchProxy(container.id, '/predict', {
       method: 'POST',
       headers: {'Accept': 'application/json',
                 'Content-Type': 'application/json'},
       body: JSON.stringify({input: 'CCO'}),
     });
   if (!resp.ok) {
     const err = (await resp.json())['datagrok-error'] ?? resp.statusText;
     throw new Error(err);
   }
   const data = await resp.json();
   ```
   Expected: `resp.ok === true`; `data` is parsed JSON.

7. **Open a WebSocket through the proxy (optional).**
   `webSocketProxy(id, path, timeout=60000)` resolves only after the
   server emits the `"CONNECTED"` handshake message (knowledge
   `DG-FACT-131`). For `on_demand` containers, raise `timeout` past 60s
   to absorb cold-start.
   ```typescript
   const ws: WebSocket = await grok.dapi.docker.dockerContainers
     .webSocketProxy(container.id, '/ws', 120_000);
   ws.addEventListener('message', (e) => console.log(e.data));
   ws.send('hello');
   setTimeout(() => ws.close(), 3000);
   ```
   Expected: `ws.readyState === WebSocket.OPEN` immediately after the
   await. On timeout the socket closes with code `4001`.

## Common failure modes

- **`grok.dapi.dockerfiles.fetchProxy is not a function`.** Article
  prose at line 111 names a path that doesn't exist (drift
  `DG-FACT-DRIFT-049`). Fix: use
  `grok.dapi.docker.dockerContainers.fetchProxy` (knowledge
  `DG-FACT-129`/`130`).
- **Env vars never reach the container.** Article shows `"env": {...}`,
  but server-side model and production packages use
  `"environmentVars": {...}` (drift `DG-FACT-DRIFT-050`). Fix: rename
  the key.
- **Idle-shutdown ignored.** `"timeout_minutes"` is silently dropped —
  not in the model (drift `DG-FACT-DRIFT-051`). Fix: rename to
  `shutdown_timeout` (knowledge `DG-FACT-132`).
- **`fetchProxy` returns HTTP 400.** Container is not in `started`
  status (knowledge `DG-FACT-135`). Fix: gate on status and call
  `dockerContainers.run(id, true)` (step 5).
- **`fetchProxy` returns HTTP 404.** Wrong friendly name or unmapped
  path inside the container. Fix: confirm name per knowledge
  `DG-FACT-134`; inspect `getContainerLogs(id)`.
- **Image build queued forever / red dot.** Multiple `EXPOSE`
  directives (only one allowed — knowledge `DG-FACT-136`) or a build
  error. Fix: image card → Property pane → **Build logs**.
- **`#{x.Other:Conn}` resolves empty.** Cross-package entity refs are
  not injected (knowledge `DG-FACT-133`). Fix: move the connection
  into the same package.
- **`webSocketProxy` rejects with code 4001.** No `"CONNECTED"`
  handshake within `timeout` (knowledge `DG-FACT-131`). Fix: raise
  `timeout` for `on_demand: true` cold starts.

## Verification

- After step 3, **Platform → Dockers** shows a card for the image and
  the container with a green or grey dot (not red).
- After step 4, `container.id` is a GUID and `container.name` ends in
  the friendly name expected per `DG-FACT-134`.
- After step 6, `resp.status === 200` and `await resp.json()` parses.
- After step 7 (if used), `ws.readyState === WebSocket.OPEN`.

## See also

- Source articles:
  - `help/develop/how-to/packages/docker-containers.md`
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts
    `DG-FACT-129..136` and drifts `DG-FACT-DRIFT-049..051`.
- Related skills:
  - `db-in-docker` (containers that need a Postgres connection at
    start time — overlaps step 2's `environmentVars`).
  - `access-data` (use `grok.dapi.fetchProxy` for non-container REST
    calls; this skill is the in-package container variant).
  - `python-functions` (lighter weight: run Python without shipping a
    full image).

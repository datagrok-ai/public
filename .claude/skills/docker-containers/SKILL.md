---
name: docker-containers
version: 0.1.0
description: |
  Ship a Dockerized service alongside a Datagrok package so the platform
  builds the image, manages the container lifecycle, and exposes it to
  the package's JS via an authenticated HTTP/WebSocket proxy. For plugin
  authors whose pipeline needs Python, CUDA, a model server, or any
  long-running side-car that doesn't fit in the browser. Produces a
  `dockerfiles/` directory, an optional `container.json` config, and the
  JS glue that resolves the container and calls it.
  Use when asked to "call a Python model server from my plugin",
  "ship a GPU image with my package", or "stream WebSocket results from
  a side-car my package owns".
triggers:
  - call a python service from my plugin
  - ship a gpu image with my package
  - stream websocket results from a side-car
  - run my own model server inside datagrok
  - have the platform manage my image lifecycle
  - proxy http calls from js to my container
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
harness-authored: true
---

# docker-containers

## When to use

Your package needs code that doesn't fit in the browser or in package
JS — a Python/CUDA model, a long-running service, a side-car HTTP API —
and you want Datagrok to build the image, schedule the container, and
route authenticated traffic to it.

## Prerequisites

- A package scaffold (e.g. `grok create <Name>`); paths are relative to
  the package root.
- `datagrok-api` available — entry point `grok.dapi.docker.dockerContainers`
  (`DG-FACT-129`). The `grok.dapi.dockerfiles.*` path named in article
  prose does NOT exist (`DG-FACT-DRIFT-049`).
- Target instance running `grok-spawner` (article line 19); exactly one
  exposed port per image (`DG-FACT-136`).

## Steps

1. **Add the Dockerfile under `dockerfiles/`.**
   Single container → `dockerfiles/Dockerfile`. Multiple containers →
   one subfolder per image, each with its own `Dockerfile`
   (`DG-FACT-134`). Exactly one `EXPOSE <port>` — the platform proxy
   routes to that single port.
   ```dockerfile
   FROM python:3.11-slim
   WORKDIR /app
   COPY . /app
   RUN pip install --no-cache-dir -r requirements.txt
   EXPOSE 5000
   CMD ["python", "app.py"]
   ```
   Expected: `<package>/dockerfiles/Dockerfile` exists. Reference:
   `packages/Admetica/dockerfiles/Dockerfile`.

2. **Add `container.json` next to the Dockerfile (optional).**
   Same directory as the `Dockerfile`. Defaults per `DG-FACT-132`
   (`cpu` 0.25, `gpu` 0, `memory` 512 MB, `on_demand` false,
   `shutdown_timeout` null, `storage` 21 GB, `shm_size` 64 MB). Use
   `environmentVars` not `env` (`DG-FACT-DRIFT-050`); `shutdown_timeout`
   not `timeout_minutes` (`DG-FACT-DRIFT-051`).
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
   `#{x.<Package>:<Entity>}` token is JSON-serialized at container-start
   time; only same-package entities resolve — cross-package credentials
   are NOT injected (`DG-FACT-133`).

3. **Publish the package.**
   First publish queues an image build; later publishes re-queue only
   when the build context changes. Image and container appear in
   **Platform → Dockers**.
   ```bash
   webpack
   grok publish dev
   ```
   Expected: exit code `0`. Status dots per `DG-FACT-135` — green =
   running, red = error, blinking grey = pending, no dot = stopped.

4. **Resolve the container by friendly name.**
   Friendly name is `<PackageName>` for a single Dockerfile, or
   `<PackageName>-<FolderName>` for multi-container packages
   (`DG-FACT-134`).
   ```typescript
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';

   const container: DG.DockerContainer = await grok.dapi.docker
     .dockerContainers.filter('mypackage').first();
   ```
   Expected: `container.id` is a non-empty GUID. Siblings on the same
   DataSource: `.list()`, `.run(id, awaitStart?)`, `.stop(id)`,
   `.getContainerLogs(id, limit?)` (`DG-FACT-129`).

5. **Start the container if it isn't already (`on_demand: true`).**
   Call `run` only when status is neither `started` nor `checking`
   (`DG-FACT-135`); pass `awaitStart=true` to block until up.
   ```typescript
   if (container.status !== 'started' && container.status !== 'checking')
     await grok.dapi.docker.dockerContainers.run(container.id, true);
   ```
   Expected: subsequent `fetchProxy` returns `200`, not HTTP 400.

6. **Send an HTTP request through the proxy.**
   Method shape mirrors the Fetch API; `params` is `RequestInit`
   (`DG-FACT-130`). Live reference:
   `packages/MolTrack/src/services/moltrack-docker-service.ts:28`.
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
   Expected: `resp.ok === true`; `data` parses as JSON. Errors: HTTP
   400 (bad status), 404 (not found), 500 (workflow) + `datagrok-error`.

7. **Open a WebSocket through the proxy (optional).**
   `webSocketProxy(id, path, timeout=60000)` resolves only after the
   server emits the `"CONNECTED"` handshake (`DG-FACT-131`). For
   `on_demand` containers, raise `timeout` past 60s to absorb cold-start.
   ```typescript
   const ws: WebSocket = await grok.dapi.docker.dockerContainers
     .webSocketProxy(container.id, '/ws', 120_000);
   ws.addEventListener('message', (e) => console.log(e.data));
   ws.send('hello');
   setTimeout(() => ws.close(), 3000);
   ```
   Expected: `ws.readyState === WebSocket.OPEN` after the await. On
   timeout the socket closes with code `4001`.

## Common failure modes

- **`grok.dapi.dockerfiles.fetchProxy is not a function`** — article
  line 111 names a path that doesn't exist (`DG-FACT-DRIFT-049`).
  Fix: use `grok.dapi.docker.dockerContainers.fetchProxy`.
- **Env vars never reach the container** — article shows `"env": {...}`,
  but the server-side model and production packages (MolTrack,
  NodeJSDemo) use `"environmentVars"` (`DG-FACT-DRIFT-050`).
- **Idle-shutdown ignored** — `"timeout_minutes"` is silently dropped
  (`DG-FACT-DRIFT-051`); rename to `shutdown_timeout`.
- **`fetchProxy` returns HTTP 400** — container is not `started`
  (`DG-FACT-135`). Fix: gate on status and call
  `dockerContainers.run(id, true)` (step 5).
- **`fetchProxy` returns HTTP 404** — wrong friendly name or unmapped
  internal path. Fix: confirm name per `DG-FACT-134`; inspect
  `getContainerLogs(id)`.
- **Image build queued forever / red dot** — multiple `EXPOSE`
  directives (only one allowed — `DG-FACT-136`) or a build error.
  Fix: image card → Property pane → **Build logs**.
- **`#{x.Other:Conn}` resolves empty** — cross-package entity refs are
  not injected (`DG-FACT-133`); move the connection into the same package.
- **`webSocketProxy` rejects with code 4001** — no `"CONNECTED"`
  handshake within `timeout` (`DG-FACT-131`); raise `timeout` for
  `on_demand` cold starts.

## Verification

- After step 3, **Platform → Dockers** shows image + container cards
  with a green or blinking grey dot (not red).
- After step 4, `container.id` is a GUID matching `DG-FACT-134`'s
  friendly-name rule.
- After step 6, `resp.status === 200` and `await resp.json()` parses.
- After step 7 (if used), `ws.readyState === WebSocket.OPEN`.

## See also

- Source: `help/develop/how-to/packages/docker-containers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-129..136`; `contract-triangles.md` — DRIFT-049/050/051.
- Related skills: `db-in-docker`, `access-data`, `python-functions`.

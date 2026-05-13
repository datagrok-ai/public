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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# docker-containers

## When to use

Your package needs code that doesn't fit in the browser or in package
JS — a Python/CUDA model, a long-running service, a side-car HTTP API —
and you want Datagrok to build the image, schedule the container, and
route authenticated traffic to it.

## Prerequisites

- Target instance running `grok-spawner` (article line 19); exactly one
 exposed port per image (`DG-FACT-136`).

## Steps

1. **Add the Dockerfile under `dockerfiles/`.**
 Single container → `dockerfiles/Dockerfile`; multiple → one
 subfolder per image (`DG-FACT-134`). Exactly one `EXPOSE <port>`
 (`DG-FACT-136`).
 ```dockerfile
 FROM python:3.11-slim
 WORKDIR /app
 COPY. /app
 RUN pip install --no-cache-dir -r requirements.txt
 EXPOSE 5000
 CMD ["python", "app.py"]
 ```

2. **Add `container.json` next to the Dockerfile (optional).**
 Field names and defaults per `DG-FACT-132` — use `environmentVars`
 (not `env`), `shutdown_timeout` (not `timeout_minutes`).
 `#{x.<Package>:<Entity>}` env-var tokens resolve at container-start
 time, same-package entities only (`DG-FACT-133`).
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

3. **Publish the package.** Image build queues on first publish;
 image+container appear in **Platform → Dockers**.
 ```bash
 webpack
 grok publish dev
 ```
 Status dots per `DG-FACT-135`.

4. **Resolve the container by friendly name** (`DG-FACT-134`).
 ```typescript
 import * as grok from 'datagrok-api/grok';
 import * as DG from 'datagrok-api/dg';

 const container: DG.DockerContainer = await grok.dapi.docker
.dockerContainers.filter('mypackage').first;
 ```
 Siblings on the DataSource: `.list`, `.run(id, awaitStart?)`,
 `.stop(id)`, `.getContainerLogs(id, limit?)` (`DG-FACT-129`).

5. **Start the container if `on_demand: true` and not already up**
 (`DG-FACT-135`).
 ```typescript
 if (container.status !== 'started' && container.status !== 'checking')
 await grok.dapi.docker.dockerContainers.run(container.id, true);
 ```

6. **Send an HTTP request through the proxy.** Mirrors Fetch API
 (`DG-FACT-130`).
 ```typescript
 const resp: Response = await grok.dapi.docker.dockerContainers
.fetchProxy(container.id, '/predict', {
 method: 'POST',
 headers: {'Accept': 'application/json',
 'Content-Type': 'application/json'},
 body: JSON.stringify({input: 'CCO'}),
 });
 if (!resp.ok) {
 const err = (await resp.json)['datagrok-error'] ?? resp.statusText;
 throw new Error(err);
 }
 const data = await resp.json;
 ```
 Errors: HTTP 400 (bad status), 404 (not found), 500 (workflow).

7. **Open a WebSocket through the proxy (optional)** —
 `webSocketProxy(id, path, timeout=60000)` resolves after `"CONNECTED"`
 handshake; raise timeout for `on_demand` cold-starts (`DG-FACT-131`).
 ```typescript
 const ws: WebSocket = await grok.dapi.docker.dockerContainers
.webSocketProxy(container.id, '/ws', 120_000);
 ws.addEventListener('message', (e) => console.log(e.data));
 ws.send('hello');
 setTimeout( => ws.close, 3000);
 ```
 On timeout the socket closes with code `4001`.

## Common failure modes

- `grok.dapi.dockerfiles.fetchProxy is not a function` — article path is wrong; use `grok.dapi.docker.dockerContainers.fetchProxy`.
- Env vars not reaching container — use `environmentVars` not `env` (`DG-FACT-132`).
- Idle-shutdown ignored — use `shutdown_timeout` not `timeout_minutes` (`DG-FACT-132`).
- `fetchProxy` HTTP 400 — container not `started`; gate on status + `run(id, true)` (`DG-FACT-135`).
- `fetchProxy` HTTP 404 — wrong friendly name (`DG-FACT-134`) or unmapped path; inspect `getContainerLogs(id)`.
- Image build stuck / red dot — multiple `EXPOSE` directives or build error (`DG-FACT-136`).
- `#{x.Other:Conn}` empty — cross-package refs aren't injected (`DG-FACT-133`).
- `webSocketProxy` rejects with 4001 — no `"CONNECTED"` handshake within timeout (`DG-FACT-131`).

## See also

- Source: `help/develop/how-to/packages/docker-containers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
 `DG-FACT-129..136`; `contract-triangles.md` — DRIFT-049/050/051.
- Related skills: `db-in-docker`, `access-data`, `python-functions`.

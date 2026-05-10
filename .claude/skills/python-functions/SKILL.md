---
name: python-functions
description: Ship Celery-backed Python functions inside a Datagrok package and call them from JS via grok.functions.call
---

# python-functions

## When to use

Your package needs to run Python server-side — an ML predictor, a
NumPy/pandas transform, a model too heavy for the browser — and you
want it invocable as a Datagrok function from JS. Triggers: "expose
Python from my package", "call a sklearn model from the plugin",
"register Python function for the JS API."

## Prerequisites

- A package scaffold (`grok create <Name>`). Paths are relative to the
  package root.
- A target Datagrok instance where `grok-spawner` runs (it injects
  `$CELERY_HOSTNAME` / `$TASK_QUEUE_NAME` per knowledge `DG-FACT-162`).
- Python 3.8+ and `requirements.in` (uv) OR `environment.yaml`
  (micromamba + uv) for deps (knowledge `DG-FACT-163`).
- Familiarity with the `docker-containers` skill — Python functions
  ship as a Datagrok docker container under the hood.

## Steps

1. **Scaffold the app folder under `dockerfiles/`.**
   Canonical layout is `dockerfiles/<app>/` — NOT `python/<app>/` as
   the article states (drift `DG-FACT-DRIFT-060`; reference
   `packages/Samples/dockerfiles/`). Each app folder owns `app.py`,
   a deps file, optional `container.json`, and a Dockerfile
   (knowledge `DG-FACT-157`).
   ```bash
   mkdir -p dockerfiles/<app>
   ```
   Expected: `<package>/dockerfiles/<app>/` exists. `<app>` is
   conventionally lowercased package shorthand (`samples`, `chem`).

2. **Write `app.py` with Celery wiring + annotated functions.**
   Celery-main is NOT auto-generated (drift `DG-FACT-DRIFT-061`).
   Hand-write `app = Celery(...)` plus per-function
   `@app.task(name=..., bind=True, base=DatagrokTask)` wrappers — header
   comments alone are not enough to dispatch (knowledge `DG-FACT-159`).
   Header grammar is the standard Datagrok script grammar (`#name:`,
   `#language: python`, `#input:`, `#output:`, `#meta.*`); per-input
   options `category`, `choices: [...]`, `range: 1-2` are honored
   (knowledge `DG-FACT-158`).
   ```python
   # dockerfiles/<app>/app.py
   from celery import Celery
   from datagrok_celery_task import DatagrokTask, Settings, get_logger

   settings = Settings()
   app = Celery(settings.celery_name, broker=settings.broker_url)
   logger = get_logger()

   #name: Add
   #language: python
   #input: int x
   #input: int y
   #output: int z
   @app.task(name='add', bind=True, base=DatagrokTask)
   def add(self, x: int, y: int) -> int:
       return x + y
   ```
   Expected: `app.py` declares module-level `app` and at least one
   `@app.task`-decorated function. Reference:
   `packages/Samples/dockerfiles/app.py` (`PyKNNTrain`, `PyKNNApply`,
   `PyKNNIsApplicable`).

3. **Declare deps — `environment.yaml` OR `requirements.in`.**
   `datagrok-celery-task` is REQUIRED at runtime (it provides
   `DatagrokTask`, `Settings`, `get_logger`); the article never mentions
   it (drift `DG-FACT-DRIFT-062`). Pin to canonical `0.1.8` alongside
   `celery` (knowledge `DG-FACT-160`). Use `environment.yaml` for
   conda-managed natives (RDKit, PyTorch); else `requirements.in` is
   lighter (knowledge `DG-FACT-163`).
   ```yaml
   # dockerfiles/<app>/environment.yaml
   name: myenv
   channels: [conda-forge, defaults]
   dependencies:
     - python=3.11
     - pip
     - numpy
     - pandas
     - pip:
       - datagrok-celery-task==0.1.8
       - celery
   ```
   Expected: `datagrok-celery-task==0.1.8` listed under `pip:`. With
   `requirements.in`, list it on its own line and validate locally with
   `uv pip install -r requirements.in`.

4. **Add `container.json` (optional).**
   Same shape as `dockerfiles/container.json`, but min/default `cpu`
   for Celery containers is **1**, not `0.25` (knowledge `DG-FACT-161`).
   Article says worker count = `cpu * 2`; production Dockerfile
   hard-codes `--concurrency=1` (drift `DG-FACT-DRIFT-064`) — leave
   concurrency unset rather than rely on either claim.
   ```json
   {"cpu": 1, "memory": 1024}
   ```
   Expected: `dockerfiles/<app>/container.json` exists. See
   `docker-containers` skill for the full field list.

5. **Ship the Dockerfile.**
   Article claims one is auto-generated (drift `DG-FACT-DRIFT-063`);
   production hand-writes it and must launch a Celery worker bound to
   the platform-injected queue. `$CELERY_HOSTNAME` / `$TASK_QUEUE_NAME`
   are set by `grok-spawner` at start time (knowledge `DG-FACT-162`).
   ```dockerfile
   FROM mambaorg/micromamba:1.5.3
   USER root
   WORKDIR /app
   COPY environment.yaml .
   RUN micromamba create -n myenv -f environment.yaml && \
       micromamba clean --all --yes
   SHELL ["micromamba", "run", "-n", "myenv", "/bin/bash", "-c"]
   COPY . .
   EXPOSE 5555
   ENTRYPOINT micromamba run -n myenv celery -A app worker \
       --loglevel=info --hostname=$CELERY_HOSTNAME -Q $TASK_QUEUE_NAME
   ```
   Expected: `dockerfiles/<app>/Dockerfile` exists. `-A app` resolves to
   the module-level `app` symbol from step 2 — that's why the file is
   `app.py`.

6. **Publish the package.**
   First publish queues the image build; one worker container deploys
   per app folder.
   ```bash
   webpack
   grok publish dev
   ```
   Expected: exit code `0`. Image and container appear in
   **Platform → Dockers** with a green dot once the build settles.

7. **Call the function from JS.**
   Use the standard dispatcher with `<Package>:<HeaderName>`.
   `<HeaderName>` is the `#name:` header value — NOT the Python `def`
   name and NOT the Celery `@app.task(name=...)` identifier (knowledge
   `DG-FACT-164`; misleading article example in drift
   `DG-FACT-DRIFT-065`).
   ```typescript
   import * as grok from 'datagrok-api/grok';

   const z: number = await grok.functions.call(
     `${_package.name}:Add`, {x: 1, y: 2});
   ```
   Expected: `z === 3`. From the Datagrok console the equivalent is
   `<PackageName>:Add(1, 2)`.

## Common failure modes

- **`Plugin:add` from copy-paste returns "function not found".** Article
  example uses literal `Plugin` (drift `DG-FACT-DRIFT-065`). Fix: use
  the real package name + `#name:` header (e.g. `Samples:PyKNNTrain`).
- **`ImportError: No module named datagrok_celery_task`.** Required
  package not declared (drift `DG-FACT-DRIFT-062`). Fix: add
  `datagrok-celery-task==0.1.8` per step 3 (knowledge `DG-FACT-160`).
- **Function visible to JS but never executes.** Header comments
  present, but `@app.task(... base=DatagrokTask)` wrapper missing
  (drift `DG-FACT-DRIFT-061`). Fix: wrap each annotated `def` per
  step 2 (knowledge `DG-FACT-159`).
- **Worker starts but no tasks land in it.** Dockerfile omits
  `-Q $TASK_QUEUE_NAME` or hard-codes a queue name. Fix: keep the env
  var unmodified (knowledge `DG-FACT-162`).
- **`container.json` `cpu: 0.25` rejected at deploy.** Celery
  containers require `cpu >= 1` (knowledge `DG-FACT-161`). Fix: bump
  to `1`+.
- **Wrote files under `python/<app>/`, platform never built them.**
  Article documents that path but no production package follows it
  (drift `DG-FACT-DRIFT-060`). Fix: move to `dockerfiles/<app>/`.

## Verification

- After step 6, **Platform → Dockers** shows a card for the app's image
  and container with a green or grey dot (not red); the function is
  visible in **Functions** as `<Package>:<HeaderName>`.
- After step 7, `await grok.functions.call(...)` resolves; check
  container logs from the Property pane on failure.

## See also

- Source articles:
  - `help/develop/how-to/packages/python-functions.md`
  - `help/develop/how-to/packages/docker-containers.md` (sibling for
    `container.json` field reference)
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts
    `DG-FACT-157..164` and drifts `DG-FACT-DRIFT-060..065`.
- Related skills:
  - `docker-containers` (lifecycle + `fetchProxy` for non-Celery
    HTTP-server containers; this skill is the function-dispatch variant).
  - `publish-packages` (how `grok publish dev` queues image builds).
  - `access-data` (post-process query results in Python before
    returning to JS).

---
name: python-functions
version: 0.1.0
description: |
  Ship a server-side Python compute (numpy / pandas / scikit-learn /
  PyTorch / RDKit) inside a Datagrok plugin so users can invoke it
  from JS via `grok.functions.call`. For plugin authors with heavy
  analysis that doesn't belong in the browser. Produces a
  `dockerfiles/<app>/` folder with `app.py`, deps, optional
  `container.json`, and a Dockerfile that runs a Celery worker bound
  to the platform's task queue.
  Use when asked to "run a numpy analysis server-side from my plugin",
  "call a trained sklearn model from JS", or "expose ML inference alongside the package".
triggers:
  - run a numpy job server-side
  - call a trained ml model from js
  - ship sklearn inference with the plugin
  - queue heavy compute from the ui
  - expose pandas analysis to a button
  - run rdkit alongside the package
allowed-tools:
  - Read
  - Write
  - Bash
harness-authored: true
---

# python-functions

## When to use

The plugin needs to run a server-side compute — ML predictor, NumPy /
pandas transform, RDKit / PyTorch job — invocable as a regular
Datagrok function from JS.

## Prerequisites

- A package scaffold (`grok create <Name>`); paths below are relative
  to the package root.
- A Datagrok instance where `grok_spawner` runs (it injects
  `$CELERY_HOSTNAME` / `$TASK_QUEUE_NAME` — knowledge `DG-FACT-162`).
- Python 3.8+ with `requirements.in` (uv) OR `environment.yaml`
  (micromamba + uv) for deps (knowledge `DG-FACT-163`).
- Familiarity with the `docker-containers` skill — Python functions
  ship as a Datagrok docker container.

## Steps

1. **Scaffold the app folder under `dockerfiles/`.**
   Each Python "app" is its own folder owning `app.py`, a deps file,
   optional `container.json`, and a Dockerfile (knowledge
   `DG-FACT-157`). Reference: `packages/Samples/dockerfiles/`.
   ```bash
   mkdir -p dockerfiles/<app>
   ```
   Expected: `<package>/dockerfiles/<app>/` exists. `<app>` is
   conventionally lowercased package shorthand (`samples`, `chem`).

2. **Write `app.py` with Celery wiring + annotated functions.**
   Celery-main is NOT auto-generated despite the article's claim
   (drift `DG-FACT-DRIFT-061`). Hand-write `app = Celery(...)` plus
   per-function `@app.task(name=..., bind=True, base=DatagrokTask)`
   wrappers — the platform parses Datagrok header comments to register
   the function but does NOT synthesize Celery glue (knowledge
   `DG-FACT-159`). Header grammar: `#name:`, `#language: python`,
   `#input:` (with `category`, `choices`, `range`), `#output:`,
   `#meta.*` (knowledge `DG-FACT-158`).
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
   `@app.task`-decorated function (ref `packages/Samples/dockerfiles/app.py`).

3. **Declare deps — `environment.yaml` OR `requirements.in`.**
   `datagrok-celery-task` is REQUIRED at runtime (supplies
   `DatagrokTask`, `Settings`, `get_logger`); article omits it
   (drift `DG-FACT-DRIFT-062`). Pin to `0.1.8` alongside `celery`
   (knowledge `DG-FACT-160`). Use `environment.yaml` for conda-managed
   natives (RDKit, PyTorch); else `requirements.in` is lighter
   (knowledge `DG-FACT-163`).
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
   `requirements.in`, list it directly and validate with
   `uv pip install -r requirements.in`.

4. **Add `container.json` (optional).**
   Same shape as `dockerfiles/container.json`, but min and default
   `cpu` for Celery containers is **1**, not `0.25` (knowledge
   `DG-FACT-161`). Article claims worker count = `cpu * 2`; the
   production Dockerfile hard-codes `--concurrency=1` (drift
   `DG-FACT-DRIFT-064`) — leave concurrency unset.
   ```json
   {"cpu": 1, "memory": 1024}
   ```
   Expected: `dockerfiles/<app>/container.json` exists. See
   `docker-containers` skill for the full field list.

5. **Ship the Dockerfile.**
   Article claims it's auto-generated (drift `DG-FACT-DRIFT-063`);
   production hand-writes it and must launch a Celery worker bound to
   the platform-injected queue. `$CELERY_HOSTNAME` /
   `$TASK_QUEUE_NAME` come from `grok_spawner` at start time
   (knowledge `DG-FACT-162`).
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
   Expected: `dockerfiles/<app>/Dockerfile` exists. `-A app` resolves
   to the module-level `app` symbol from step 2.

6. **Publish the package.** First publish queues the image build; one
   worker container deploys per app folder.
   ```bash
   webpack
   grok publish dev
   ```
   Expected: exit `0`. Image and container appear in **Platform →
   Dockers** with a green dot once the build settles.

7. **Call the function from JS.**
   Use the standard dispatcher with `<Package>:<HeaderName>`.
   `<HeaderName>` is the `#name:` header value — NOT the Python `def`
   name and NOT the Celery `@app.task(name=...)` identifier
   (knowledge `DG-FACT-164`).
   ```typescript
   import * as grok from 'datagrok-api/grok';

   const z: number = await grok.functions.call(
     `${_package.name}:Add`, {x: 1, y: 2});
   ```
   Expected: `z === 3`.

## Common failure modes

- **`<YourPackage>:add` returns "function not found".** Article shows
  a placeholder; substitute the real package name and `#name:` header
  (e.g. `Samples:PyKNNTrain`).
- **`ImportError: No module named datagrok_celery_task`.** Missing
  required dep (drift `DG-FACT-DRIFT-062`). Fix: add
  `datagrok-celery-task==0.1.8` (knowledge `DG-FACT-160`).
- **Function visible to JS but never executes.** Header comments
  present, `@app.task(... base=DatagrokTask)` wrapper missing (drift
  `DG-FACT-DRIFT-061`). Fix: wrap each annotated `def` per step 2.
- **Worker starts but no tasks land in it.** Dockerfile omits or
  hard-codes the queue. Fix: keep `-Q $TASK_QUEUE_NAME` unmodified
  (knowledge `DG-FACT-162`).
- **`cpu: 0.25` rejected at deploy.** Celery containers require
  `cpu >= 1` (knowledge `DG-FACT-161`). Fix: bump to `1`+.

## Verification

- After step 6, **Platform → Dockers** shows a card for the app's
  image and container with a green or grey dot (not red); the function
  appears in **Functions** as `<Package>:<HeaderName>`.
- After step 7, `grok.functions.call(...)` resolves; on failure, check
  container logs from the Property pane.

## See also

- Source: `help/develop/how-to/packages/python-functions.md`,
  `help/develop/how-to/packages/docker-containers.md` (sibling).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-157..164`, drifts `DG-FACT-DRIFT-061..064`.
- Related skills: `docker-containers`, `publish-packages`, `access-data`.

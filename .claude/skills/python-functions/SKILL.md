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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# python-functions

## When to use

The plugin needs to run a server-side compute — ML predictor, NumPy /
pandas transform, RDKit / PyTorch job — invocable as a regular
Datagrok function from JS.

## Prerequisites

- A Datagrok instance where `grok_spawner` runs (it injects
  `$CELERY_HOSTNAME` / `$TASK_QUEUE_NAME` — knowledge `DG-FACT-162`).
- Python 3.8+ with `requirements.in` (uv) OR `environment.yaml`
  (micromamba + uv) for deps (knowledge `DG-FACT-163`).
- Familiarity with the `docker-containers` skill — Python functions
  ship as a Datagrok docker container.

## Steps

1. **Scaffold the app folder under `dockerfiles/`** (`DG-FACT-157`).
   ```bash
   mkdir -p dockerfiles/<app>
   ```

2. **Write `app.py` with Celery wiring + annotated functions.**
   Hand-write the Celery boilerplate; the platform parses header
   comments but does NOT auto-generate Celery glue (`DG-FACT-159`,
   `DG-FACT-DRIFT-061`). Header grammar at `DG-FACT-158`.
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

3. **Declare deps — `environment.yaml` OR `requirements.in`.**
   `datagrok-celery-task==0.1.8` is REQUIRED (`DG-FACT-160`,
   `DG-FACT-DRIFT-062`). Use `environment.yaml` for conda-managed
   natives (RDKit, PyTorch); `requirements.in` for pure-pip
   (`DG-FACT-163`).
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

4. **Add `container.json` (optional).** Celery containers require
   `cpu >= 1` (`DG-FACT-161`); leave concurrency unset
   (`DG-FACT-DRIFT-064`).
   ```json
   {"cpu": 1, "memory": 1024}
   ```

5. **Ship the Dockerfile.** Hand-written (`DG-FACT-DRIFT-063`); must
   launch a Celery worker bound to the platform-injected queue
   (`DG-FACT-162` — `$CELERY_HOSTNAME`/`$TASK_QUEUE_NAME` from
   `grok_spawner`).
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

6. **Publish the package.** First publish queues the image build; one
   worker container deploys per app folder.
   ```bash
   webpack
   grok publish dev
   ```

7. **Call the function from JS** with `<Package>:<HeaderName>` —
   `<HeaderName>` is `#name:`, NOT the Python `def` or Celery task name
   (`DG-FACT-164`).
   ```typescript
   import * as grok from 'datagrok-api/grok';

   const z: number = await grok.functions.call(
     `${_package.name}:Add`, {x: 1, y: 2});
   ```

## Common failure modes

- `<YourPackage>:add` "function not found" — substitute real package + `#name:` header.
- `ImportError: No module named datagrok_celery_task` — add `datagrok-celery-task==0.1.8` (`DG-FACT-160`).
- Function visible but never executes — missing `@app.task(... base=DatagrokTask)` wrapper (`DG-FACT-DRIFT-061`).
- Worker starts but no tasks land — keep `-Q $TASK_QUEUE_NAME` unmodified (`DG-FACT-162`).
- `cpu: 0.25` rejected — Celery containers need `cpu >= 1` (`DG-FACT-161`).

## See also

- Source: `help/develop/how-to/packages/python-functions.md`,
  `help/develop/how-to/packages/docker-containers.md` (sibling).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-157..164`, drifts `DG-FACT-DRIFT-061..064`.
- Related skills: `docker-containers`, `publish-packages`, `access-data`.

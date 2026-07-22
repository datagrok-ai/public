---
title: "JS server functions (meta.queue)"
description: Run package JS/TS functions server-side as Celery tasks in Docker-based Node workers managed automatically by Datagrok.
keywords:
  - Celery task
  - meta.queue
  - Node worker
  - server function
  - RabbitMQ worker
  - docker container
---

## Overview

A function declared in your package's `src/package.ts` normally runs in the browser. Annotate it with
`meta.queue` and it becomes a **server function**: on publish, Datagrok registers it as a Docker
function, spins up a Docker-based Node worker container, and routes every call through the platform's
task queue (RabbitMQ + Celery protocol) — exactly like [Python Docker functions](python-functions.md),
but for JavaScript/TypeScript.

The worker loads your published package bundle via the js-api Node runtime, so the function body runs
with the familiar `grok` / `DG` APIs available (server-side, no UI).

## Declaring a function

Comment annotation:

```typescript
//name: heavyCompute
//meta.queue: true
//input: dataframe df
//output: dataframe result
export function heavyCompute(df: DG.DataFrame): DG.DataFrame {
  // runs inside the Node worker container, not the browser
  return df;
}
```

Or with a decorator:

```typescript
@grok.decorators.func({meta: {queue: true}})
```

Calling it is unchanged:

```typescript
const result = await grok.functions.call('MyPackage:heavyCompute', {df});
```

## Container options

* `//meta.queue: true` — the platform auto-creates a worker container from the stock
  `datagrok/celery_node_worker` image. Optionally add `queue/container.json` in the package root to
  configure resources: `{"cpu": 1, "memory": 1024, "gpu": 0, "on_demand": true, "shutdown_timeout": 60}`.
  The `queue/` folder is reserved for this purpose.
* `//meta.queue: <dirName>` — bind the function to your own container defined in
  `dockerfiles/<dirName>/`. Base your Dockerfile on the worker image:

  ```dockerfile
  FROM datagrok/celery_node_worker:bleeding-edge
  EXPOSE 8000
  ```

  A container cannot mix Python Celery tasks (`.py` files) and JS queue functions.

## Inside the function

* All package code and bundled npm dependencies are available (the worker loads the published webpack
  bundle). `grok.dapi.*` calls run with the calling user's session token.
* Report progress with the worker-provided global (a no-op in the browser):

  ```typescript
  (globalThis as any).DG_TASK_PROGRESS?.(50, 'half way');
  ```

  When the call was cancelled, this throws, letting long loops stop cooperatively.

## Limitations

* One output parameter per function (same as Python Celery tasks).
* DataFrames are transferred as CSV (no Parquet yet).
* Cancellation is cooperative: the call completes as `Canceled` immediately, but a running function
  body only observes it at the next `DG_TASK_PROGRESS` call.
* `meta.queue` is only recognized in `package.ts` (not detectors or test entrypoints).
* Package `init` functions run once per worker process, with the first caller's token.

See also:

* [Creating Python Docker Apps](python-functions.md)
* [Docker containers in packages](docker-containers.md)

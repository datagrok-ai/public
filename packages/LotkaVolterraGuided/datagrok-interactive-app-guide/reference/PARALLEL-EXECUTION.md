# Parallel Execution Guide

Reference for distributing independent computations across multiple web workers.

> **End-to-end example:** function `runOptimization()` in `../example/code/src/levins/app.ts` — worker pool with task queue, progress bar, partial error handling.

For single-worker patterns, see `WORKER-GUIDE.md`.

## Worker Count

```typescript
import {MIN_WORKERS_COUNT, WORKERS_COUNT_DOWNSHIFT} from './worker-utils/worker-defs';

const workerCount = Math.max(MIN_WORKERS_COUNT, navigator.hardwareConcurrency - WORKERS_COUNT_DOWNSHIFT);
```

## Fan-out / Fan-in Pattern

```typescript
async function runParallel<TInput, TOutput>(
  inputs: TInput[],
  workerUrl: URL,
): Promise<TOutput[]> {
  const nWorkers = Math.min(
    Math.max(MIN_WORKERS_COUNT, navigator.hardwareConcurrency - WORKERS_COUNT_DOWNSHIFT),
    inputs.length,
  );

  // Distribute inputs round-robin
  const chunks: TInput[][] = Array.from({length: nWorkers}, () => []);
  for (let i = 0; i < inputs.length; i++)
    chunks[i % nWorkers].push(inputs[i]);

  const promises = chunks.map((chunk) =>
    new Promise<TOutput[]>((resolve, reject) => {
      const worker = new Worker(workerUrl);
      worker.postMessage(chunk);
      worker.onmessage = (e) => {
        worker.terminate();
        if (e.data.success)
          resolve(e.data.data);
        else
          reject(new Error(e.data.error));
      };
      worker.onerror = (e) => {
        worker.terminate();
        reject(new Error(e.message));
      };
    }),
  );

  const results = await Promise.all(promises);
  return results.flat();
}
```

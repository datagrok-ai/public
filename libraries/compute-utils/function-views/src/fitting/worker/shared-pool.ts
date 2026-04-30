// Singleton fitting worker pool.
//
// The default `runWithEphemeralPool` flow constructs and disposes a fresh
// `WorkerPool` per fit, throwing away the per-worker compile cache and
// paying the worker spin-up + per-slot setup compile every time. For the
// common Compute2 workload (15 fits/click, several clicks per session,
// few unique script bodies) those costs amortize to zero from the second
// fit on if the pool persists. This module owns that long-lived pool.
//
// Lazy: spins workers on first request. Browser-tab lifetime — never
// disposed during normal operation; the tab's GC reclaims it on close.
// Tests must call `disposeSharedFittingPool()` to tear it down between
// scenarios.

import {WorkerPool, defaultPoolSize, WorkerPoolOptions} from './pool';

let sharedPool: WorkerPool | null = null;

export function getSharedFittingPool(opts?: WorkerPoolOptions): WorkerPool {
  if (sharedPool) return sharedPool;
  sharedPool = new WorkerPool(defaultPoolSize(), opts);
  return sharedPool;
}

export function disposeSharedFittingPool(): void {
  if (!sharedPool) return;
  sharedPool.dispose();
  sharedPool = null;
}

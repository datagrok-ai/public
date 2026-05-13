// Tab-lifetime singleton fitting pool: amortizes worker spin-up and per-slot
// compile cost across fits. Lazy on first request; tests must call
// `disposeSharedFittingPool()` between scenarios.

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

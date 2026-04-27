// Long-lived worker pool with a work-stealing seed queue.
//
// Lazy creation: workers spin up on the first dispatch. Each worker stays
// alive across tasks so the in-worker compile cache (func-call-shim.ts) is
// reused, paying the body's parse cost exactly once.
//
// Per-seed dispatch (not pre-batched batches): a slow seed on one worker no
// longer bottlenecks the whole fit because faster workers steal the next
// pending seed as soon as they finish.

import type {WorkerReply, WorkerTask} from './serialize';

export type PoolTask = WorkerTask;
export type PoolReply = WorkerReply;

interface Pending {
  task: PoolTask;
  resolve: (r: PoolReply) => void;
  transferables: Transferable[];
}

type Slot = {
  worker: Worker;
  busy: boolean;
  current: Pending | null;
};

export class WorkerPool {
  private slots: Slot[] = [];
  private queue: Pending[] = [];
  private disposed = false;
  private nextTaskId = 1;

  constructor(private readonly size: number) {}

  // Construct the underlying workers lazily, on first dispatch. Lets the
  // pool be allocated cheaply (e.g. per-FittingView) without paying the
  // worker spin-up cost until a fit actually runs.
  private ensureWorkers(): void {
    if (this.slots.length > 0 || this.disposed) return;
    for (let i = 0; i < this.size; ++i) {
      const w = new Worker(new URL('./fitting.worker.ts', import.meta.url));
      const slot: Slot = {worker: w, busy: false, current: null};
      w.onmessage = (ev: MessageEvent<PoolReply>) => this.onMessage(slot, ev.data);
      w.onerror = (ev: ErrorEvent) => this.onError(slot, ev);
      this.slots.push(slot);
    }
  }

  private onMessage(slot: Slot, reply: PoolReply): void {
    const pending = slot.current;
    slot.current = null;
    slot.busy = false;
    pending?.resolve(reply);
    this.pump();
  }

  private onError(slot: Slot, ev: ErrorEvent): void {
    const pending = slot.current;
    slot.current = null;
    slot.busy = false;
    if (pending) {
      pending.resolve({
        kind: 'failure',
        taskId: pending.task.taskId,
        message: `worker error: ${ev.message ?? 'unknown'}`,
        failKind: 'other',
        seed: pending.task.seed,
      });
    }
    this.pump();
  }

  private pump(): void {
    if (this.disposed) return;
    while (this.queue.length > 0) {
      const slot = this.slots.find((s) => !s.busy);
      if (!slot) return;
      const pending = this.queue.shift()!;
      slot.busy = true;
      slot.current = pending;
      slot.worker.postMessage(pending.task, pending.transferables);
    }
  }

  dispatch(args: {
    task: Omit<PoolTask, 'taskId'>;
    transferables: Transferable[];
  }): Promise<PoolReply> {
    if (this.disposed) return Promise.reject(new Error('pool disposed'));
    this.ensureWorkers();
    const taskId = this.nextTaskId++;
    const task = {...args.task, taskId} as PoolTask;
    return new Promise<PoolReply>((resolve) => {
      this.queue.push({task, resolve, transferables: args.transferables});
      this.pump();
    });
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    for (const slot of this.slots)
      try {slot.worker.terminate();} catch {/* ignore */}

    this.slots = [];
    // Drain pending — resolve as failures so callers don't hang.
    for (const pending of this.queue) {
      pending.resolve({
        kind: 'failure',
        taskId: pending.task.taskId,
        message: 'pool disposed before task ran',
        failKind: 'other',
        seed: pending.task.seed,
      });
    }
    this.queue = [];
  }
}

export function defaultPoolSize(): number {
  const hc = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 2;
  return Math.max(2, hc - 1);
}

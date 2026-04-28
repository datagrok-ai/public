// Worker pool with two-phase dispatch (setup + run-seed) plus session
// teardown.
//
// `setupAll(setup)` ships an immutable bundle to every slot once per fit;
// `dispatchRun(run)` ships only the seed and routes to a slot whose primed
// session set contains the run's sessionId; `dropSession(id)` releases the
// session's per-worker state at fit teardown.
//
// Lazy creation: workers spin up on the first `setupAll`. Each worker stays
// alive across tasks and across sessions, so the in-worker compile cache is
// reused.

import type {
  FitSessionSetup,
  RunSeed,
  DropSession,
  WorkerInbound,
  SetupAck,
  WorkerSuccess,
  WorkerFailure,
  SessionId,
} from './wire-types';

export type RunReply = WorkerSuccess | WorkerFailure;

interface PendingRun {
  run: RunSeed;
  resolve: (r: RunReply) => void;
  transferables: Transferable[];
}

interface PendingSetup {
  sessionId: SessionId;
  resolve: (ack: SetupAck) => void;
}

type Slot = {
  worker: Worker;
  busy: boolean;
  current: PendingRun | null;
  sessions: Set<SessionId>;
  // Per-slot pending setup ack indexed by sessionId — at most one per
  // session per slot in flight.
  pendingSetups: Map<SessionId, PendingSetup>;
};

/** Bounds on every promise returned from the pool, so that a worker that
 *  silently wedges (infinite loop in user JS, dead under structured-clone
 *  errors, etc.) can't hang a fit indefinitely. */
export interface WorkerPoolOptions {
  /** Reject setupAll's per-slot ack after this duration. Default: 20000. */
  setupTimeoutMs?: number;
  /** Reject a single dispatched seed (one Nelder-Mead run from one
   *  starting point) after this duration. Default: 60000 (1 min). */
  runTimeoutMs?: number;
}

export class WorkerPool {
  private slots: Slot[] = [];
  private queue: PendingRun[] = [];
  private disposed = false;
  private nextTaskId = 1;
  private readonly setupTimeoutMs: number;
  private readonly runTimeoutMs: number;

  constructor(private readonly size: number, opts: WorkerPoolOptions = {}) {
    this.setupTimeoutMs = opts.setupTimeoutMs ?? 20_000;
    this.runTimeoutMs = opts.runTimeoutMs ?? 60_000;
  }

  // Construct workers lazily, on first setup. Lets the pool be allocated
  // cheaply (e.g. per-FittingView) without paying the worker spin-up cost
  // until a fit actually runs.
  private ensureWorkers(): void {
    if (this.slots.length > 0 || this.disposed) return;
    for (let i = 0; i < this.size; ++i) {
      const w = new Worker(new URL('./fitting.worker.ts', import.meta.url));
      const slot: Slot = {
        worker: w,
        busy: false,
        current: null,
        sessions: new Set(),
        pendingSetups: new Map(),
      };
      w.onmessage = (ev: MessageEvent<WorkerInbound>) => this.onMessage(slot, ev.data);
      w.onerror = (ev: ErrorEvent) => this.onError(slot, ev);
      this.slots.push(slot);
    }
  }

  private onMessage(slot: Slot, msg: WorkerInbound): void {
    if (msg.kind === 'setup-ack') {
      const pending = slot.pendingSetups.get(msg.sessionId);
      if (!pending) return;
      slot.pendingSetups.delete(msg.sessionId);
      if (msg.ok) slot.sessions.add(msg.sessionId);
      pending.resolve(msg);
      return;
    }
    // success / failure for a run-seed
    const pending = slot.current;
    slot.current = null;
    slot.busy = false;
    pending?.resolve(msg);
    this.pump();
  }

  private onError(slot: Slot, ev: ErrorEvent): void {
    const pending = slot.current;
    slot.current = null;
    slot.busy = false;
    if (pending) {
      pending.resolve({
        kind: 'failure',
        taskId: pending.run.taskId,
        message: `worker error: ${ev.message ?? 'unknown'}`,
        failKind: 'other',
        seed: pending.run.seed,
      });
    }
    // Reject any in-flight setup acks for this slot too.
    for (const ps of slot.pendingSetups.values()) {
      ps.resolve({
        kind: 'setup-ack',
        sessionId: ps.sessionId,
        ok: false,
        message: `worker error: ${ev.message ?? 'unknown'}`,
      });
    }
    slot.pendingSetups.clear();
    this.pump();
  }

  private pump(): void {
    if (this.disposed) return;
    // Walk the queue and place each pending run on a slot that has the
    // session primed AND is idle. Items that can't be placed yet stay in
    // the queue; they'll be retried on the next pump (next reply).
    let i = 0;
    while (i < this.queue.length) {
      const pending = this.queue[i];
      const slot = this.slots.find(
        (s) => !s.busy && s.sessions.has(pending.run.sessionId));
      if (!slot) { ++i; continue; }
      this.queue.splice(i, 1);
      slot.busy = true;
      slot.current = pending;
      // Per-run timeout: terminates the slot if the worker wedges (e.g.
      // infinite loop in user JS body) so the dispatchRun promise resolves
      // instead of hanging forever.
      const timer = setTimeout(() => {
        if (slot.current !== pending) return;
        slot.current = null;
        slot.busy = false;
        this.removeSlot(slot, 'run timed out');
        pending.resolve({
          kind: 'failure', taskId: pending.run.taskId,
          message: `run timed out after ${this.runTimeoutMs}ms`,
          failKind: 'other', seed: pending.run.seed,
        });
      }, this.runTimeoutMs);
      const originalResolve = pending.resolve;
      pending.resolve = (r) => { clearTimeout(timer); originalResolve(r); };
      slot.worker.postMessage(pending.run, pending.transferables);
    }
    // If the pool has lost all workers (timeouts or onerror) and there's
    // still queued work, callers would otherwise wait forever. Drain.
    if (this.slots.length === 0 && this.queue.length > 0) {
      const stuck = this.queue;
      this.queue = [];
      for (const pending of stuck) {
        pending.resolve({
          kind: 'failure', taskId: pending.run.taskId,
          message: 'no workers available', failKind: 'other',
          seed: pending.run.seed,
        });
      }
    }
  }

  private removeSlot(slot: Slot, reason: string): void {
    const idx = this.slots.indexOf(slot);
    if (idx < 0) return;
    this.slots.splice(idx, 1);
    for (const ps of slot.pendingSetups.values()) {
      ps.resolve({
        kind: 'setup-ack', sessionId: ps.sessionId, ok: false,
        message: `worker removed: ${reason}`,
      });
    }
    slot.pendingSetups.clear();
    if (slot.current) {
      const cur = slot.current;
      slot.current = null;
      slot.busy = false;
      cur.resolve({
        kind: 'failure', taskId: cur.run.taskId,
        message: `worker removed: ${reason}`,
        failKind: 'other', seed: cur.run.seed,
      });
    }
    try { slot.worker.terminate(); } catch { /* ignore */ }
    this.pump();
  }

  // Send `setup` to every slot; resolve when all acks come back ok, reject
  // on the first ok=false. Setup is sent without transferables — every slot
  // receives a structured-clone copy. Transferring to one slot would detach
  // the underlying buffer before sibling slots could clone it; copying is
  // simpler and still avoids the per-seed re-encoding cost the split was
  // designed to fix.
  //
  // Each ack is bounded by `setupTimeoutMs`; on timeout the slot is
  // terminated and removed from the pool, and setupAll fails.
  async setupAll(setup: FitSessionSetup, _transferables: Transferable[]): Promise<void> {
    if (this.disposed) throw new Error('pool disposed');
    this.ensureWorkers();
    const acks: Promise<SetupAck>[] = [];
    for (const slot of [...this.slots]) {
      acks.push(new Promise<SetupAck>((resolve) => {
        const timer = setTimeout(() => {
          if (!slot.pendingSetups.has(setup.sessionId)) return;
          slot.pendingSetups.delete(setup.sessionId);
          this.removeSlot(slot, 'setup timed out');
          resolve({
            kind: 'setup-ack', sessionId: setup.sessionId, ok: false,
            message: `setup timed out after ${this.setupTimeoutMs}ms`,
          });
        }, this.setupTimeoutMs);
        slot.pendingSetups.set(setup.sessionId, {
          sessionId: setup.sessionId,
          resolve: (ack) => { clearTimeout(timer); resolve(ack); },
        });
        slot.worker.postMessage(setup);
      }));
    }
    const results = await Promise.all(acks);
    const failed = results.find((a) => !a.ok);
    if (failed && !failed.ok)
      throw new Error(`worker setup failed: ${failed.message}`);
  }

  // Fire-and-forget: tell each slot that has the session loaded to drop it.
  // Workers process drops in receipt order, so any earlier run-seed for the
  // session has already been replied to before drop is observed.
  dropSession(sessionId: SessionId): void {
    if (this.disposed) return;
    const drop: DropSession = {kind: 'drop-session', sessionId};
    for (const slot of this.slots) {
      if (!slot.sessions.has(sessionId)) continue;
      slot.sessions.delete(sessionId);
      try { slot.worker.postMessage(drop); } catch { /* ignore */ }
    }
  }

  dispatchRun(args: {
    run: Omit<RunSeed, 'taskId'>;
    transferables: Transferable[];
  }): Promise<RunReply> {
    if (this.disposed) return Promise.reject(new Error('pool disposed'));
    if (!this.slots.some((s) => s.sessions.has(args.run.sessionId))) {
      return Promise.reject(new Error(
        `no slot is primed for session ${args.run.sessionId}; call setupAll first`));
    }
    const taskId = this.nextTaskId++;
    const run: RunSeed = {...args.run, taskId} as RunSeed;
    return new Promise<RunReply>((resolve) => {
      this.queue.push({run, resolve, transferables: args.transferables});
      this.pump();
    });
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    for (const slot of this.slots) {
      // Drain in-flight setup acks so setupAll's Promise.all can settle.
      for (const ps of slot.pendingSetups.values()) {
        ps.resolve({
          kind: 'setup-ack', sessionId: ps.sessionId, ok: false,
          message: 'pool disposed before setup completed',
        });
      }
      slot.pendingSetups.clear();
      // Drain a run already dispatched but not yet replied.
      if (slot.current) {
        slot.current.resolve({
          kind: 'failure', taskId: slot.current.run.taskId,
          message: 'pool disposed mid-run', failKind: 'other',
          seed: slot.current.run.seed,
        });
        slot.current = null;
      }
      try { slot.worker.terminate(); } catch { /* ignore */ }
    }
    this.slots = [];
    // Drain queued runs that were never dispatched.
    for (const pending of this.queue) {
      pending.resolve({
        kind: 'failure',
        taskId: pending.run.taskId,
        message: 'pool disposed before task ran',
        failKind: 'other',
        seed: pending.run.seed,
      });
    }
    this.queue = [];
  }
}

export function defaultPoolSize(): number {
  const hc = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 2;
  return Math.max(2, hc - 1);
}

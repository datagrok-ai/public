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

// Slot run-state as a discriminated union. The previous encoding was a
// (busy: boolean, current: PendingRun | null, manualTimer) triple with the
// invariant `busy === (current !== null)` maintained by hand at each
// transition site. This makes that invariant structural, and binds the
// run-timeout `setTimeout` handle to the same state object, so exiting
// `running` requires producing the timer to clear — orphan timers become
// unrepresentable.
type SlotRunState =
  | {phase: 'idle'}
  | {phase: 'running'; run: PendingRun; runTimer: ReturnType<typeof setTimeout>};

type Slot = {
  worker: Worker;
  runState: SlotRunState;
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

interface ActiveSetup {
  setup: FitSessionSetup;
  // Per-session re-prime budget: incremented every time a replacement slot
  // is primed for this session in `removeSlot`. Capped so a script body
  // that fails every worker can't put the pool into a respawn-and-retry
  // loop. Once exhausted, slot replacements still spawn but skip
  // re-priming for this session — the in-flight fit runs at degraded
  // parallelism for the rest of its life.
  replacementAttempts: number;
}

export class WorkerPool {
  private static readonly MAX_REPLACEMENT_REPRIMES = 1;

  private slots: Slot[] = [];
  private queue: PendingRun[] = [];
  // Setups for sessions that have completed `setupAll` and not yet been
  // dropped. Used to re-prime replacement slots so an active fit can
  // keep using full pool parallelism after a worker dies.
  private activeSetups: Map<SessionId, ActiveSetup> = new Map();
  private disposed = false;
  private nextTaskId = 1;
  private readonly setupTimeoutMs: number;
  private readonly runTimeoutMs: number;

  constructor(private readonly size: number, opts: WorkerPoolOptions = {}) {
    this.setupTimeoutMs = opts.setupTimeoutMs ?? 20_000;
    this.runTimeoutMs = opts.runTimeoutMs ?? 60_000;
  }

  // Lazy: workers spin up on the first setupAll so a pool can be allocated
  // without paying the worker spin-up cost until a fit actually runs.
  private ensureWorkers(): void {
    if (this.slots.length > 0 || this.disposed) return;
    for (let i = 0; i < this.size; ++i)
      this.slots.push(this.spawnSlot());
  }

  private spawnSlot(): Slot {
    const w = new Worker(new URL('./fitting.worker.ts', import.meta.url));
    const slot: Slot = {
      worker: w,
      runState: {phase: 'idle'},
      sessions: new Set(),
      pendingSetups: new Map(),
    };
    w.onmessage = (ev: MessageEvent<WorkerInbound>) => this.onMessage(slot, ev.data);
    w.onerror = (ev: ErrorEvent) => this.onError(slot, ev);
    return slot;
  }

  private onMessage(slot: Slot, msg: WorkerInbound): void {
    if (msg.kind === 'setup-ack') {
      const pending = slot.pendingSetups.get(msg.sessionId);
      if (!pending) return;
      slot.pendingSetups.delete(msg.sessionId);
      if (msg.ok) {
        slot.sessions.add(msg.sessionId);
        // Pump in case runs were queued while this slot was being primed
        // (e.g., a replacement re-prime mid-fit) — without this, queued
        // runs would only retry on the next reply, which may never come
        // if the new slot is the only one able to serve the session.
        this.pump();
      }
      pending.resolve(msg);
      return;
    }
    // success / failure for a run-seed
    if (slot.runState.phase !== 'running') {
      // Stray reply (worker double-posted, or a reply arrived after the
      // timeout already fired and removed the slot). Drop it.
      return;
    }
    const {run: pending, runTimer} = slot.runState;
    clearTimeout(runTimer);
    slot.runState = {phase: 'idle'};
    pending.resolve(msg);
    this.pump();
  }

  private onError(slot: Slot, ev: ErrorEvent): void {
    if (slot.runState.phase === 'running') {
      const {run: pending, runTimer} = slot.runState;
      clearTimeout(runTimer);
      slot.runState = {phase: 'idle'};
      pending.resolve({
        kind: 'failure',
        taskId: pending.run.taskId,
        seedIndex: pending.run.seedIndex,
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
    // Place each pending run on an idle slot that has its session primed.
    // Unplaceable items stay queued; pump runs again on every setup-ack ok
    // and every run reply, so they get retried as soon as a slot is
    // primed-and-idle for that session.
    let i = 0;
    while (i < this.queue.length) {
      const pending = this.queue[i];
      const slot = this.slots.find(
        (s) => s.runState.phase === 'idle' && s.sessions.has(pending.run.sessionId));
      if (!slot) { ++i; continue; }
      this.queue.splice(i, 1);
      // Per-run timeout: terminates the slot if the worker wedges (e.g.
      // infinite loop in user JS body) so the dispatchRun promise resolves
      // instead of hanging forever. Bound the timer to this run via
      // slot.runState — if the slot transitions out of running for any
      // other reason (reply, error, removeSlot, dispose), those paths
      // clearTimeout it; the callback's identity check is a defensive
      // guard against a fired-then-cancelled race.
      const runTimer = setTimeout(() => {
        if (slot.runState.phase !== 'running' || slot.runState.run !== pending) return;
        slot.runState = {phase: 'idle'};
        this.removeSlot(slot, 'run timed out');
        pending.resolve({
          kind: 'failure', taskId: pending.run.taskId,
          seedIndex: pending.run.seedIndex,
          message: `run timed out after ${this.runTimeoutMs}ms`,
          failKind: 'other', seed: pending.run.seed,
        });
      }, this.runTimeoutMs);
      slot.runState = {phase: 'running', run: pending, runTimer};
      slot.worker.postMessage(pending.run, pending.transferables);
    }
    // Per-session stall drain: if a session has lost every primed slot
    // AND its replacement-reprime budget is exhausted AND no setup-ack
    // is in flight for it on any slot, no setup-ack-ok or run-reply will
    // ever arrive to retrigger pump for it. Without this branch, queued
    // runs for the session would sit forever — the executor's
    // `Promise.all(runSlot)` would hang. The hasInFlightSetup guard is
    // load-bearing: removeSlot increments replacementAttempts BEFORE
    // primeSlot's ack arrives, so the predicate would otherwise fire
    // during a recoverable transient state.
    if (this.queue.length > 0) {
      for (const [sessionId, entry] of this.activeSetups) {
        if (entry.replacementAttempts >= WorkerPool.MAX_REPLACEMENT_REPRIMES &&
            this.primedSlotCount(sessionId) === 0 &&
            !this.hasInFlightSetup(sessionId)) {
          this.drainSessionQueue(sessionId, 'session lost all primed slots');
        }
      }
    }
    // If the pool has lost all workers (timeouts or onerror) and there's
    // still queued work, callers would otherwise wait forever. Drain.
    if (this.slots.length === 0 && this.queue.length > 0) {
      const stuck = this.queue;
      this.queue = [];
      for (const pending of stuck) {
        pending.resolve({
          kind: 'failure', taskId: pending.run.taskId,
          seedIndex: pending.run.seedIndex,
          message: 'no workers available', failKind: 'other',
          seed: pending.run.seed,
        });
      }
    }
  }

  // Number of live slots whose `slot.sessions` Set contains the given
  // sessionId. Derived state — slot count is bounded (≤ pool size,
  // typically 2–8) and Set.has is O(1), so the linear scan is cheap.
  private primedSlotCount(sessionId: SessionId): number {
    let n = 0;
    for (const slot of this.slots) if (slot.sessions.has(sessionId)) ++n;
    return n;
  }

  // True iff any slot has a pending setup-ack for this session. Used by
  // the stall-drain predicate to avoid firing during a transient state
  // where a replacement re-prime is in flight.
  private hasInFlightSetup(sessionId: SessionId): boolean {
    for (const slot of this.slots) {
      if (slot.pendingSetups.has(sessionId)) return true;
    }
    return false;
  }

  // Resolve every queued run for the given session as a failure with the
  // supplied reason. Used by the per-session stall drain in `pump` and
  // available for future use cases (e.g. session abort).
  private drainSessionQueue(sessionId: SessionId, reason: string): void {
    const remaining: PendingRun[] = [];
    for (const pending of this.queue) {
      if (pending.run.sessionId !== sessionId) {
        remaining.push(pending);
        continue;
      }
      pending.resolve({
        kind: 'failure',
        taskId: pending.run.taskId,
        seedIndex: pending.run.seedIndex,
        message: reason,
        failKind: 'other',
        seed: pending.run.seed,
      });
    }
    this.queue = remaining;
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
    if (slot.runState.phase === 'running') {
      const {run: cur, runTimer} = slot.runState;
      clearTimeout(runTimer);
      slot.runState = {phase: 'idle'};
      cur.resolve({
        kind: 'failure', taskId: cur.run.taskId,
        seedIndex: cur.run.seedIndex,
        message: `worker removed: ${reason}`,
        failKind: 'other', seed: cur.run.seed,
      });
    }
    try { slot.worker.terminate(); } catch { /* ignore */ }
    if (!this.disposed) {
      // Spawn a fresh replacement so a long-lived pool doesn't drift toward
      // zero workers under repeated errors/timeouts.
      const fresh = this.spawnSlot();
      this.slots.push(fresh);
      // Re-prime the replacement for every active session so the in-flight
      // fit can keep using full pool parallelism. Capped per session at
      // MAX_REPLACEMENT_REPRIMES to break a pathological respawn loop on
      // a script body that breaks every worker — once the budget is
      // exhausted, replacements still spawn but skip re-priming for that
      // session and the fit runs at degraded parallelism.
      for (const entry of this.activeSetups.values()) {
        if (entry.replacementAttempts >= WorkerPool.MAX_REPLACEMENT_REPRIMES) continue;
        ++entry.replacementAttempts;
        // Fire-and-forget — primeSlot's own timeout path handles failure
        // (calls removeSlot, which respects the cap).
        void this.primeSlot(fresh, entry.setup);
      }
    }
    this.pump();
  }

  // Post a setup to one slot and resolve when its ack arrives or the
  // setup timeout fires. Shared by `setupAll` (initial priming of every
  // slot) and `removeSlot` (re-priming a freshly spawned replacement).
  // On timeout we remove the slot, which itself may spawn another
  // replacement — bounded by `MAX_REPLACEMENT_REPRIMES` per session.
  private primeSlot(slot: Slot, setup: FitSessionSetup): Promise<SetupAck> {
    return new Promise<SetupAck>((resolve) => {
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
    });
  }

  // Send `setup` to every slot; resolve when all acks come back ok, reject
  // on the first ok=false. Setup is sent without transferables — every slot
  // gets a structured-clone copy. Transferring to one slot would detach the
  // underlying buffer before sibling slots could clone it.
  //
  // Each ack is bounded by `setupTimeoutMs`; on timeout the slot is
  // terminated and removed from the pool, and setupAll fails.
  async setupAll(setup: FitSessionSetup, _transferables: Transferable[]): Promise<void> {
    if (this.disposed) throw new Error('pool disposed');
    this.ensureWorkers();
    const results = await Promise.all(
      [...this.slots].map((slot) => this.primeSlot(slot, setup)));
    const failed = results.find((a) => !a.ok);
    if (failed && !failed.ok)
      throw new Error(`worker setup failed: ${failed.message}`);
    // Remember the setup so a slot replacement (after onerror / timeout)
    // can re-prime a fresh worker without involving the executor layer.
    this.activeSetups.set(setup.sessionId, {setup, replacementAttempts: 0});
  }

  // Fire-and-forget: tell each slot that has the session loaded to drop it.
  // Workers process drops in receipt order, so any earlier run-seed for the
  // session has already been replied to before drop is observed.
  dropSession(sessionId: SessionId): void {
    if (this.disposed) return;
    this.activeSetups.delete(sessionId);
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
    // Reject only if the session was never set up. If it WAS set up but
    // no slot currently has it primed (e.g. a replacement re-prime is in
    // flight), let the run sit in the queue — `pump` is called on every
    // setup-ack ok and on every run reply, so the run will dispatch
    // when a primed slot is available again.
    if (!this.activeSetups.has(args.run.sessionId)) {
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
    this.activeSetups.clear();
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
      if (slot.runState.phase === 'running') {
        const {run: cur, runTimer} = slot.runState;
        clearTimeout(runTimer);
        cur.resolve({
          kind: 'failure', taskId: cur.run.taskId,
          seedIndex: cur.run.seedIndex,
          message: 'pool disposed mid-run', failKind: 'other',
          seed: cur.run.seed,
        });
        slot.runState = {phase: 'idle'};
      }
      try { slot.worker.terminate(); } catch { /* ignore */ }
    }
    this.slots = [];
    // Drain queued runs that were never dispatched.
    for (const pending of this.queue) {
      pending.resolve({
        kind: 'failure',
        taskId: pending.run.taskId,
        seedIndex: pending.run.seedIndex,
        message: 'pool disposed before task ran',
        failKind: 'other',
        seed: pending.run.seed,
      });
    }
    this.queue = [];
  }

  // Test-only: read live slot count without exposing the array.
  _slotsForTest(): number {
    return this.slots.length;
  }

  // Test-only: how many slots currently have the session in their primed
  // set. Used by white-box stall tests to observe progress through the
  // remove → reprime cycle without reaching into the internals via casts.
  _primedSlotCountForTest(sessionId: SessionId): number {
    return this.primedSlotCount(sessionId);
  }
}

export function defaultPoolSize(): number {
  const hc = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 2;
  return Math.max(2, hc - 1);
}

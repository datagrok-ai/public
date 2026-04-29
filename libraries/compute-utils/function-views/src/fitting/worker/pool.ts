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
import {Query} from './query';
import {makeSetupFailure} from './failures';

export type RunReply = WorkerSuccess | WorkerFailure;
export {Query} from './query';

interface PendingSetup {
  sessionId: SessionId;
  resolve: (ack: SetupAck) => void;
}

// Slot run-state as a discriminated union. The running variant holds a
// reference to the in-flight Query, which owns its own run-timeout timer;
// transitioning out of `running` happens via takeRunningQuery() and lets
// the caller settle the query (which clears the timer idempotently).
type SlotRunState =
  | {phase: 'idle'}
  | {phase: 'running'; query: Query};

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
  /** Per-session cap on replacement re-primes. Each setupAll'd session
   *  starts at 0 and increments every time `removeSlot` re-primes a
   *  fresh replacement slot for it. Once the cap is hit, further slot
   *  replacements still spawn but skip re-priming for this session,
   *  trading parallelism for safety against a script body that breaks
   *  every worker. Default 3 — covers a few transient deaths (browser
   *  process eviction, OOM hiccup, GC pause crossing setupTimeoutMs)
   *  without enabling a runaway respawn loop, since each failed
   *  reprime takes ≤ setupTimeoutMs to surface. */
  maxReplacementReprimes?: number;
}

interface ActiveSetup {
  setup: FitSessionSetup;
  // Per-session re-prime budget: incremented every time a replacement slot
  // is primed for this session in `removeSlot`. Capped by
  // maxReplacementReprimes so a script body that fails every worker
  // can't put the pool into a respawn-and-retry loop.
  replacementAttempts: number;
}

export class WorkerPool {
  private slots: Slot[] = [];
  private queue: Query[] = [];
  // Setups for sessions that have completed `setupAll` and not yet been
  // dropped. Used to re-prime replacement slots so an active fit can
  // keep using full pool parallelism after a worker dies.
  private activeSetups: Map<SessionId, ActiveSetup> = new Map();
  private disposed = false;
  private nextTaskId = 1;
  private readonly setupTimeoutMs: number;
  private readonly runTimeoutMs: number;
  private readonly maxReplacementReprimes: number;

  constructor(private readonly size: number, opts: WorkerPoolOptions = {}) {
    this.setupTimeoutMs = opts.setupTimeoutMs ?? 20_000;
    this.runTimeoutMs = opts.runTimeoutMs ?? 60_000;
    this.maxReplacementReprimes = opts.maxReplacementReprimes ?? 3;
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
    const {query} = slot.runState;
    slot.runState = {phase: 'idle'};
    query.settle(msg);
    this.pump();
  }

  private onError(slot: Slot, ev: ErrorEvent): void {
    const message = `worker error: ${ev.message ?? 'unknown'}`;
    if (slot.runState.phase === 'running') {
      const {query} = slot.runState;
      slot.runState = {phase: 'idle'};
      query.fail(message);
    }
    // Reject any in-flight setup acks for this slot too.
    for (const ps of slot.pendingSetups.values())
      ps.resolve(makeSetupFailure(ps.sessionId, message));
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
      const query = this.queue[i];
      const slot = this.slots.find(
        (s) => s.runState.phase === 'idle' && s.sessions.has(query.sessionId));
      if (!slot) { ++i; continue; }
      this.queue.splice(i, 1);
      // Per-run timeout: terminates the slot if the worker wedges so the
      // caller's promise resolves instead of hanging. Query.settle is
      // idempotent, so a reply arriving after the timeout fires is
      // harmless (and vice versa) — no defensive identity check needed.
      const runTimer = setTimeout(() => {
        if (query.phase !== 'running') return;
        if (slot.runState.phase === 'running' && slot.runState.query === query)
          slot.runState = {phase: 'idle'};
        this.removeSlot(slot, 'run timed out');
        query.fail(`run timed out after ${this.runTimeoutMs}ms`);
      }, this.runTimeoutMs);
      query.startRunning(runTimer);
      slot.runState = {phase: 'running', query};
      slot.worker.postMessage(query.run, query.transferables);
    }
    // Drain sessions whose queues can't make progress. pump only retries
    // on setup-ack-ok and run-reply, neither of which can arrive when no
    // slot is primed and the reprime budget is exhausted. hasInFlightSetup
    // is load-bearing — replacementAttempts increments synchronously
    // before primeSlot's ack returns.
    if (this.queue.length > 0) {
      for (const [sessionId, entry] of this.activeSetups) {
        if (entry.replacementAttempts >= this.maxReplacementReprimes &&
            this.primedSlotCount(sessionId) === 0 &&
            !this.hasInFlightSetup(sessionId)) {
          this.drainQueue('session lost all primed slots',
            (q) => q.sessionId === sessionId);
        }
      }
    }
    // If the pool has lost all workers (timeouts or onerror) and there's
    // still queued work, callers would otherwise wait forever. Drain.
    if (this.slots.length === 0 && this.queue.length > 0)
      this.drainQueue('no workers available', () => true);
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

  // Resolve every queued run that matches `match` as a failure with the
  // supplied reason. Used by pump's per-session stall drain and the
  // total-pool drain, and by dispose to clear the queue at teardown.
  private drainQueue(reason: string, match: (q: Query) => boolean): void {
    const remaining: Query[] = [];
    for (const query of this.queue) {
      if (!match(query)) { remaining.push(query); continue; }
      query.fail(reason);
    }
    this.queue = remaining;
  }

  private removeSlot(slot: Slot, reason: string): void {
    const idx = this.slots.indexOf(slot);
    if (idx < 0) return;
    this.slots.splice(idx, 1);
    const message = `worker removed: ${reason}`;
    for (const ps of slot.pendingSetups.values())
      ps.resolve(makeSetupFailure(ps.sessionId, message));
    slot.pendingSetups.clear();
    if (slot.runState.phase === 'running') {
      const {query} = slot.runState;
      slot.runState = {phase: 'idle'};
      query.fail(message);
    }
    try { slot.worker.terminate(); } catch { /* ignore */ }
    if (!this.disposed) {
      // Spawn a fresh replacement so a long-lived pool doesn't drift toward
      // zero workers under repeated errors/timeouts.
      const fresh = this.spawnSlot();
      this.slots.push(fresh);
      // Re-prime the replacement for every active session so the in-flight
      // fit can keep using full pool parallelism. Capped per session by
      // maxReplacementReprimes — once exhausted, replacements still spawn
      // but skip re-priming for that session and the fit runs at degraded
      // parallelism (and pump's stall-drain branch fails any queued runs
      // for it to keep the executor from hanging).
      for (const entry of this.activeSetups.values()) {
        if (entry.replacementAttempts >= this.maxReplacementReprimes) continue;
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
        resolve(makeSetupFailure(setup.sessionId,
          `setup timed out after ${this.setupTimeoutMs}ms`));
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
      this.queue.push(new Query(run, args.transferables, resolve));
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
        ps.resolve(makeSetupFailure(ps.sessionId,
          'pool disposed before setup completed'));
      }
      slot.pendingSetups.clear();
      // Drain a run already dispatched but not yet replied.
      if (slot.runState.phase === 'running') {
        const {query} = slot.runState;
        slot.runState = {phase: 'idle'};
        query.fail('pool disposed mid-run');
      }
      try { slot.worker.terminate(); } catch { /* ignore */ }
    }
    this.slots = [];
    // Drain queued runs that were never dispatched.
    this.drainQueue('pool disposed before task ran', () => true);
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

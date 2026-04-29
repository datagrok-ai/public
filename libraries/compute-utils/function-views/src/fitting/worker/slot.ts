// One worker slot in the pool.
//
// A Slot owns a Worker plus the state needed to dispatch one run at a time
// and to track multiple concurrent setup-acks (one per session being
// primed on this slot). The state machine is:
//
//   idle ─ claim(query) ─→ running(query) ─ takeRunningQuery() ─→ idle
//
// Setup tracking is independent of run state — a slot can be running a
// query for session A while a setup-ack for session B is in flight on the
// same slot.
//
// The pool drives transitions via method calls; reaching into slot fields
// is an anti-pattern (and unrepresentable for test code, since `runState`
// is private). Test fixtures that need to inject in-flight state use the
// `_setRunningForTest` test-only method.

import type {FitSessionSetup, WorkerInbound, SetupAck, SessionId} from './wire-types';
import {Query} from './query';
import {makeSetupFailure} from './failures';

interface PendingSetup {
  sessionId: SessionId;
  resolve: (ack: SetupAck) => void;
  setupTimer: ReturnType<typeof setTimeout>;
}

type RunState = {phase: 'idle'} | {phase: 'running'; query: Query};

export class Slot {
  private runState: RunState = {phase: 'idle'};
  /** Sessions that have been successfully primed on this slot's worker. */
  readonly sessions = new Set<SessionId>();
  /** Per-session setup-acks currently awaiting reply. */
  readonly pendingSetups = new Map<SessionId, PendingSetup>();
  private readonly worker: Worker;

  constructor(
    onMessage: (msg: WorkerInbound) => void,
    onError: (ev: ErrorEvent) => void,
  ) {
    // The literal `new Worker(new URL('./fitting.worker.ts', import.meta.url))`
    // pattern is what webpack's worker bundler recognizes — passing the URL
    // in as a parameter breaks static analysis and the worker never lands in
    // its own bundle chunk. So Slot owns the entry path directly.
    this.worker = new Worker(new URL('./fitting.worker.ts', import.meta.url));
    this.worker.onmessage = (ev: MessageEvent<WorkerInbound>) => onMessage(ev.data);
    this.worker.onerror = onError;
  }

  // ---- queries ----

  isIdle(): boolean { return this.runState.phase === 'idle'; }
  hasSession(id: SessionId): boolean { return this.sessions.has(id); }
  hasPendingSetup(id: SessionId): boolean { return this.pendingSetups.has(id); }
  currentQuery(): Query | null {
    return this.runState.phase === 'running' ? this.runState.query : null;
  }

  // ---- run dispatch ----

  /**
   * Atomic transition idle → running. The caller hands in the query and a
   * run-timeout callback; we arm the timer (handing it to the query so it
   * survives idempotent settlement), transition state, and post the message
   * to the worker.
   */
  claim(query: Query, runTimeoutMs: number, onTimeout: () => void): void {
    const runTimer = setTimeout(onTimeout, runTimeoutMs);
    query.startRunning(runTimer);
    this.runState = {phase: 'running', query};
    this.worker.postMessage(query.run, query.transferables);
  }

  /**
   * Take the running query for resolution, transitioning back to idle.
   * Returns null if the slot is already idle (stray reply or already-handled
   * timeout / removeSlot path); caller treats null as "drop on the floor".
   */
  takeRunningQuery(): Query | null {
    if (this.runState.phase !== 'running') return null;
    const {query} = this.runState;
    this.runState = {phase: 'idle'};
    return query;
  }

  // ---- setup ----

  /**
   * Post a setup; resolve when the ack arrives or when the timeout fires.
   * On timeout, the slot is no longer usable for this session — the pool's
   * onTimeout callback typically removes the slot.
   */
  prime(setup: FitSessionSetup, setupTimeoutMs: number, onTimeout: () => void): Promise<SetupAck> {
    return new Promise<SetupAck>((resolve) => {
      const setupTimer = setTimeout(() => {
        if (!this.pendingSetups.has(setup.sessionId)) return;
        this.pendingSetups.delete(setup.sessionId);
        onTimeout();
        resolve(makeSetupFailure(setup.sessionId,
          `setup timed out after ${setupTimeoutMs}ms`));
      }, setupTimeoutMs);
      this.pendingSetups.set(setup.sessionId, {
        sessionId: setup.sessionId,
        resolve: (ack) => { clearTimeout(setupTimer); resolve(ack); },
        setupTimer,
      });
      this.worker.postMessage(setup);
    });
  }

  /**
   * Pop the pending setup for this session (if any), clearing its timer.
   * Returns the entry so the caller can resolve it with the actual ack.
   * Returns null for stray acks.
   */
  consumeSetupAck(sessionId: SessionId): PendingSetup | null {
    const ps = this.pendingSetups.get(sessionId);
    if (!ps) return null;
    this.pendingSetups.delete(sessionId);
    clearTimeout(ps.setupTimer);
    return ps;
  }

  /** Record that the worker has acked a successful setup for this session. */
  recordPrimed(sessionId: SessionId): void {
    this.sessions.add(sessionId);
  }

  /**
   * Tell the worker to drop per-session state for `sessionId` and remove
   * the session from the primed set. No-op if the session wasn't primed.
   */
  dropSession(sessionId: SessionId): void {
    if (!this.sessions.has(sessionId)) return;
    this.sessions.delete(sessionId);
    try {
      this.worker.postMessage({kind: 'drop-session', sessionId});
    } catch { /* ignore */ }
  }

  // ---- teardown ----

  /**
   * Fail every in-flight setup-ack on this slot with the given message.
   * Used by onError (worker crashed) and removeSlot (slot being torn down).
   */
  failPendingSetups(message: string): void {
    for (const ps of this.pendingSetups.values()) {
      clearTimeout(ps.setupTimer);
      ps.resolve(makeSetupFailure(ps.sessionId, message));
    }
    this.pendingSetups.clear();
  }

  /** Terminate the underlying Worker. Safe to call multiple times. */
  terminate(): void {
    try { this.worker.terminate(); } catch { /* ignore */ }
  }

  // ---- test-only ----

  /**
   * Inject a query into the running slot for white-box tests. The query
   * must already be primed (e.g. via `query.startRunning(noopTimer)` or by
   * being constructed in the test). Calling this on a non-idle slot
   * silently overwrites the prior state — tests are responsible for slot
   * sequencing.
   */
  _setRunningForTest(query: Query): void {
    this.runState = {phase: 'running', query};
  }
}

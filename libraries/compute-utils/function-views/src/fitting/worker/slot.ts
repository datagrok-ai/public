// One worker slot: owns a Worker, a single run-state, and a per-session
// pending-setup map. All transitions are pool-driven.

import type {FitSessionSetup, RunDispatch, WorkerInbound, SetupAck, SessionId} from './wire-types';
import {RunJob} from './run-job';
import type {RunReply} from './pool';
import {makeSetupFailure, makeSetupTimeoutFailure} from './failures';

interface PendingSetup {
  sessionId: SessionId;
  resolve: (ack: SetupAck) => void;
  setupTimer: ReturnType<typeof setTimeout>;
}

type RunState = {phase: 'idle'} | {phase: 'running'; job: RunJob};

export class Slot {
  private runState: RunState = {phase: 'idle'};
  readonly sessions = new Set<SessionId>();
  readonly pendingSetups = new Map<SessionId, PendingSetup>();
  private readonly worker: Worker;

  constructor(
    onMessage: (msg: WorkerInbound) => void,
    onError: (ev: ErrorEvent) => void,
  ) {
    // Webpack's worker bundler only recognizes this exact `new Worker(new URL(...))`
    // literal; parameterizing the URL breaks static analysis.
    this.worker = new Worker(new URL('./fitting.worker.ts', import.meta.url));
    this.worker.onmessage = (ev: MessageEvent<WorkerInbound>) => onMessage(ev.data);
    this.worker.onerror = onError;
  }

  isIdle(): boolean { return this.runState.phase === 'idle'; }
  hasSession(id: SessionId): boolean { return this.sessions.has(id); }
  hasPendingSetup(id: SessionId): boolean { return this.pendingSetups.has(id); }

  claim(job: RunJob, runTimeoutMs: number, onTimeout: () => void): void {
    const runTimer = setTimeout(onTimeout, runTimeoutMs);
    job.startRunning(runTimer);
    this.runState = {phase: 'running', job};
    const run: RunDispatch = {kind: 'run-dispatch', ...job.spec};
    this.worker.postMessage(run, job.transferables);
  }

  // Settle the in-flight job with this reply, transition to idle.
  // Returns false (no-op) on stray reply (slot already idle).
  deliverReply(reply: RunReply): boolean {
    const job = this.takeRunning();
    if (!job) return false;
    job.settle(reply);
    return true;
  }

  // Fail the in-flight job (if any) with this message, transition to idle.
  failRunning(message: string): void {
    this.takeRunning()?.fail(message);
  }

  private takeRunning(): RunJob | null {
    if (this.runState.phase !== 'running') return null;
    const {job} = this.runState;
    this.runState = {phase: 'idle'};
    return job;
  }

  prime(setup: FitSessionSetup, setupTimeoutMs: number): Promise<SetupAck> {
    return new Promise<SetupAck>((resolve) => {
      const setupTimer = setTimeout(() => {
        if (!this.pendingSetups.has(setup.sessionId)) return;
        this.pendingSetups.delete(setup.sessionId);
        resolve(makeSetupTimeoutFailure(setup.sessionId, setupTimeoutMs));
      }, setupTimeoutMs);
      this.pendingSetups.set(setup.sessionId, {
        sessionId: setup.sessionId,
        resolve: (ack) => { clearTimeout(setupTimer); resolve(ack); },
        setupTimer,
      });
      this.worker.postMessage(setup);
    });
  }

  // Returns null for stray acks (already timed out / never registered).
  completeSetupAck(sessionId: SessionId, ok: boolean): PendingSetup | null {
    const ps = this.pendingSetups.get(sessionId);
    if (!ps) return null;
    this.pendingSetups.delete(sessionId);
    clearTimeout(ps.setupTimer);
    if (ok) this.sessions.add(sessionId);
    return ps;
  }

  dropSession(sessionId: SessionId): void {
    if (!this.sessions.has(sessionId)) return;
    this.sessions.delete(sessionId);
    try {
      this.worker.postMessage({kind: 'drop-session', sessionId});
    } catch { /* ignore */ }
  }

  failPendingSetups(message: string): void {
    for (const ps of this.pendingSetups.values()) {
      clearTimeout(ps.setupTimer);
      ps.resolve(makeSetupFailure(ps.sessionId, message));
    }
    this.pendingSetups.clear();
  }

  terminate(): void {
    try { this.worker.terminate(); } catch { /* ignore */ }
  }

  // ---- test-only ----

  _setRunningForTest(job: RunJob): void {
    this.runState = {phase: 'running', job};
  }
}

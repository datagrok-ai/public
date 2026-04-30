// One worker slot. Owns a Worker, a single run-state (idle | running),
// and a per-session pending-setup map (multiple sessions can be primed
// concurrently). All transitions go through pool callbacks; tests inject
// state via the `_*ForTest` hooks below.

import type {FitSessionSetup, RunDispatch, WorkerInbound, SetupAck, SessionId} from './wire-types';
import {RunJob} from './run-job';
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
    // The literal `new Worker(new URL('./fitting.worker.ts', import.meta.url))`
    // pattern is what webpack's worker bundler recognizes; passing the URL
    // through a parameter would break static analysis and fail to emit the
    // worker chunk.
    this.worker = new Worker(new URL('./fitting.worker.ts', import.meta.url));
    this.worker.onmessage = (ev: MessageEvent<WorkerInbound>) => onMessage(ev.data);
    this.worker.onerror = onError;
  }

  isIdle(): boolean { return this.runState.phase === 'idle'; }
  hasSession(id: SessionId): boolean { return this.sessions.has(id); }
  hasPendingSetup(id: SessionId): boolean { return this.pendingSetups.has(id); }
  currentJob(): RunJob | null {
    return this.runState.phase === 'running' ? this.runState.job : null;
  }

  claim(job: RunJob, runTimeoutMs: number, onTimeout: () => void): void {
    const runTimer = setTimeout(onTimeout, runTimeoutMs);
    job.startRunning(runTimer);
    this.runState = {phase: 'running', job};
    const run: RunDispatch = {kind: 'run-dispatch', ...job.spec};
    this.worker.postMessage(run, job.transferables);
  }

  takeRunningJob(): RunJob | null {
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

  // Pop the pending setup for `sessionId`, clear its timer, and (on ok)
  // mark the session primed. Returns null for stray acks.
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

  /** Inject a running job; tests are responsible for slot sequencing. */
  _setRunningForTest(job: RunJob): void {
    this.runState = {phase: 'running', job};
  }
}

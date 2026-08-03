// Worker pool with two-phase dispatch: setupAll → dispatchRun → dropSession.
// Workers spin lazily on first setupAll and stay alive across fits so the
// in-worker compile cache is reused.

import type {
  FitSessionSetup,
  JobSpec,
  SetupAck,
  WorkerInbound,
  WorkerSuccess,
  WorkerFailure,
  SessionId,
} from './wire-types';
import {RunJob} from './run-job';
import {Slot} from './slot';

export type RunReply = WorkerSuccess | WorkerFailure;
export {RunJob} from './run-job';

/** Bounds on every promise from the pool, so a wedged worker can't hang
 *  a fit indefinitely. */
export interface WorkerPoolOptions {
  /** Reject setupAll's per-slot ack after this duration. Default: 20000. */
  setupTimeoutMs?: number;
  /** Reject one Nelder-Mead run after this duration. Default: 60000. */
  runTimeoutMs?: number;
  /** Per-session cap on replacement re-primes. Guards against a script body
   *  that breaks every worker by skipping reprime once exhausted. Default: 3. */
  maxReplacementReprimes?: number;
}

interface ActiveSetup {
  setup: FitSessionSetup;
  replacementAttempts: number;
}

class JobQueue {
  private items: RunJob[] = [];

  enqueue(job: RunJob): void { this.items.push(job); }
  isEmpty(): boolean { return this.items.length === 0; }

  // Walk once. If `claim` returns true, remove that item; else leave it.
  tryDispatch(claim: (job: RunJob) => boolean): void {
    let i = 0;
    while (i < this.items.length) {
      if (claim(this.items[i])) this.items.splice(i, 1);
      else ++i;
    }
  }

  // Remove and fail every queued item matching `match`.
  drain(reason: string, match: (job: RunJob) => boolean): void {
    const kept: RunJob[] = [];
    for (const job of this.items) {
      if (match(job)) job.fail(reason);
      else kept.push(job);
    }
    this.items = kept;
  }
}

export class WorkerPool {
  private slots: Slot[] = [];
  private queue = new JobQueue();
  private activeSetups: Map<SessionId, ActiveSetup> = new Map();
  private disposed = false;
  private readonly setupTimeoutMs: number;
  private readonly runTimeoutMs: number;
  private readonly maxReplacementReprimes: number;

  constructor(private readonly size: number, opts: WorkerPoolOptions = {}) {
    this.setupTimeoutMs = opts.setupTimeoutMs ?? 20_000;
    this.runTimeoutMs = opts.runTimeoutMs ?? 60_000;
    this.maxReplacementReprimes = opts.maxReplacementReprimes ?? 3;
  }

  private ensureWorkers(): void {
    if (this.slots.length > 0 || this.disposed) return;
    for (let i = 0; i < this.size; ++i)
      this.slots.push(this.spawnSlot());
  }

  private spawnSlot(): Slot {
    const slot: Slot = new Slot(
      (msg) => this.onMessage(slot, msg),
      (ev) => this.onError(slot, ev),
    );
    return slot;
  }

  private onMessage(slot: Slot, msg: WorkerInbound): void {
    if (msg.kind === 'setup-ack') {
      const ps = slot.completeSetupAck(msg.sessionId, msg.ok);
      if (!ps) return;
      // Pump on every ok ack — without this, runs queued while a replacement
      // was being primed would only retry on the next reply, which may
      // never come if this slot is the only one able to serve the session.
      if (msg.ok) this.pump();
      ps.resolve(msg);
      return;
    }
    if (!slot.deliverReply(msg)) return;   // stray reply (post-timeout race)
    this.pump();
  }

  private onError(slot: Slot, ev: ErrorEvent): void {
    const reason = ev.message ?? 'unknown';
    console.warn(`worker error: ${reason}`);
    // A crash is only script-attributable when the worker dies mid-run (e.g.
    // OOM); setup-phase script errors come back as ok:false acks, not onerror.
    const culprit = slot.runningSessionId() ?? undefined;
    this.removeSlot(slot, `worker errored: ${reason}`, culprit);
  }

  private pump(): void {
    if (this.disposed) return;
    this.queue.tryDispatch((job) => {
      const slot = this.slots.find((s) => s.isIdle() && s.hasSession(job.spec.sessionId));
      if (!slot) return false;
      slot.claim(job, this.runTimeoutMs, () => {
        if (job.phase !== 'running') return;
        slot.failRunning(`run timed out after ${this.runTimeoutMs}ms`);
        this.removeSlot(slot, 'run timed out', job.spec.sessionId);
      });
      return true;
    });
    // Stall drain: with no primed slot and the reprime budget exhausted,
    // pump won't fire again. `hasInFlightSetup` guards against draining
    // mid-reprime (replacementAttempts increments before the ack lands).
    if (!this.queue.isEmpty()) {
      for (const [sessionId, entry] of this.activeSetups) {
        if (entry.replacementAttempts >= this.maxReplacementReprimes &&
            this.primedSlotCount(sessionId) === 0 &&
            !this.hasInFlightSetup(sessionId)) {
          this.queue.drain('session lost all primed slots',
            (q) => q.spec.sessionId === sessionId);
        }
      }
    }
    if (this.slots.length === 0 && !this.queue.isEmpty())
      this.queue.drain('no workers available', () => true);
  }

  private primedSlotCount(sessionId: SessionId): number {
    let n = 0;
    for (const slot of this.slots) if (slot.hasSession(sessionId)) ++n;
    return n;
  }

  private hasInFlightSetup(sessionId: SessionId): boolean {
    for (const slot of this.slots) {
      if (slot.hasPendingSetup(sessionId)) return true;
    }
    return false;
  }

  // `culprit` is the session whose run/setup caused this removal (undefined for
  // an unattributed infra crash). Only the culprit's reprime budget is charged,
  // so a misbehaving session can't drain an innocent concurrent session's budget.
  private removeSlot(slot: Slot, reason: string, culprit?: SessionId): void {
    const idx = this.slots.indexOf(slot);
    if (idx < 0) return;
    this.slots.splice(idx, 1);
    const message = `worker removed: ${reason}`;
    slot.failPendingSetups(message);
    slot.failRunning(message);
    slot.terminate();
    if (!this.disposed) {
      const fresh = this.spawnSlot();
      this.slots.push(fresh);
      // Re-prime every active session that hasn't been given up on. Only the
      // culprit is charged + cap-gated; that cap stops the timeout-recurse loop
      // for a script that breaks its own workers.
      for (const [sid, entry] of this.activeSetups) {
        if (entry.replacementAttempts >= this.maxReplacementReprimes) continue;
        if (sid === culprit) ++entry.replacementAttempts;
        void this.primeSlot(fresh, entry.setup);
      }
    }
    this.pump();
  }

  // Prime a slot and remove it on setup-timeout. Caller decides what to do
  // with the returned ack (await + throw, or fire-and-forget).
  private async primeSlot(slot: Slot, setup: FitSessionSetup): Promise<SetupAck> {
    const ack = await slot.prime(setup, this.setupTimeoutMs);
    if (ack.ok === false && ack.timedOut) this.removeSlot(slot, 'setup timed out', setup.sessionId);
    return ack;
  }

  // No transferables: each slot needs its own structured-clone copy.
  // Transferring would detach the buffer before siblings could clone it.
  async setupAll(setup: FitSessionSetup): Promise<void> {
    if (this.disposed) throw new Error('pool disposed');
    this.ensureWorkers();
    const acks = await Promise.all([...this.slots].map((s) => this.primeSlot(s, setup)));
    for (const ack of acks) {
      if (ack.ok === false) throw new Error(`worker setup failed: ${ack.message}`);
    }
    this.activeSetups.set(setup.sessionId, {setup, replacementAttempts: 0});
  }

  // Fire-and-forget; workers process drops in receipt order, so earlier
  // dispatches for the session are replied to before drop is observed.
  dropSession(sessionId: SessionId): void {
    if (this.disposed) return;
    this.activeSetups.delete(sessionId);
    for (const slot of this.slots) slot.dropSession(sessionId);
  }

  dispatchRun(spec: JobSpec): Promise<RunReply> {
    if (this.disposed) return Promise.reject(new Error('pool disposed'));
    // If the session was set up but no slot has it primed (mid-reprime), the
    // run sits in the queue — pump retries on every ack-ok / run-reply.
    if (!this.activeSetups.has(spec.sessionId)) {
      return Promise.reject(new Error(
        `no slot is primed for session ${spec.sessionId}; call setupAll first`));
    }
    return new Promise<RunReply>((resolve) => {
      this.queue.enqueue(new RunJob(spec, [spec.seed.buffer], resolve));
      this.pump();
    });
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.activeSetups.clear();
    for (const slot of this.slots) {
      slot.failPendingSetups('pool disposed before setup completed');
      slot.failRunning('pool disposed mid-run');
      slot.terminate();
    }
    this.slots = [];
    this.queue.drain('pool disposed before task ran', () => true);
  }

  // ---- test-only: live slot count, white-box hooks for replacement / drain tests ----

  _slotsForTest(): number { return this.slots.length; }
  _ensureWorkersForTest(): void { this.ensureWorkers(); }
  _slotsArrayForTest(): readonly Slot[] { return this.slots; }
  _primedSlotCountForTest(sessionId: SessionId): number { return this.primedSlotCount(sessionId); }
  _removeSlotForTest(slot: Slot, reason: string, culprit?: SessionId): void {
    this.removeSlot(slot, reason, culprit);
  }
  _replacementAttemptsForTest(sessionId: SessionId): number | undefined {
    return this.activeSetups.get(sessionId)?.replacementAttempts;
  }
}

export function defaultPoolSize(): number {
  const hc = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 2;
  return Math.max(2, hc - 1);
}

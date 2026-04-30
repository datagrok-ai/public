// Worker pool with two-phase dispatch. `setupAll` ships an immutable bundle
// to every slot once per fit; `dispatchRun` ships a seed and routes to a
// slot whose primed-session set contains it; `dropSession` releases the
// session at fit teardown. Workers are spun lazily on first setupAll and
// stay alive across fits so the in-worker compile cache is reused.
//
// State is split across Slot (one Worker + run-state + setup tracking),
// RunJob (one queued/in-flight dispatch + its result promise), and
// WorkerPool (slot allocator + session registry + pump).

import type {
  FitSessionSetup,
  JobSpec,
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

export class WorkerPool {
  private slots: Slot[] = [];
  private queue: RunJob[] = [];
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
    // Run reply. takeRunningJob returns null for a stray reply (worker
    // double-posted, or a reply arriving after the timeout removed the slot).
    const job = slot.takeRunningJob();
    if (!job) return;
    job.settle(msg);
    this.pump();
  }

  private onError(slot: Slot, ev: ErrorEvent): void {
    const reason = ev.message ?? 'unknown';
    console.warn(`worker error: ${reason}`);
    this.removeSlot(slot, `worker errored: ${reason}`);
  }

  private pump(): void {
    if (this.disposed) return;
    let i = 0;
    while (i < this.queue.length) {
      const job = this.queue[i];
      const slot = this.slots.find((s) => s.isIdle() && s.hasSession(job.spec.sessionId));
      if (!slot) { ++i; continue; }
      this.queue.splice(i, 1);
      // RunJob.settle is idempotent, so a reply landing after the timeout
      // (or vice versa) is harmless.
      slot.claim(job, this.runTimeoutMs, () => {
        if (job.phase !== 'running') return;
        if (slot.currentJob() === job) slot.takeRunningJob();
        this.removeSlot(slot, 'run timed out');
        job.fail(`run timed out after ${this.runTimeoutMs}ms`);
      });
    }
    // Stall drain: pump only retries on setup-ack-ok and run-reply, neither
    // of which arrives once a session has lost all primed slots and the
    // reprime budget is exhausted. `hasInFlightSetup` is load-bearing —
    // replacementAttempts increments synchronously before primeSlot acks.
    if (this.queue.length > 0) {
      for (const [sessionId, entry] of this.activeSetups) {
        if (entry.replacementAttempts >= this.maxReplacementReprimes &&
            this.primedSlotCount(sessionId) === 0 &&
            !this.hasInFlightSetup(sessionId)) {
          this.drainQueue('session lost all primed slots',
            (q) => q.spec.sessionId === sessionId);
        }
      }
    }
    if (this.slots.length === 0 && this.queue.length > 0)
      this.drainQueue('no workers available', () => true);
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

  private drainQueue(reason: string, match: (q: RunJob) => boolean): void {
    const remaining: RunJob[] = [];
    for (const job of this.queue) {
      if (!match(job)) { remaining.push(job); continue; }
      job.fail(reason);
    }
    this.queue = remaining;
  }

  private removeSlot(slot: Slot, reason: string): void {
    const idx = this.slots.indexOf(slot);
    if (idx < 0) return;
    this.slots.splice(idx, 1);
    const message = `worker removed: ${reason}`;
    slot.failPendingSetups(message);
    const inflight = slot.takeRunningJob();
    if (inflight) inflight.fail(message);
    slot.terminate();
    if (!this.disposed) {
      const fresh = this.spawnSlot();
      this.slots.push(fresh);
      // Re-prime the replacement for every active session, capped per session
      // by maxReplacementReprimes. If the reprime itself times out, the .then
      // recurses through removeSlot — the cap is what stops the loop.
      for (const entry of this.activeSetups.values()) {
        if (entry.replacementAttempts >= this.maxReplacementReprimes) continue;
        ++entry.replacementAttempts;
        void fresh.prime(entry.setup, this.setupTimeoutMs).then((ack) => {
          if (!ack.ok && ack.timedOut) this.removeSlot(fresh, 'setup timed out');
        });
      }
    }
    this.pump();
  }

  // Send `setup` to every slot and resolve when all acks come back ok; throw
  // on the first failure. Setup is sent without transferables: each slot
  // gets its own structured-clone copy because transferring would detach
  // the buffer before sibling slots could clone it.
  async setupAll(setup: FitSessionSetup): Promise<void> {
    if (this.disposed) throw new Error('pool disposed');
    this.ensureWorkers();
    const targets = [...this.slots];
    const acks = await Promise.all(
      targets.map((slot) => slot.prime(setup, this.setupTimeoutMs)));
    for (let i = 0; i < acks.length; ++i) {
      const ack = acks[i];
      if (!ack.ok) {
        if (ack.timedOut) this.removeSlot(targets[i], 'setup timed out');
        throw new Error(`worker setup failed: ${ack.message}`);
      }
    }
    this.activeSetups.set(setup.sessionId, {setup, replacementAttempts: 0});
  }

  // Fire-and-forget. Workers process drops in receipt order, so any earlier
  // run-seed for the session is replied to before the drop is observed.
  dropSession(sessionId: SessionId): void {
    if (this.disposed) return;
    this.activeSetups.delete(sessionId);
    for (const slot of this.slots) slot.dropSession(sessionId);
  }

  dispatchRun(spec: JobSpec): Promise<RunReply> {
    if (this.disposed) return Promise.reject(new Error('pool disposed'));
    // Reject only if the session was never set up. If it was set up but no
    // slot currently has it primed (e.g. a reprime is in flight), the run
    // sits in the queue — pump fires on every ack-ok / run-reply.
    if (!this.activeSetups.has(spec.sessionId)) {
      return Promise.reject(new Error(
        `no slot is primed for session ${spec.sessionId}; call setupAll first`));
    }
    return new Promise<RunReply>((resolve) => {
      this.queue.push(new RunJob(spec, [spec.seed.buffer], resolve));
      this.pump();
    });
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.activeSetups.clear();
    for (const slot of this.slots) {
      slot.failPendingSetups('pool disposed before setup completed');
      const inflight = slot.takeRunningJob();
      if (inflight) inflight.fail('pool disposed mid-run');
      slot.terminate();
    }
    this.slots = [];
    this.drainQueue('pool disposed before task ran', () => true);
  }

  // ---- test-only: live slot count, white-box hooks for replacement / drain tests ----

  _slotsForTest(): number { return this.slots.length; }
  _ensureWorkersForTest(): void { this.ensureWorkers(); }
  _slotsArrayForTest(): readonly Slot[] { return this.slots; }
  _primedSlotCountForTest(sessionId: SessionId): number { return this.primedSlotCount(sessionId); }
  _removeSlotForTest(slot: Slot, reason: string): void { this.removeSlot(slot, reason); }
  _replacementAttemptsForTest(sessionId: SessionId): number | undefined {
    return this.activeSetups.get(sessionId)?.replacementAttempts;
  }
}

export function defaultPoolSize(): number {
  const hc = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 2;
  return Math.max(2, hc - 1);
}

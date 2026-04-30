// One queued or in-flight worker dispatch.
//
// `Query` owns the resolve callback that wakes the dispatchRun caller, plus
// the per-run timeout handle once the query is dispatched. Settling is
// idempotent so a stray reply arriving after the timeout fired (or a
// removeSlot fired after the reply landed) doesn't double-resolve.

import type {RunSeed} from './wire-types';
import {makeRunFailure} from './failures';
import type {RunReply} from './pool';

type QueryState =
  | {phase: 'queued'}
  | {phase: 'running'; runTimer: ReturnType<typeof setTimeout>}
  | {phase: 'settled'};

export class Query {
  private state: QueryState = {phase: 'queued'};

  constructor(
    readonly run: RunSeed,
    readonly transferables: Transferable[],
    private readonly resolveCb: (r: RunReply) => void,
  ) {}

  get sessionId(): number { return this.run.sessionId; }
  get taskId(): number { return this.run.taskId; }
  get phase(): QueryState['phase'] { return this.state.phase; }

  /** idle → running. Caller hands in the timer it just armed. */
  startRunning(runTimer: ReturnType<typeof setTimeout>): void {
    this.state = {phase: 'running', runTimer};
  }

  /**
   * Resolve the caller's promise and clear the run timer if armed. Idempotent —
   * second-and-later calls are no-ops, so race-prone callers (timeout vs reply,
   * removeSlot vs reply) don't need a defensive identity check.
   */
  settle(reply: RunReply): void {
    if (this.state.phase === 'settled') return;
    if (this.state.phase === 'running') clearTimeout(this.state.runTimer);
    this.state = {phase: 'settled'};
    this.resolveCb(reply);
  }

  /** Convenience for the failure path — settle with a generic-other failure. */
  fail(message: string): void {
    this.settle(makeRunFailure(this.run, message));
  }
}

// One queued or in-flight worker dispatch. `settle` is idempotent — late
// timeouts and late replies race, both can hit safely.

import type {JobSpec} from './wire-types';
import {makeRunFailure} from './failures';
import type {RunReply} from './pool';

type JobState =
  | {phase: 'queued'}
  | {phase: 'running'; runTimer: ReturnType<typeof setTimeout>}
  | {phase: 'settled'};

export class RunJob {
  private state: JobState = {phase: 'queued'};

  constructor(
    readonly spec: JobSpec,
    readonly transferables: Transferable[],
    private readonly resolveCb: (r: RunReply) => void,
  ) {}

  get phase(): JobState['phase'] { return this.state.phase; }

  startRunning(runTimer: ReturnType<typeof setTimeout>): void {
    this.state = {phase: 'running', runTimer};
  }

  settle(reply: RunReply): void {
    if (this.state.phase === 'settled') return;
    if (this.state.phase === 'running') clearTimeout(this.state.runTimer);
    this.state = {phase: 'settled'};
    this.resolveCb(reply);
  }

  fail(message: string): void {
    this.settle(makeRunFailure(this.spec.seedIndex, this.spec.seed, message));
  }
}

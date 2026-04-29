// Failure-shape builders.
//
// Every site in the worker pool that fails a queued run or a pending
// setup-ack constructs the same wire-types literal. Centralizing the
// builders here keeps the contract in one place and makes "forgot a
// field" a type error rather than a silent runtime bug.

import type {RunSeed, SetupAck, WorkerFailure, SessionId} from './wire-types';

export function makeRunFailure(run: RunSeed, message: string): WorkerFailure {
  return {
    kind: 'failure',
    taskId: run.taskId,
    seedIndex: run.seedIndex,
    message,
    failKind: 'other',
    seed: run.seed,
  };
}

export function makeSetupFailure(sessionId: SessionId, message: string): SetupAck {
  return {kind: 'setup-ack', sessionId, ok: false, message};
}

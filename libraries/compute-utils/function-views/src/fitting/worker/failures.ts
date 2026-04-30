// Failure-shape builders. Centralized so a missing field is a type error.

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

export function makeSetupTimeoutFailure(sessionId: SessionId, ms: number): SetupAck {
  return {kind: 'setup-ack', sessionId, ok: false, message: `setup timed out after ${ms}ms`, timedOut: true};
}

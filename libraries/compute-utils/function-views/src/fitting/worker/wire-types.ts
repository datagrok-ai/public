// Wire types for main ↔ worker messages.
//
// Two-message protocol: a fit's immutable bundle (script body, fixed inputs,
// targets, NM settings) ships once per worker as `FitSessionSetup`; per-seed
// dispatch sends only `{sessionId, seed}` via `RunSeed`. Worker memory is
// freed at fit teardown via `DropSession`.
//
// Pure type declarations + the LOSS string-enum re-export. Zero runtime
// imports beyond `../constants`. Kept separate from serialize.ts so the
// worker entry can import these types without pulling serialize.ts's
// DG-coupled dependencies into the worker bundle.

import type {LOSS} from '../constants';
import type {ValueBoundsData} from '../optimizer-misc';

export type SessionId = number;

export type SerializedScalarTarget = {
  kind: 'scalar';
  propName: string;
  type: 'int' | 'bigint' | 'float';
  target: number;
};

export type SerializedDataFrameTarget = {
  kind: 'dataframe';
  propName: string;
  argName: string;
  arrowIPC: Uint8Array;
  funcColNames: string[];
};

export type SerializedOutputTarget = SerializedScalarTarget | SerializedDataFrameTarget;

// Outbound — sent once per worker per fit.
export type FitSessionSetup = {
  kind: 'setup-fit';
  sessionId: SessionId;
  fnSource: string;
  paramList: string[];
  outputParamNames: string[];
  lossType: LOSS;
  fixedInputs: Record<string, number | string | boolean | null>;
  fixedDataFrames: Record<string, Uint8Array>;
  variedInputNames: string[];
  bounds: Record<string, ValueBoundsData>;
  outputTargets: SerializedOutputTarget[];
  nmSettings: [string, number][];
  threshold?: number;
};

// Outbound — sent once per seed.
export type RunSeed = {
  kind: 'run-seed';
  taskId: number;
  sessionId: SessionId;
  seed: Float64Array;
};

// Outbound — sent once per fit at teardown. Lets the worker free its
// per-session state. Fire-and-forget, no ack expected. Harmless on
// ephemeral pools (the worker is about to be terminated anyway); required
// on long-lived pools so memory doesn't grow per fit.
export type DropSession = {
  kind: 'drop-session';
  sessionId: SessionId;
};

export type WorkerOutbound = FitSessionSetup | RunSeed | DropSession;

// Inbound — setup acknowledgement. `ok: false` carries a compile or
// Arrow-decode error.
export type SetupAck =
  | { kind: 'setup-ack'; sessionId: SessionId; ok: true }
  | { kind: 'setup-ack'; sessionId: SessionId; ok: false; message: string };

// Inbound — run results.
export type WorkerSuccess = {
  kind: 'success';
  taskId: number;
  point: Float64Array;
  cost: number;
  iterCosts: number[];
  iterCount: number;
};

export type WorkerFailure = {
  kind: 'failure';
  taskId: number;
  message: string;
  failKind: 'inconsistent' | 'other';
  seed: Float64Array;
};

export type WorkerInbound = SetupAck | WorkerSuccess | WorkerFailure;

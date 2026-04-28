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

// Sparse type tag for non-scalar fixedInputs that need worker-side
// reification. Plain scalars (number, string, boolean, null) omit
// themselves from the map entirely — only Dayjs and Date get tagged so
// the worker knows to wrap the ISO string back into the original object.
export type FixedInputKind = 'dayjs' | 'date';

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
  fixedInputTypes: Record<string, FixedInputKind>;
  fixedDataFrames: Record<string, Uint8Array>;
  variedInputNames: string[];
  bounds: Record<string, ValueBoundsData>;
  outputTargets: SerializedOutputTarget[];
  nmSettings: [string, number][];
  threshold?: number;
};

// Outbound — sent once per seed. `seedIndex` is the position in the per-fit
// `params` array; the worker echoes it on every reply so the executor can
// place results in deterministic seed-index order regardless of completion
// order across workers.
export type RunSeed = {
  kind: 'run-seed';
  taskId: number;
  sessionId: SessionId;
  seedIndex: number;
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

// Inbound — run results. `seedIndex` echoes the dispatch's `RunSeed.seedIndex`
// so the executor can record replies by index rather than completion order.
export type WorkerSuccess = {
  kind: 'success';
  taskId: number;
  seedIndex: number;
  point: Float64Array;
  cost: number;
  iterCosts: number[];
  iterCount: number;
};

export type WorkerFailure = {
  kind: 'failure';
  taskId: number;
  seedIndex: number;
  message: string;
  failKind: 'inconsistent' | 'other';
  seed: Float64Array;
};

export type WorkerInbound = SetupAck | WorkerSuccess | WorkerFailure;

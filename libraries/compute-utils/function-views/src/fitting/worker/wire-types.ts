// Wire types for main ↔ worker messages. Pure type declarations — kept
// separate from serialize.ts so the worker entry can import these without
// pulling serialize.ts's DG-coupled deps into the worker bundle.

import type {LOSS} from '../constants';
import type {ValueBoundsData} from '../optimizer-misc';

export type SessionId = number;

// Sparse tag for non-scalar fixedInputs needing worker-side reification.
// Plain scalars are absent from the map; only Dayjs/Date carry a tag.
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

// Outbound — sent once per seed. The worker echoes `seedIndex` on every
// reply so the executor can place results in deterministic order.
export type RunSeed = {
  kind: 'run-seed';
  taskId: number;
  sessionId: SessionId;
  seedIndex: number;
  seed: Float64Array;
};

// Outbound — fire-and-forget at teardown. Required on long-lived pools so
// per-session memory doesn't grow per fit.
export type DropSession = {
  kind: 'drop-session';
  sessionId: SessionId;
};

export type WorkerOutbound = FitSessionSetup | RunSeed | DropSession;

// Inbound — setup acknowledgement. `ok: false` carries a compile or
// Arrow-decode error. `timedOut` is set only by the local pool path (never
// posted by the worker) so the pool can distinguish ack-timeouts from
// worker-reported failures.
export type SetupAck =
  | { kind: 'setup-ack'; sessionId: SessionId; ok: true }
  | { kind: 'setup-ack'; sessionId: SessionId; ok: false; message: string; timedOut?: true };

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

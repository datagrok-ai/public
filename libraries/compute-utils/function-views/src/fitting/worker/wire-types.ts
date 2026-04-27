// Wire types for main ↔ worker messages.
//
// Pure type declarations + the LOSS string-enum re-export. Zero runtime
// imports beyond `../constants` (which itself has none). Kept in a separate
// file from serialize.ts so the worker entry can import these types without
// pulling serialize.ts's DG-coupled dependencies into the worker bundle.

import type {LOSS} from '../constants';
import type {ValueBoundsData} from '../optimizer-misc';

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

export type ObjectiveTask = {
  kind: 'objective';
  taskId: number;
  source: string;
  seed: Float64Array;
  nmSettings: [string, number][];
  threshold?: number;
};

export type FitTask = {
  kind: 'fit';
  taskId: number;
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
  seed: Float64Array;
};

export type WorkerTask = ObjectiveTask | FitTask;

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

export type WorkerReply = WorkerSuccess | WorkerFailure;

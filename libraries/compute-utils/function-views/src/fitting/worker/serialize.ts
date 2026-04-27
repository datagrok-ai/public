// Wire format and serializer for main → worker dispatch.
//
// Two task kinds cross the boundary:
//
//   - 'objective' — generic Float64Array → number objective. Used by tests
//     that fit pure-math objectives without any DG.Func / FuncCall plumbing.
//   - 'fit'       — full fitting payload. Mirrors the closure built by
//     cost-functions.ts:makeConstFunction, but pieces it apart so the
//     worker can rebuild the cost function locally from the source body and
//     deserialized DataFrames.
//
// DataFrames cross as Arrow IPC bytes via @datagrok-libraries/arrow toFeather;
// the worker decodes via arrow-to-lite.
//
// Day.js dates → ISO string. Reconstructed in the worker if needed.

import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {toFeather} from '@datagrok-libraries/arrow';
import {LOSS} from '../constants';
import type {OutputTargetItem, ValueBoundsData} from '../optimizer-misc';
import type {
  FitTask, ObjectiveTask, SerializedDataFrameTarget, SerializedOutputTarget,
  SerializedScalarTarget, WorkerFailure, WorkerReply, WorkerSuccess, WorkerTask,
} from './wire-types';

export type {
  FitTask, ObjectiveTask, SerializedDataFrameTarget, SerializedOutputTarget,
  SerializedScalarTarget, WorkerFailure, WorkerReply, WorkerSuccess, WorkerTask,
};

function serializeFixedInputs(
  fixedInputs: Record<string, any>,
): {scalars: Record<string, any>; dataFrames: Record<string, Uint8Array>; transferables: ArrayBuffer[]} {
  const scalars: Record<string, any> = {};
  const dataFrames: Record<string, Uint8Array> = {};
  const transferables: ArrayBuffer[] = [];
  for (const [name, value] of Object.entries(fixedInputs)) {
    if (value instanceof DG.DataFrame) {
      const bytes = toFeather(value);
      if (!bytes) throw new Error(`failed to serialize fixed dataframe input "${name}"`);
      dataFrames[name] = bytes;
      transferables.push(bytes.buffer as ArrayBuffer);
      continue;
    }
    if (dayjs.isDayjs(value)) {
      scalars[name] = value.toISOString();
      continue;
    }
    if (value instanceof Date) {
      scalars[name] = value.toISOString();
      continue;
    }
    scalars[name] = value;
  }
  return {scalars, dataFrames, transferables};
}

function serializeOutputTargets(
  targets: OutputTargetItem[],
): {targets: SerializedOutputTarget[]; transferables: ArrayBuffer[]} {
  const out: SerializedOutputTarget[] = [];
  const transferables: ArrayBuffer[] = [];
  for (const t of targets) {
    if (t.type === DG.TYPE.DATA_FRAME) {
      const bytes = toFeather(t.target as DG.DataFrame);
      if (!bytes) throw new Error(`failed to serialize target dataframe for "${t.propName}"`);
      transferables.push(bytes.buffer as ArrayBuffer);
      out.push({
        kind: 'dataframe',
        propName: t.propName,
        argName: t.argName,
        arrowIPC: bytes,
        funcColNames: t.cols.map((c) => c.name),
      });
    } else {
      let scalarType: 'int' | 'bigint' | 'float' = 'float';
      if (t.type === DG.TYPE.INT) scalarType = 'int';
      else if (t.type === DG.TYPE.BIG_INT) scalarType = 'bigint';
      out.push({kind: 'scalar', propName: t.propName, type: scalarType, target: t.target as number});
    }
  }
  return {targets: out, transferables};
}

export function buildObjectiveTask(args: {
  taskId: number;
  source: string;
  seed: Float64Array;
  nmSettings: Map<string, number>;
  threshold?: number;
}): {task: ObjectiveTask; transferables: Transferable[]} {
  const task: ObjectiveTask = {
    kind: 'objective',
    taskId: args.taskId,
    source: args.source,
    seed: args.seed,
    nmSettings: Array.from(args.nmSettings.entries()),
    threshold: args.threshold,
  };
  return {task, transferables: [args.seed.buffer]};
}

export function buildFitTask(args: {
  taskId: number;
  fnSource: string;
  paramList: string[];
  outputParamNames: string[];
  lossType: LOSS;
  fixedInputs: Record<string, any>;
  variedInputNames: string[];
  bounds: Record<string, ValueBoundsData>;
  outputTargets: OutputTargetItem[];
  nmSettings: Map<string, number>;
  threshold?: number;
  seed: Float64Array;
}): {task: FitTask; transferables: Transferable[]} {
  const fixed = serializeFixedInputs(args.fixedInputs);
  const targets = serializeOutputTargets(args.outputTargets);
  const task: FitTask = {
    kind: 'fit',
    taskId: args.taskId,
    fnSource: args.fnSource,
    paramList: args.paramList,
    outputParamNames: args.outputParamNames,
    lossType: args.lossType,
    fixedInputs: fixed.scalars,
    fixedDataFrames: fixed.dataFrames,
    variedInputNames: args.variedInputNames,
    bounds: args.bounds,
    outputTargets: targets.targets,
    nmSettings: Array.from(args.nmSettings.entries()),
    threshold: args.threshold,
    seed: args.seed,
  };
  return {
    task,
    transferables: [args.seed.buffer, ...fixed.transferables, ...targets.transferables],
  };
}

// Cap on accepted body size — guards against pasted blobs and runaway parse
// time. Plan §1.3 calls for 64 KB.
export const MAX_FN_SOURCE_BYTES = 64 * 1024;

export function checkSourceSize(source: string): void {
  if (source.length > MAX_FN_SOURCE_BYTES)
    throw new Error(`function source exceeds ${MAX_FN_SOURCE_BYTES} bytes (got ${source.length})`);
}

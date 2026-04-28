// Serializer for main → worker dispatch.
//
// Two messages cross the boundary:
//
//   - 'setup-fit' — full immutable bundle (script body, fixed inputs and
//     fixed DataFrames, output targets, bounds, NM settings, threshold).
//     Sent once per worker per fit.
//   - 'run-seed' — seed + sessionId pointer. Sent once per seed.
//
// DataFrames cross as Arrow IPC bytes via @datagrok-libraries/arrow toFeather;
// the worker decodes via arrow-to-lite once when it receives the setup.
//
// Day.js dates → ISO string. Reconstructed in the worker if needed.

import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {toFeather} from '@datagrok-libraries/arrow';
import {LOSS} from '../constants';
import type {OutputTargetItem, ValueBoundsData} from '../optimizer-misc';
import type {
  FitSessionSetup, RunSeed, SerializedDataFrameTarget, SerializedOutputTarget,
  SerializedScalarTarget, SessionId, SetupAck, WorkerFailure, WorkerInbound, WorkerOutbound,
  WorkerSuccess,
} from './wire-types';

export type {
  FitSessionSetup, RunSeed, SerializedDataFrameTarget, SerializedOutputTarget,
  SerializedScalarTarget, SessionId, SetupAck, WorkerFailure, WorkerInbound, WorkerOutbound,
  WorkerSuccess,
};

function serializeFixedInputs(
  fixedInputs: Record<string, any>,
): {
  scalars: Record<string, any>;
  dataFrames: Record<string, Uint8Array>;
  types: Record<string, 'dayjs' | 'date'>;
  transferables: ArrayBuffer[];
} {
  const scalars: Record<string, any> = {};
  const dataFrames: Record<string, Uint8Array> = {};
  const types: Record<string, 'dayjs' | 'date'> = {};
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
      types[name] = 'dayjs';
      continue;
    }
    if (value instanceof Date) {
      scalars[name] = value.toISOString();
      types[name] = 'date';
      continue;
    }
    scalars[name] = value;
  }
  return {scalars, dataFrames, types, transferables};
}

// Const VALUES live in fixedInputs/fixedDataFrames; setup.bounds keeps the
// structural shape (which entries are const vs varying, formula bound
// strings) only. Replacing const value with `null` keeps `getAccData`'s
// constValues count and ordering intact while making the wire safe.
function stripBoundsConstValues(
  bounds: Record<string, ValueBoundsData>,
): Record<string, ValueBoundsData> {
  const out: Record<string, ValueBoundsData> = {};
  for (const [name, entry] of Object.entries(bounds)) {
    if (entry.type === 'const')
      out[name] = {type: 'const', value: null};
    else
      out[name] = entry;
  }
  return out;
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

export function buildSetup(args: {
  sessionId: SessionId;
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
}): {setup: FitSessionSetup; transferables: Transferable[]} {
  const fixed = serializeFixedInputs(args.fixedInputs);
  const targets = serializeOutputTargets(args.outputTargets);
  const setup: FitSessionSetup = {
    kind: 'setup-fit',
    sessionId: args.sessionId,
    fnSource: args.fnSource,
    paramList: args.paramList,
    outputParamNames: args.outputParamNames,
    lossType: args.lossType,
    fixedInputs: fixed.scalars,
    fixedInputTypes: fixed.types,
    fixedDataFrames: fixed.dataFrames,
    variedInputNames: args.variedInputNames,
    bounds: stripBoundsConstValues(args.bounds),
    outputTargets: targets.targets,
    nmSettings: Array.from(args.nmSettings.entries()),
    threshold: args.threshold,
  };
  return {
    setup,
    transferables: [...fixed.transferables, ...targets.transferables],
  };
}

export function buildRunSeed(args: {
  sessionId: SessionId;
  taskId: number;
  seedIndex: number;
  seed: Float64Array;
}): {run: RunSeed; transferables: Transferable[]} {
  const run: RunSeed = {
    kind: 'run-seed',
    taskId: args.taskId,
    sessionId: args.sessionId,
    seedIndex: args.seedIndex,
    seed: args.seed,
  };
  return {run, transferables: [args.seed.buffer]};
}

// Cap on accepted body size — guards against pasted blobs and runaway
// parse time.
export const MAX_FN_SOURCE_BYTES = 64 * 1024;

export function checkSourceSize(source: string): void {
  if (source.length > MAX_FN_SOURCE_BYTES)
    throw new Error(`function source exceeds ${MAX_FN_SOURCE_BYTES} bytes (got ${source.length})`);
}

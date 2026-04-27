// Fitting worker entry. Long-lived; receives one WorkerTask per message,
// runs Nelder-Mead on it, posts WorkerReply. Stays alive across messages so
// the compile cache (func-call-shim.ts) hits across tasks.
//
// Imports are deliberately DG-free:
//   - `optimizer-nelder-mead` exports only `optimizeNM`; its type-only import
//     of optimizer-misc is erased by tsc, so the compiled JS has no DG ref.
//   - `dg-shim`, `arrow-to-lite`, `cost-math`, `func-call-shim`, `serialize`
//     are pure TypeScript with no datagrok-api/* runtime imports.
//
// Anything that touches DG.* (fitting-utils.ts, optimizer-misc.ts runtime,
// cost-functions.ts) is intentionally NOT imported here — it would crash a
// worker the moment the bundle loaded.

import {optimizeNM} from '../optimizer-nelder-mead';
import {LOSS} from '../constants';
import {arrowIpcToLite} from './arrow-to-lite';
import {createWorkerDG} from './dg-shim';
import type {LiteDataFrame} from './types';
import {ColLike, DfLike, getErrors, InconsistentTablesError} from './cost-math';
import {compileBody, createWorkerFuncCall} from './func-call-shim';
import type {
  ObjectiveTask,
  FitTask,
  WorkerTask,
  WorkerReply,
  SerializedOutputTarget,
  SerializedDataFrameTarget,
} from './wire-types';

// Mirror of the makeBoundsChecker logic from optimizer-sampler.ts, inlined to
// avoid pulling in formulas-resolver (which compiles user formula expressions
// via runtime helpers we don't need in-worker since pure-numeric bounds cover
// the worker fitting use case). Worker bodies that contain formula bounds
// would not pass canHandle() on the main thread, so in practice only value
// bounds reach this path.
function makeWorkerBoundsChecker(
  bounds: Record<string, any>,
  variedInputNames: string[],
): (point: Float64Array) => boolean {
  const ranges: Array<{idx: number; min: number; max: number}> = [];
  variedInputNames.forEach((name, idx) => {
    const b = bounds[name];
    if (!b || b.type !== 'changing') return;
    if (b.bottom?.type !== 'value' || b.top?.type !== 'value') {
      // Formula bound — give up and accept all points (canHandle should have
      // already routed this fit to the main arm).
      ranges.push({idx, min: -Infinity, max: Infinity});
      return;
    }
    ranges.push({idx, min: b.bottom.value, max: b.top.value});
  });
  return (point: Float64Array) => {
    for (const {idx, min, max} of ranges) {
      const v = point[idx];
      if (v == null || v < min || v > max) return false;
    }
    return true;
  };
}

// Adapter that lets cost-math.getErrors read a target DG.DataFrame
// reconstituted from its serialized form. The serialized form keeps only the
// arg column plus the listed func columns; we wrap them in a structural
// DfLike.
type TargetEntry = {
  propName: string;
  argName: string;
  argCol: ColLike;
  funcCols: ColLike[];
  df: DfLike;
};

function buildTargetEntries(targets: SerializedOutputTarget[]): {
  scalars: Array<{propName: string; target: number}>;
  dataFrames: TargetEntry[];
} {
  const scalars: Array<{propName: string; target: number}> = [];
  const dataFrames: TargetEntry[] = [];
  for (const t of targets) {
    if (t.kind === 'scalar') {
      scalars.push({propName: t.propName, target: t.target});
      continue;
    }
    const dft = t as SerializedDataFrameTarget;
    const lite = arrowIpcToLite(dft.arrowIPC);
    const argCol = lite.col(dft.argName);
    if (!argCol)
      throw new InconsistentTablesError(`target dataframe is missing arg column "${dft.argName}"`);
    const funcCols: ColLike[] = [];
    for (const cn of dft.funcColNames) {
      const c = lite.col(cn);
      if (!c)
        throw new InconsistentTablesError(`target dataframe is missing func column "${cn}"`);
      funcCols.push(c as ColLike);
    }
    dataFrames.push({
      propName: dft.propName,
      argName: dft.argName,
      argCol: argCol as ColLike,
      funcCols,
      df: lite as DfLike,
    });
  }
  return {scalars, dataFrames};
}

function reifyFixedDataFrames(blobs: Record<string, Uint8Array>): Record<string, LiteDataFrame> {
  const out: Record<string, LiteDataFrame> = {};
  for (const [name, bytes] of Object.entries(blobs))
    out[name] = arrowIpcToLite(bytes);
  return out;
}

function buildCostFunc(task: FitTask): (x: Float64Array) => number | undefined {
  const targets = buildTargetEntries(task.outputTargets);
  const fixedDfs = reifyFixedDataFrames(task.fixedDataFrames);
  const fc = createWorkerFuncCall({
    source: task.fnSource,
    paramList: task.paramList,
    outputNames: task.outputParamNames,
    fixedInputs: {...task.fixedInputs, ...fixedDfs},
    variedInputNames: task.variedInputNames,
  });
  const variedNames = task.variedInputNames;
  const checker = makeWorkerBoundsChecker(task.bounds, variedNames);
  const useRmse = task.lossType === LOSS.RMSE;

  return (x: Float64Array): number | undefined => {
    if (!checker(x)) return undefined;
    const varied: Record<string, number> = {};
    for (let i = 0; i < variedNames.length; ++i)
      varied[variedNames[i]] = x[i];
    fc.call(varied);

    if (!useRmse) {
      let mad = 0;
      for (const s of targets.scalars) {
        const sim = fc.getParamValue(s.propName) as number;
        mad = Math.max(mad, Math.abs(s.target - sim));
      }
      for (const df of targets.dataFrames) {
        const sim = fc.getParamValue(df.propName) as DfLike;
        const errs = getErrors(df.argCol, df.funcCols, sim, false);
        for (let i = 0; i < errs.length; ++i)
          mad = Math.max(mad, Math.abs(errs[i]));
      }
      return mad;
    }

    let sumSq = 0;
    let count = 0;
    for (const s of targets.scalars) {
      const sim = fc.getParamValue(s.propName) as number;
      const cur = s.target;
      sumSq += ((cur - sim) / (cur !== 0 ? cur : 1)) ** 2;
      ++count;
    }
    for (const df of targets.dataFrames) {
      const sim = fc.getParamValue(df.propName) as DfLike;
      const errs = getErrors(df.argCol, df.funcCols, sim, true);
      for (let i = 0; i < errs.length; ++i) {
        sumSq += errs[i] ** 2;
        ++count;
      }
    }
    return Math.sqrt(sumSq / count);
  };
}

function buildObjective(task: ObjectiveTask): (x: Float64Array) => number | undefined {
  const fn = compileBody(`return (${task.source}).call(null, x);`, ['x'], []) as
    (dg: unknown, x: Float64Array) => any;
  const dg = createWorkerDG();
  return (x: Float64Array): number | undefined => {
    const v = fn(dg, x);
    if (v == null || Number.isNaN(v)) return undefined;
    return v as number;
  };
}

async function run(task: WorkerTask): Promise<WorkerReply> {
  const settings = new Map<string, number>(task.nmSettings);
  try {
    const cost: (x: Float64Array) => number | undefined =
      task.kind === 'objective' ? buildObjective(task) : buildCostFunc(task);
    const objective = async (x: Float64Array): Promise<number | undefined> => cost(x);
    const ext = await optimizeNM(objective, task.seed, settings, task.threshold);
    return {
      kind: 'success',
      taskId: task.taskId,
      point: ext.point,
      cost: ext.cost,
      iterCosts: ext.iterCosts,
      iterCount: ext.iterCount,
    };
  } catch (e: any) {
    const msg = e instanceof Error ? e.message : String(e);
    const failKind = (e && e.name === 'InconsistentTables') ? 'inconsistent' : 'other';
    return {
      kind: 'failure',
      taskId: task.taskId,
      message: msg,
      failKind,
      seed: task.seed,
    };
  }
}

(self as any).onmessage = async (event: MessageEvent<WorkerTask>) => {
  const reply = await run(event.data);
  // Float64Arrays inside reply structured-clone fine; success path can also
  // transfer point.buffer to avoid the extra copy.
  const transferables: Transferable[] = [];
  if (reply.kind === 'success') transferables.push(reply.point.buffer);
  else transferables.push(reply.seed.buffer);
  (self as any).postMessage(reply, transferables);
};

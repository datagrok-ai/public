// Fitting worker entry. Long-lived; receives `FitSessionSetup` once per fit
// to prime per-session state, then `RunSeed` per seed; `DropSession`
// releases per-session state at fit teardown.
//
// Imports are deliberately DG-free:
//   - `optimizer-nelder-mead` exports only `optimizeNM`; its type-only import
//     of optimizer-misc is erased by tsc, so the compiled JS has no DG ref.
//   - `bounds-checker` is the shared (with the main-arm sampler) bounds
//     checker; it imports only `formulas-resolver` (also DG-free) and uses
//     `import type` for the bound shapes.
//   - `dg-shim`, `arrow-to-lite`, `cost-math`, `func-call-shim` are pure
//     TypeScript with no datagrok-api/* runtime imports.
//
// Anything that touches DG.* (fitting-utils.ts, optimizer-misc.ts runtime,
// cost-functions.ts) is intentionally NOT imported here — it would crash a
// worker the moment the bundle loaded.

import {optimizeNM} from '../optimizer-nelder-mead';
import {LOSS} from '../constants';
import {makeBoundsChecker} from '../bounds-checker';
import {arrowIpcToLite} from './arrow-to-lite';
import type {LiteDataFrame} from './types';
import {ColLike, DfLike, getErrors, InconsistentTablesError} from './cost-math';
import {createWorkerFuncCall, WorkerFuncCall} from './func-call-shim';
import type {
  FitSessionSetup,
  RunSeed,
  DropSession,
  WorkerOutbound,
  SetupAck,
  WorkerSuccess,
  WorkerFailure,
  SerializedOutputTarget,
  SerializedDataFrameTarget,
  SessionId,
} from './wire-types';

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

// Per-fit state, primed once on `setup-fit` and reused across all
// `run-seed` calls for that session.
type Session = {
  costFunc: (x: Float64Array) => number | undefined;
  nmSettings: Map<string, number>;
  threshold?: number;
  // Hold a strong reference to the funcCall so V8 keeps the compiled
  // body and decoded LiteDataFrames alive across runs.
  fc: WorkerFuncCall;
};

const sessions: Map<SessionId, Session> = new Map();

function buildCostFunc(setup: FitSessionSetup): {
  cost: (x: Float64Array) => number | undefined;
  fc: WorkerFuncCall;
} {
  const targets = buildTargetEntries(setup.outputTargets);
  const fixedDfs = reifyFixedDataFrames(setup.fixedDataFrames);
  const fc = createWorkerFuncCall({
    source: setup.fnSource,
    paramList: setup.paramList,
    outputNames: setup.outputParamNames,
    fixedInputs: {...setup.fixedInputs, ...fixedDfs},
    variedInputNames: setup.variedInputNames,
  });
  const variedNames = setup.variedInputNames;
  const checker = makeBoundsChecker(setup.bounds, variedNames);
  const useRmse = setup.lossType === LOSS.RMSE;

  const cost = (x: Float64Array): number | undefined => {
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
  return {cost, fc};
}

function handleSetup(setup: FitSessionSetup): SetupAck {
  try {
    const {cost, fc} = buildCostFunc(setup);
    sessions.set(setup.sessionId, {
      costFunc: cost,
      nmSettings: new Map<string, number>(setup.nmSettings),
      threshold: setup.threshold,
      fc,
    });
    return {kind: 'setup-ack', sessionId: setup.sessionId, ok: true};
  } catch (e: any) {
    const msg = e instanceof Error ? e.message : String(e);
    return {kind: 'setup-ack', sessionId: setup.sessionId, ok: false, message: msg};
  }
}

async function handleRun(run: RunSeed): Promise<WorkerSuccess | WorkerFailure> {
  const session = sessions.get(run.sessionId);
  if (!session) {
    return {
      kind: 'failure',
      taskId: run.taskId,
      message: `unknown session ${run.sessionId}`,
      failKind: 'other',
      seed: run.seed,
    };
  }
  try {
    const objective = async (x: Float64Array): Promise<number | undefined> => session.costFunc(x);
    const ext = await optimizeNM(objective, run.seed, session.nmSettings, session.threshold);
    return {
      kind: 'success',
      taskId: run.taskId,
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
      taskId: run.taskId,
      message: msg,
      failKind,
      seed: run.seed,
    };
  }
}

function handleDrop(drop: DropSession): void {
  sessions.delete(drop.sessionId);
}

(self as any).onmessage = async (event: MessageEvent<WorkerOutbound>) => {
  const msg = event.data;
  if (msg.kind === 'setup-fit') {
    const ack = handleSetup(msg);
    (self as any).postMessage(ack);
    return;
  }
  if (msg.kind === 'run-seed') {
    const reply = await handleRun(msg);
    const transferables: Transferable[] = [];
    if (reply.kind === 'success') transferables.push(reply.point.buffer);
    else transferables.push(reply.seed.buffer);
    (self as any).postMessage(reply, transferables);
    return;
  }
  if (msg.kind === 'drop-session') {
    handleDrop(msg);
    return;
  }
};

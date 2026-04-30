// Fitting worker entry. Receives `FitSessionSetup` once per fit, then
// `RunDispatch` per seed; `DropSession` releases per-session state at teardown.
// All imports are DG-free — anything that touches DG.* would crash a worker
// the moment the bundle loaded.

// Webpack externalizes the bare 'dayjs' specifier to the platform global by
// exact match, so we import via the deep ESM path to force dayjs into the
// worker chunk. Mirroring onto globalThis lets user script bodies that write
// `dayjs(...)` resolve as on the main thread.
import dayjs from 'dayjs/esm/index.js';
(globalThis as any).dayjs = dayjs;
import {optimizeNM} from '../optimizer-nelder-mead';
import {LOSS} from '../constants';
import {makeBoundsChecker} from '../bounds-checker';
import {arrowIpcToLite} from '../../../../webworkers/dg-lite/arrow-to-lite';
import type {LiteDataFrame} from '../../../../webworkers/dg-lite/types';
import {ColLike, DfLike, FrameTarget, ScalarTarget, accumulateLoss, InconsistentTablesError}
  from './cost-math';
import {createWorkerFuncCall, WorkerFuncCall}
  from '../../../../webworkers/script-runner/func-call-shim';
import type {
  FitSessionSetup,
  RunDispatch,
  DropSession,
  WorkerOutbound,
  SetupAck,
  WorkerSuccess,
  WorkerFailure,
  SerializedOutputTarget,
  SerializedDataFrameTarget,
  SessionId,
  FixedInputKind,
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
  for (const [name, bytes] of Object.entries(blobs)) {
    const df = arrowIpcToLite(bytes);
    df.newId();
    out[name] = df;
  }
  return out;
}

// Mirror of serialize.serializeFixedInputs: turn tagged ISO strings back
// into Dayjs / Date so formulas and script bodies see the same shapes the
// main arm does.
function reifyFixedInputs(
  scalars: Record<string, any>,
  types: Record<string, FixedInputKind>,
): Record<string, any> {
  const out: Record<string, any> = {};
  for (const [name, value] of Object.entries(scalars)) {
    const kind = types[name];
    if (kind === 'dayjs') out[name] = dayjs(value as string);
    else if (kind === 'date') out[name] = new Date(value as string);
    else out[name] = value;
  }
  return out;
}

type Session = {
  costFunc: (x: Float64Array) => number | undefined;
  nmSettings: Map<string, number>;
  threshold?: number;
  // Strong reference to the funcCall so V8 keeps the compiled body and
  // decoded LiteDataFrames alive across runs.
  fc: WorkerFuncCall;
};

const sessions: Map<SessionId, Session> = new Map();

function buildCostFunc(setup: FitSessionSetup): {
  cost: (x: Float64Array) => number | undefined;
  fc: WorkerFuncCall;
} {
  const targets = buildTargetEntries(setup.outputTargets);
  // Shared by FuncCall body and bounds-checker formula context.
  const merged: Record<string, any> = {
    ...reifyFixedInputs(setup.fixedInputs, setup.fixedInputTypes),
    ...reifyFixedDataFrames(setup.fixedDataFrames),
  };
  const fc = createWorkerFuncCall({
    source: setup.fnSource,
    paramList: setup.paramList,
    outputNames: setup.outputParamNames,
    fixedInputs: merged,
    variedInputNames: setup.variedInputNames,
  });
  const variedNames = setup.variedInputNames;
  const checker = makeBoundsChecker(setup.bounds, variedNames, merged);
  const useRmse = setup.lossType === LOSS.RMSE;

  const cost = (x: Float64Array): number | undefined => {
    if (!checker(x)) return undefined;
    const varied: Record<string, number> = {};
    for (let i = 0; i < variedNames.length; ++i)
      varied[variedNames[i]] = x[i];
    fc.call(varied);

    const scalars: ScalarTarget[] = targets.scalars.map((s) =>
      ({target: s.target, sim: fc.getParamValue(s.propName) as number}));
    const frames: FrameTarget[] = targets.dataFrames.map((d) =>
      ({argCol: d.argCol, funcCols: d.funcCols, simDf: fc.getParamValue(d.propName) as DfLike}));
    return accumulateLoss(useRmse, scalars, frames);
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

async function handleRun(run: RunDispatch): Promise<WorkerSuccess | WorkerFailure> {
  const session = sessions.get(run.sessionId);
  if (!session) {
    return {
      kind: 'failure',
      taskId: run.taskId,
      seedIndex: run.seedIndex,
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
      seedIndex: run.seedIndex,
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
      seedIndex: run.seedIndex,
      message: msg,
      failKind,
      seed: run.seed,
    };
  }
}

function handleDrop(drop: DropSession): void {
  sessions.delete(drop.sessionId);
}

// compute-utils' tsconfig lib has no webworker, so `self` is typed as Window.
// Cast through a typed shim instead of sprinkling `(self as any)` everywhere.
type OutboundReply = SetupAck | WorkerSuccess | WorkerFailure;

interface FittingWorkerScope {
  onmessage: ((ev: MessageEvent<WorkerOutbound>) => unknown) | null;
  postMessage(message: OutboundReply, transfer?: Transferable[]): void;
}

const ctx = self as unknown as FittingWorkerScope;

// Convert a structured-clone failure on the reply into a non-transferable
// failure reply, so the parent's dispatchRun doesn't hang on a worker
// `unhandledrejection` that wouldn't reliably reach `worker.onerror`.
function safePostMessage(msg: OutboundReply, transferables?: Transferable[]): void {
  try {
    if (transferables && transferables.length) ctx.postMessage(msg, transferables);
    else ctx.postMessage(msg);
  } catch (e) {
    const taskId = (msg.kind === 'success' || msg.kind === 'failure') ? msg.taskId : 0;
    const seedIndex = (msg.kind === 'success' || msg.kind === 'failure') ? msg.seedIndex : -1;
    const fallback: WorkerFailure = {
      kind: 'failure',
      taskId,
      seedIndex,
      message: `reply send failed: ${e instanceof Error ? e.message : String(e)}`,
      failKind: 'other',
      seed: new Float64Array(0),
    };
    try { ctx.postMessage(fallback); } catch { /* unrecoverable */ }
  }
}

ctx.onmessage = async (event: MessageEvent<WorkerOutbound>) => {
  const msg = event.data;
  if (msg.kind === 'setup-fit') {
    const ack = handleSetup(msg);
    safePostMessage(ack);
    return;
  }
  if (msg.kind === 'run-dispatch') {
    let reply: WorkerSuccess | WorkerFailure;
    try {
      reply = await handleRun(msg);
    } catch (e) {
      // handleRun has its own try/catch; this is belt-and-braces.
      reply = {
        kind: 'failure', taskId: msg.taskId, seedIndex: msg.seedIndex,
        message: e instanceof Error ? e.message : String(e),
        failKind: 'other', seed: msg.seed,
      };
    }
    const transferables: Transferable[] = [];
    if (reply.kind === 'success') transferables.push(reply.point.buffer);
    else transferables.push(reply.seed.buffer);
    safePostMessage(reply, transferables);
    return;
  }
  if (msg.kind === 'drop-session') {
    handleDrop(msg);
    return;
  }
};

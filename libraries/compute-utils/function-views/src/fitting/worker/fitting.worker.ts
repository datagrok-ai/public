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
//   - `dayjs` is pulled in to reify ISO-string-encoded Dayjs inputs back to
//     Dayjs objects on this side, matching what the main arm's FuncCall and
//     formula context see.
//   - `cost-math` is pure TypeScript with no datagrok-api/* runtime imports.
//   - `webworkers/dg-lite/*` and `webworkers/script-runner/*` (the Lite DG
//     shim, Arrow→Lite deserializer, and JS-script body compiler) are also
//     pure TypeScript with no datagrok-api/* runtime imports.
//
// Anything that touches DG.* (fitting-utils.ts, optimizer-misc.ts runtime,
// cost-functions.ts) is intentionally NOT imported here — it would crash a
// worker the moment the bundle loaded.

import dayjs from 'dayjs';
// Mirror `window.dayjs` (the main-thread platform global, see
// packages/LibTests/webpack.config.js externals) so script bodies that read
// dayjs by bare name resolve via worker globalThis instead of ReferenceError.
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
  RunSeed,
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
  for (const [name, bytes] of Object.entries(blobs))
    out[name] = arrowIpcToLite(bytes);
  return out;
}

// Mirror of serialize.serializeFixedInputs's tagging logic: turn ISO
// strings tagged `'dayjs'` back into Dayjs objects and `'date'` back into
// Date objects. Plain scalars pass through. The result is what
// formulas + script bodies see — keeping observable behavior identical
// to the main arm.
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
  // One reified map feeds two consumers (FuncCall body + bounds-checker
  // formula context). Built once so both arms see equivalent values.
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

async function handleRun(run: RunSeed): Promise<WorkerSuccess | WorkerFailure> {
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

// `self` is typed as `Window` here because compute-utils' tsconfig
// lib is ["es2023", "dom"] (no webworker). The worker scope's
// postMessage signature differs from Window's (no targetOrigin, transfer
// is the second arg, not third). Cast through a small typed shim so the
// rest of this file stays type-safe instead of using `(self as any)`.
type OutboundReply = SetupAck | WorkerSuccess | WorkerFailure;

interface FittingWorkerScope {
  onmessage: ((ev: MessageEvent<WorkerOutbound>) => unknown) | null;
  postMessage(message: OutboundReply, transfer?: Transferable[]): void;
}

const ctx = self as unknown as FittingWorkerScope;

// Wrap the outbound postMessage so that a structured-clone failure on the
// reply (detached transferable, exotic value past the type system, etc.)
// turns into a non-transferable failure reply instead of an unhandled
// rejection in this async handler — the parent's pending dispatchRun would
// otherwise hang because worker `unhandledrejection` doesn't reliably
// propagate to `worker.onerror`.
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
  if (msg.kind === 'run-seed') {
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

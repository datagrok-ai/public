// Executor — the dispatch layer optimizer.ts branches into. MainExecutor is
// the in-process seed loop; WorkerExecutor primes a per-fit session in every
// pool slot via `setupAll` and dispatches one `RunDispatch` per starting point.
// Both arms share `EarlyStopTracker`, `buildFailsDataFrame`, and `sampleSeeds`.
//
// canHandle() requires a JS-language DG.Script with `//meta.workerSafe: true`
// in its header — the contract is that the script body won't reach for
// `grok.*` / `ui.*` or any other main-thread global.

import * as DG from 'datagrok-api/dg';
import {Extremum, OptimizationResult, ValueBoundsData, OutputTargetItem}
  from '../optimizer-misc';
import {optimizeNM} from '../optimizer-nelder-mead';
import {ReproSettings, EarlyStoppingSettings, LOSS} from '../constants';
import {buildFailsDataFrame, getInputsData, sampleSeeds} from '../fitting-utils';
import {EarlyStopTracker} from '../early-stop-tracker';
import {WorkerPool, defaultPoolSize, RunReply} from './pool';
import {getSharedFittingPool} from './shared-pool';
import {buildSetup} from './serialize';
import {DeferredProgressIndicator} from './deferred-progress';
import type {SessionId} from './wire-types';

export type ExecutorMode = 'main' | 'worker';
export type ExecutorChoice = 'auto' | ExecutorMode;

export interface ExecutorArgs {
  // Main-arm objective; used when the worker arm can't be picked.
  objectiveFunc: (x: Float64Array) => Promise<number | undefined>;
  inputsBounds: Record<string, ValueBoundsData>;
  samplesCount: number;
  settings: Map<string, number>;
  reproSettings: ReproSettings;
  earlyStoppingSettings: EarlyStoppingSettings;

  // Worker-arm payload. WorkerExecutor disassembles these into a session
  // setup; absent → falls back to main.
  func?: DG.Func;
  outputTargets?: OutputTargetItem[];
  lossType?: LOSS;
}

export interface Executor {
  run(args: ExecutorArgs): Promise<OptimizationResult>;
}

// ---- canHandle ------------------------------------------------------------

// Authors opt in via `//meta.workerSafe: true` in the script header; the
// platform strips `meta.`, so it lands at `func.options['workerSafe']`.
const WORKER_SAFE_OPTION = 'workerSafe';

function isWorkerSafe(func: DG.Func): boolean {
  const v = (func.options as Record<string, unknown> | undefined)?.[WORKER_SAFE_OPTION];
  return v === true || v === 'true';
}

export function canHandle(args: ExecutorArgs): boolean {
  if (typeof navigator === 'undefined' || (navigator.hardwareConcurrency ?? 0) < 2)
    return false;
  if (!(args.func && args.func instanceof DG.Script)) return false;
  const script = args.func as DG.Script;
  if (script.language !== 'javascript') return false;
  if (!isWorkerSafe(script)) return false;
  if (!args.outputTargets || !args.lossType) return false;
  return true;
}

// ---- Main arm -------------------------------------------------------------

export class MainExecutor implements Executor {
  async run(args: ExecutorArgs): Promise<OptimizationResult> {
    const params = sampleSeeds(args.samplesCount, args.inputsBounds, args.reproSettings);

    const tracker = new EarlyStopTracker(args.earlyStoppingSettings);
    const warnings: string[] = [];
    const failedInitPoint: Float64Array[] = [];
    // Pass the threshold into NM so it short-circuits at cost ≤ threshold.
    // Without this, the two arms diverge at the same seed.
    const nmThreshold = args.earlyStoppingSettings.useEarlyStopping ?
      args.earlyStoppingSettings.costFuncThreshold :
      undefined;

    let percentage = 0;
    const pi = DG.TaskBarProgressIndicator.create(`Fitting... (${percentage}%)`, {cancelable: true});

    for (let i = 0; i < args.samplesCount; ++i) {
      try {
        if ((pi as any).canceled) break;

        const extremum = await optimizeNM(args.objectiveFunc, params[i], args.settings, nmThreshold);
        if (tracker.accept(extremum)) break;

        percentage = Math.floor(100 * (i + 1) / args.samplesCount);
        pi.update(percentage, `Fitting... (${percentage}%)`);
      } catch (e) {
        // Match by `name` (not instanceof) — worker/cost-math defines its
        // own InconsistentTablesError; see cost-math.ts for why.
        if (e instanceof Error && e.name === 'InconsistentTables') {
          pi.close();
          throw new Error(`Inconsistent dataframes: ${e.message}`);
        }
        warnings.push((e instanceof Error) ? e.message : 'Platform issue');
        failedInitPoint.push(params[i]);
      }
    }

    pi.close();

    return {extremums: tracker.finalize(), fails: buildFailsDataFrame(failedInitPoint, warnings)};
  }
}

// ---- Worker arm -----------------------------------------------------------

function extractFnSource(func: DG.Func): {body: string; paramList: string[]; outputNames: string[]} {
  const script = func as DG.Script;
  const body = script.script;
  const inputParams = Array.from(func.inputs).map((p: any) => p.name as string);
  const outputParams = Array.from(func.outputs).map((p: any) => p.name as string);
  return {body, paramList: inputParams, outputNames: outputParams};
}

let nextSessionId: SessionId = 1;

export class WorkerExecutor implements Executor {
  constructor(private readonly pool: WorkerPool) {}

  async run(args: ExecutorArgs): Promise<OptimizationResult> {
    const params = sampleSeeds(args.samplesCount, args.inputsBounds, args.reproSettings);

    let percentage = 0;
    const pi = DeferredProgressIndicator.create(`Fitting... (${percentage}%)`, {cancelable: true});

    const settings = args.earlyStoppingSettings;
    const useES = settings.useEarlyStopping;
    const threshold = useES ? settings.costFuncThreshold! : Infinity;
    const stopAfter = settings.stopAfter;

    type SeedResult =
      | {kind: 'pending'}
      | {kind: 'success', extremum: Extremum}
      | {kind: 'failure-warn', seed: Float64Array, message: string}
      | {kind: 'failure-fatal'};
    const PENDING: SeedResult = {kind: 'pending'};
    const results: SeedResult[] = new Array(args.samplesCount).fill(PENDING);

    let inconsistentMessage: string | null = null;
    // Halts new dispatches; in-flight replies still drain via Promise.all.
    let stopRequested = false;
    let validCount = 0;
    let completedCount = 0;

    const sessionId = nextSessionId++;
    const fn = args.func!;
    const {body, paramList, outputNames} = extractFnSource(fn);
    const {variedInputNames, fixedInputs} = getInputsData(args.inputsBounds);
    const built = buildSetup({
      sessionId,
      fnSource: body,
      paramList,
      outputParamNames: outputNames,
      lossType: args.lossType!,
      fixedInputs,
      variedInputNames,
      bounds: args.inputsBounds,
      outputTargets: args.outputTargets!,
      nmSettings: args.settings,
      threshold: useES ? threshold : undefined,
    });

    try {
      try {
        await this.pool.setupAll(built.setup);
      } catch (e: any) {
        pi.close();
        const msg = e instanceof Error ? e.message : String(e);
        if (/InconsistentTables|argument column|output dataframe/i.test(msg))
          throw new Error(`Inconsistent dataframes: ${msg}`);
        throw e;
      }

      // Placement is keyed by `reply.seedIndex` — echoed by the worker and
      // synthesized failures alike, so completion order doesn't matter.
      const handleReply = (reply: RunReply): void => {
        ++completedCount;
        percentage = Math.floor(100 * completedCount / args.samplesCount);
        pi.update(percentage, `Fitting... (${percentage}%)`);

        const idx = reply.seedIndex;
        if (reply.kind === 'failure') {
          if (reply.failKind === 'inconsistent') {
            if (inconsistentMessage == null) inconsistentMessage = reply.message;
            stopRequested = true;
            results[idx] = {kind: 'failure-fatal'};
            return;
          }
          console.warn(`fitting worker seed ${idx} failed: ${reply.message}`);
          results[idx] = {kind: 'failure-warn', seed: reply.seed, message: reply.message};
          return;
        }

        const extremum: Extremum = {
          point: reply.point,
          cost: reply.cost,
          iterCosts: reply.iterCosts,
          iterCount: reply.iterCount,
        };
        results[idx] = {kind: 'success', extremum};
        // Stop on global valid count (not contiguous-prefix) — matches main arm.
        if (useES && extremum.cost <= threshold) {
          ++validCount;
          if (validCount >= stopAfter) stopRequested = true;
        }
      };

      let nextIdx = 0;
      const runSlot = async (): Promise<void> => {
        while (nextIdx < args.samplesCount) {
          if ((pi as any).canceled) {
            stopRequested = true;
            break;
          }
          if (stopRequested) break;
          const idx = nextIdx++;
          const seed = new Float64Array(params[idx]);
          const reply = await this.pool.dispatchRun({sessionId, seedIndex: idx, seed});
          handleReply(reply);
        }
      };

      const slots = Math.min(args.samplesCount, defaultPoolSize());
      await Promise.all(Array.from({length: slots}, () => runSlot()));

      pi.close();

      if (inconsistentMessage)
        throw new Error(`Inconsistent dataframes: ${inconsistentMessage}`);

      // Index-ordered finalize through the shared EarlyStopTracker. Breaking
      // on the first `pending` handles cancellation and stop-with-gap.
      const tracker = new EarlyStopTracker(settings);
      const failedInitPoint: Float64Array[] = [];
      const warnings: string[] = [];

      for (let i = 0; i < args.samplesCount; ++i) {
        const r = results[i];
        if (r.kind === 'pending') break;
        if (r.kind === 'failure-fatal') break;
        if (r.kind === 'failure-warn') {
          warnings.push(r.message);
          failedInitPoint.push(r.seed);
          continue;
        }
        if (tracker.accept(r.extremum)) break;
      }

      return {
        extremums: tracker.finalize(),
        fails: buildFailsDataFrame(failedInitPoint, warnings),
      };
    } finally {
      // Long-lived pools accumulate per-session state otherwise.
      this.pool.dropSession(sessionId);
    }
  }
}

// ---- Pool lifecycle -------------------------------------------------------

// Short-lived pool: created on demand, terminated when the fit returns.
export async function runWithEphemeralPool(args: ExecutorArgs): Promise<OptimizationResult> {
  const pool = new WorkerPool(defaultPoolSize());
  try {
    return await new WorkerExecutor(pool).run(args);
  } finally {
    pool.dispose();
  }
}

// Long-lived pool: reuses the singleton from shared-pool.ts so worker spin-up
// and per-slot compile amortize across fits.
export async function runWithSharedPool(args: ExecutorArgs): Promise<OptimizationResult> {
  const pool = getSharedFittingPool();
  return new WorkerExecutor(pool).run(args);
}

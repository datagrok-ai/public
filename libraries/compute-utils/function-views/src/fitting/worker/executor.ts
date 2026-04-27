// Executor — the dispatch layer optimizer.ts branches into.
//
// MainExecutor wraps the existing seed loop verbatim. WorkerExecutor primes
// a per-fit session in every pool slot via `setupAll`, then dispatches one
// `RunSeed` per starting point so a slow seed doesn't bottleneck the whole
// fit.
//
// canHandle() is the canary that decides whether a fit can run on the worker
// arm at all: requires a JS-language DG.Script whose body doesn't reach for
// `grok.*` or `ui.*` (since neither exists in worker context).

import * as DG from 'datagrok-api/dg';
import {Extremum, OptimizationResult, ValueBoundsData, OutputTargetItem}
  from '../optimizer-misc';
import {optimizeNM} from '../optimizer-nelder-mead';
import {sampleParamsWithFormulaBounds} from '../optimizer-sampler';
import {seededRandom} from '../fitting-utils';
import {ReproSettings, EarlyStoppingSettings, LOSS} from '../constants';
import {getInputsData} from '../fitting-utils';
import {WorkerPool, defaultPoolSize, RunReply} from './pool';
import {buildSetup, buildRunSeed} from './serialize';
import type {SessionId} from './wire-types';

export type ExecutorMode = 'main' | 'worker';
export type ExecutorChoice = 'auto' | ExecutorMode;

export interface ExecutorArgs {
  // Always present — the legacy main-arm objective. Used when the worker
  // arm can't be picked.
  objectiveFunc: (x: Float64Array) => Promise<number | undefined>;
  inputsBounds: Record<string, ValueBoundsData>;
  samplesCount: number;
  settings: Map<string, number>;
  reproSettings: ReproSettings;
  earlyStoppingSettings: EarlyStoppingSettings;

  // Worker-arm payload. When present, WorkerExecutor disassembles this into
  // a session setup. When absent, the worker arm falls back to main.
  func?: DG.Func;
  outputTargets?: OutputTargetItem[];
  lossType?: LOSS;
}

export interface Executor {
  run(args: ExecutorArgs): Promise<OptimizationResult>;
}

// ---- canHandle ------------------------------------------------------------

const FORBIDDEN_API_RE = /(^|[^.\w])(grok|ui)\./;

export function canHandle(args: ExecutorArgs): boolean {
  if (typeof navigator === 'undefined' || (navigator.hardwareConcurrency ?? 0) < 2)
    return false;
  if (!(args.func && args.func instanceof DG.Script)) return false;
  const script = args.func as DG.Script;
  if (script.language !== 'javascript') return false;
  if (FORBIDDEN_API_RE.test(script.script)) return false;
  if (!args.outputTargets || !args.lossType) return false;
  return true;
}

// ---- Main arm -------------------------------------------------------------

export class MainExecutor implements Executor {
  async run(args: ExecutorArgs): Promise<OptimizationResult> {
    const rand = args.reproSettings.reproducible ?
      seededRandom(args.reproSettings.seed) : Math.random;
    const params = sampleParamsWithFormulaBounds(args.samplesCount, args.inputsBounds, rand);

    const validExtremums: Extremum[] = [];
    const allExtremums: Extremum[] = [];
    const warnings: string[] = [];
    const failedInitPoint: Float64Array[] = [];
    let failsCount = 0;
    let failsDF: DG.DataFrame | null = null;

    let percentage = 0;
    const pi = DG.TaskBarProgressIndicator.create(`Fitting... (${percentage}%)`, {cancelable: true});

    const useEarlyStopping = args.earlyStoppingSettings.useEarlyStopping;
    const useAboveThresholdPoints = args.earlyStoppingSettings.useAboveThresholdPoints;
    const threshold = useEarlyStopping ? args.earlyStoppingSettings.costFuncThreshold : undefined;
    const maxValidPoints = args.earlyStoppingSettings.stopAfter;
    let validPointsCount = 0;

    for (let i = 0; i < args.samplesCount; ++i) {
      try {
        if ((pi as any).canceled) break;

        const extremum = await optimizeNM(args.objectiveFunc, params[i], args.settings);

        if (useEarlyStopping) {
          if (extremum.cost <= threshold!) {
            validExtremums.push(extremum);
            ++validPointsCount;
          } else
            allExtremums.push(extremum);
          if (validPointsCount >= maxValidPoints) break;
        } else
          validExtremums.push(extremum);

        percentage = Math.floor(100 * (i + 1) / args.samplesCount);
        pi.update(percentage, `Fitting... (${percentage}%)`);
      } catch (e) {
        // Discriminate by `name` (not instanceof) so both the original
        // optimizer-misc.InconsistentTables and worker/cost-math's
        // InconsistentTablesError are caught — see cost-math.ts for why
        // the worker arm can't share the optimizer-misc class.
        if (e instanceof Error && e.name === 'InconsistentTables') {
          pi.close();
          throw new Error(`Inconsistent dataframes: ${e.message}`);
        }
        ++failsCount;
        warnings.push((e instanceof Error) ? e.message : 'Platform issue');
        failedInitPoint.push(params[i]);
      }
    }

    pi.close();

    if (failsCount > 0) {
      const dim = params[0]?.length ?? 0;
      const raw = new Array<Float64Array>(dim);
      for (let i = 0; i < dim; ++i) raw[i] = new Float64Array(failsCount);
      failedInitPoint.forEach((point, idx) => point.forEach((val, jdx) => raw[jdx][idx] = val));
      failsDF = DG.DataFrame.fromColumns(raw.map((arr, idx) =>
        DG.Column.fromFloat64Array(`arg${idx}`, arr)));
      failsDF.columns.add(DG.Column.fromStrings('Issue', warnings));
    }

    if (useEarlyStopping && useAboveThresholdPoints && maxValidPoints > validExtremums.length)
      validExtremums.push(...allExtremums);

    return {extremums: validExtremums, fails: failsDF};
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
    const rand = args.reproSettings.reproducible ?
      seededRandom(args.reproSettings.seed) : Math.random;
    const params = sampleParamsWithFormulaBounds(args.samplesCount, args.inputsBounds, rand);

    const useEarlyStopping = args.earlyStoppingSettings.useEarlyStopping;
    const useAboveThresholdPoints = args.earlyStoppingSettings.useAboveThresholdPoints;
    const threshold = useEarlyStopping ? args.earlyStoppingSettings.costFuncThreshold : undefined;
    const maxValidPoints = args.earlyStoppingSettings.stopAfter;

    let percentage = 0;
    const pi = DG.TaskBarProgressIndicator.create(`Fitting... (${percentage}%)`, {cancelable: true});

    const validExtremums: Extremum[] = [];
    const allExtremums: Extremum[] = [];
    const warnings: string[] = [];
    const failedInitPoint: Float64Array[] = [];
    let failsCount = 0;
    let failsDF: DG.DataFrame | null = null;
    let validPointsCount = 0;
    let inconsistentMessage: string | null = null;
    let stopRequested = false;
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
      threshold,
    });

    try {
      try {
        await this.pool.setupAll(built.setup, built.transferables);
      } catch (e: any) {
        pi.close();
        const msg = e instanceof Error ? e.message : String(e);
        if (/InconsistentTables|argument column|output dataframe/i.test(msg))
          throw new Error(`Inconsistent dataframes: ${msg}`);
        throw e;
      }

      const handleReply = (reply: RunReply): void => {
        ++completedCount;
        percentage = Math.floor(100 * completedCount / args.samplesCount);
        pi.update(percentage, `Fitting... (${percentage}%)`);

        if (reply.kind === 'failure') {
          if (reply.failKind === 'inconsistent') {
            if (inconsistentMessage == null) inconsistentMessage = reply.message;
            stopRequested = true;
            return;
          }
          ++failsCount;
          warnings.push(reply.message);
          failedInitPoint.push(reply.seed);
          return;
        }

        const extremum: Extremum = {
          point: reply.point,
          cost: reply.cost,
          iterCosts: reply.iterCosts,
          iterCount: reply.iterCount,
        };
        if (!useEarlyStopping) {
          validExtremums.push(extremum);
          return;
        }
        if (extremum.cost <= threshold!) {
          if (validPointsCount >= maxValidPoints) return;
          validExtremums.push(extremum);
          ++validPointsCount;
          if (validPointsCount >= maxValidPoints) stopRequested = true;
        } else
          allExtremums.push(extremum);
      };

      let nextIdx = 0;
      const runSlot = async (): Promise<void> => {
        while (!stopRequested && nextIdx < args.samplesCount) {
          if ((pi as any).canceled) {
            stopRequested = true;
            break;
          }
          const idx = nextIdx++;
          const seed = new Float64Array(params[idx]);
          const dispatched = buildRunSeed({sessionId, taskId: 0, seed});
          const reply = await this.pool.dispatchRun({
            run: dispatched.run, transferables: dispatched.transferables});
          handleReply(reply);
        }
      };

      const slots = Math.min(args.samplesCount, defaultPoolSize());
      await Promise.all(Array.from({length: slots}, () => runSlot()));

      pi.close();

      if (inconsistentMessage)
        throw new Error(`Inconsistent dataframes: ${inconsistentMessage}`);

      if (failsCount > 0) {
        const dim = params[0]?.length ?? 0;
        const raw = new Array<Float64Array>(dim);
        for (let i = 0; i < dim; ++i) raw[i] = new Float64Array(failsCount);
        failedInitPoint.forEach((point, idx) => point.forEach((val, jdx) => raw[jdx][idx] = val));
        failsDF = DG.DataFrame.fromColumns(raw.map((arr, idx) =>
          DG.Column.fromFloat64Array(`arg${idx}`, arr)));
        failsDF.columns.add(DG.Column.fromStrings('Issue', warnings));
      }

      if (useEarlyStopping && useAboveThresholdPoints && maxValidPoints > validExtremums.length)
        validExtremums.push(...allExtremums);

      return {extremums: validExtremums, fails: failsDF};
    } finally {
      // Always release session state in workers, regardless of whether the
      // pool is ephemeral (we're about to dispose it anyway) or long-lived
      // (so per-session state doesn't accumulate across fits).
      this.pool.dropSession(sessionId);
    }
  }
}

// ---- Pool lifecycle (per-call default) ------------------------------------

// Short-lived pool: created on demand, terminated when the fit returns.
// Long-lived callers can construct their own WorkerPool and pass it to
// WorkerExecutor directly; `dropSession` (called inside `run`) keeps
// per-session state from accumulating across fits.
export async function runWithEphemeralPool(args: ExecutorArgs): Promise<OptimizationResult> {
  const pool = new WorkerPool(defaultPoolSize());
  try {
    return await new WorkerExecutor(pool).run(args);
  } finally {
    pool.dispose();
  }
}

// Executor — the dispatch layer optimizer.ts branches into.
//
// MainExecutor wraps the existing seed loop verbatim (preserved from
// optimizer.ts pre-Phase-3c). WorkerExecutor splits the seed loop across a
// pool, building one WorkerTask per starting point so a slow seed doesn't
// bottleneck the whole fit.
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
import {WorkerPool, defaultPoolSize} from './pool';
import {buildFitTask, buildObjectiveTask, FitTask, ObjectiveTask} from './serialize';
import type {WorkerReply} from './wire-types';

export type ExecutorMode = 'main' | 'worker';
export type ExecutorChoice = 'auto' | ExecutorMode;

export interface ExecutorArgs {
  // Always present — the legacy main-arm objective. Used when the worker arm
  // can't be picked.
  objectiveFunc: (x: Float64Array) => Promise<number | undefined>;
  inputsBounds: Record<string, ValueBoundsData>;
  samplesCount: number;
  settings: Map<string, number>;
  reproSettings: ReproSettings;
  earlyStoppingSettings: EarlyStoppingSettings;

  // Worker-arm payload. When present, WorkerExecutor disassembles this into a
  // FitTask. When absent, the worker arm falls back to the main arm (or, for
  // pure-objective tests, uses `objectiveSource` below).
  func?: DG.Func;
  outputTargets?: OutputTargetItem[];
  lossType?: LOSS;

  // Optional — generic-objective worker path. Lets test fixtures send a raw
  // `(x: Float64Array) => number` body to the worker without a DG.Func.
  objectiveSource?: string;
}

export interface Executor {
  run(args: ExecutorArgs): Promise<OptimizationResult>;
}

// ---- canHandle ------------------------------------------------------------

const FORBIDDEN_API_RE = /(^|[^.\w])(grok|ui)\./;

export function canHandle(args: ExecutorArgs): boolean {
  if (typeof navigator === 'undefined' || (navigator.hardwareConcurrency ?? 0) < 2)
    return false;

  // Generic-objective path — no DG.Func, just a JS source string.
  if (args.objectiveSource != null) {
    if (FORBIDDEN_API_RE.test(args.objectiveSource)) return false;
    return true;
  }

  if (!(args.func && args.func instanceof DG.Script)) return false;
  const script = args.func as DG.Script;
  if (script.language !== 'javascript') return false;
  if (FORBIDDEN_API_RE.test(script.script)) return false;
  if (!args.outputTargets || !args.lossType) return false;
  return true;
}

// ---- Main arm -------------------------------------------------------------

// Verbatim port of the seed loop that used to live in optimizer.ts. Kept in
// one place so the worker arm and main arm share the exact same control flow
// for early-stopping, fails-DF assembly, and cancel checks.
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

    // Pre-build the per-task payload bits that don't change per-seed.
    let buildTask: (seed: Float64Array) =>
      {task: ObjectiveTask | FitTask; transferables: Transferable[]};

    if (args.objectiveSource != null) {
      buildTask = (seed) => buildObjectiveTask({
        taskId: 0,
        source: args.objectiveSource!,
        seed,
        nmSettings: args.settings,
        threshold,
      });
    } else {
      const fn = args.func!;
      const {body, paramList, outputNames} = extractFnSource(fn);
      const {variedInputNames, fixedInputs} = getInputsData(args.inputsBounds);
      buildTask = (seed) => buildFitTask({
        taskId: 0,
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
        seed,
      });
    }

    // Per-slot pull-based dispatch: instead of queueing all seeds up front,
    // poolSize "slots" each pull the next available seed, await its reply,
    // and pull the next. Two wins over the upfront-queue version:
    //   1. Early stop honours samplesCount budget — when validPointsCount
    //      hits stopAfter we set stopRequested and slots exit at their next
    //      iteration. In-flight seeds run to completion (single-task overrun)
    //      but no further seeds are dispatched.
    //   2. Progress bar ticks per-reply instead of staying at 0% until the
    //      very last worker finishes.
    let nextIdx = 0;

    const handleReply = (reply: WorkerReply): void => {
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
        // Guard against the in-flight overshoot: if validPointsCount already
        // hit maxValidPoints (because a sibling slot's reply arrived first
        // and tripped stopRequested), drop this extra valid extremum to keep
        // the main-arm invariant validExtremums.length <= maxValidPoints.
        if (validPointsCount >= maxValidPoints) return;
        validExtremums.push(extremum);
        ++validPointsCount;
        if (validPointsCount >= maxValidPoints) stopRequested = true;
      } else
        allExtremums.push(extremum);
    };

    const runSlot = async (): Promise<void> => {
      while (!stopRequested && nextIdx < args.samplesCount) {
        // Same cancellation contract as MainExecutor: poll the indicator at
        // each seed boundary. In-flight NM runs in workers finish on their
        // own (mirroring the main arm's "current await optimizeNM completes
        // before we honour cancel" semantics). Latency is bounded by one
        // NM run — slots run in parallel so wall-clock cancel matches main.
        if ((pi as any).canceled) {
          stopRequested = true;
          break;
        }
        const idx = nextIdx++;
        const seed = new Float64Array(params[idx]);
        const built = buildTask(seed);
        const reply = await this.pool.dispatch({task: built.task, transferables: built.transferables});
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
  }
}

// ---- Pool lifecycle (per-call default) ------------------------------------

// Short-lived pool: created on demand, terminated when the fit returns.
// Reusing a long-lived pool tied to a view's lifecycle would save the
// per-fit worker spin-up cost but is left as a future optimization.
export async function runWithEphemeralPool(args: ExecutorArgs): Promise<OptimizationResult> {
  const pool = new WorkerPool(defaultPoolSize());
  try {
    return await new WorkerExecutor(pool).run(args);
  } finally {
    pool.dispose();
  }
}

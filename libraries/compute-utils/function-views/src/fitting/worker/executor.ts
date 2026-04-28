// Executor — the dispatch layer optimizer.ts branches into.
//
// MainExecutor is the synchronous in-process seed loop. WorkerExecutor
// primes a per-fit session in every pool slot via `setupAll`, then
// dispatches one `RunSeed` per starting point so a slow seed doesn't
// bottleneck the whole fit. Both arms share `EarlyStopTracker` for the
// valid/above-threshold bookkeeping and `buildFailsDataFrame` /
// `sampleSeeds` for the surrounding scaffolding.
//
// canHandle() is the canary that decides whether a fit can run on the worker
// arm at all: requires a JS-language DG.Script whose body doesn't reach for
// `grok.*` or `ui.*` (since neither exists in worker context).

import * as DG from 'datagrok-api/dg';
import {Extremum, OptimizationResult, ValueBoundsData, OutputTargetItem}
  from '../optimizer-misc';
import {optimizeNM} from '../optimizer-nelder-mead';
import {ReproSettings, EarlyStoppingSettings, LOSS} from '../constants';
import {buildFailsDataFrame, getInputsData, sampleSeeds} from '../fitting-utils';
import {EarlyStopTracker} from '../early-stop-tracker';
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
    const params = sampleSeeds(args.samplesCount, args.inputsBounds, args.reproSettings);

    const tracker = new EarlyStopTracker(args.earlyStoppingSettings);
    const warnings: string[] = [];
    const failedInitPoint: Float64Array[] = [];
    // Match WorkerExecutor: when early stopping is on, pass the threshold
    // into each NM run so it short-circuits at cost ≤ threshold. Without
    // this, the main arm refines past the threshold and the two arms
    // produce divergent extremums for the same seed.
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
        // Discriminate by `name` (not instanceof) so both the original
        // optimizer-misc.InconsistentTables and worker/cost-math's
        // InconsistentTablesError are caught — see cost-math.ts for why
        // the worker arm can't share the optimizer-misc class.
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
    const pi = DG.TaskBarProgressIndicator.create(`Fitting... (${percentage}%)`, {cancelable: true});

    const settings = args.earlyStoppingSettings;
    const useES = settings.useEarlyStopping;
    const threshold = useES ? settings.costFuncThreshold! : Infinity;
    const stopAfter = settings.stopAfter;

    // Per-index reply buffer. Replies land in `results[idx]` so finalization
    // can walk in seed-index order (matching MainExecutor) regardless of
    // worker completion order.
    type Slot =
      | {kind: 'pending'}
      | {kind: 'success', extremum: Extremum}
      | {kind: 'failure-warn', seed: Float64Array, message: string}
      | {kind: 'failure-fatal'};
    const PENDING: Slot = {kind: 'pending'};
    const results: Slot[] = new Array(args.samplesCount).fill(PENDING);

    let inconsistentMessage: string | null = null;
    // `stopRequested` halts new dispatches (no more seeds pulled from
    // `nextIdx`). In-flight replies still land in `results[idx]` because
    // each runSlot awaits its dispatchRun before exiting the loop, so
    // `Promise.all(slots)` IS the drain.
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
        await this.pool.setupAll(built.setup, built.transferables);
      } catch (e: any) {
        pi.close();
        const msg = e instanceof Error ? e.message : String(e);
        if (/InconsistentTables|argument column|output dataframe/i.test(msg))
          throw new Error(`Inconsistent dataframes: ${msg}`);
        throw e;
      }

      const handleReply = (idx: number, reply: RunReply): void => {
        ++completedCount;
        percentage = Math.floor(100 * completedCount / args.samplesCount);
        pi.update(percentage, `Fitting... (${percentage}%)`);

        if (reply.kind === 'failure') {
          if (reply.failKind === 'inconsistent') {
            if (inconsistentMessage == null) inconsistentMessage = reply.message;
            stopRequested = true;
            results[idx] = {kind: 'failure-fatal'};
            return;
          }
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
        // Stop dispatching once the global valid count crosses stopAfter.
        // We don't act on contiguous-prefix valid count here — that would
        // dispatch more seeds than necessary. Tracking the global count
        // matches main-arm semantics: once enough valids exist anywhere,
        // the eventual first-M-in-index-order winners are pinned down.
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
          const dispatched = buildRunSeed({sessionId, taskId: 0, seedIndex: idx, seed});
          const reply = await this.pool.dispatchRun({
            run: dispatched.run, transferables: dispatched.transferables});
          handleReply(idx, reply);
        }
      };

      const slots = Math.min(args.samplesCount, defaultPoolSize());
      await Promise.all(Array.from({length: slots}, () => runSlot()));

      pi.close();

      if (inconsistentMessage)
        throw new Error(`Inconsistent dataframes: ${inconsistentMessage}`);

      // Finalize: walk `results` in seed-index order, replaying replies
      // through `EarlyStopTracker` so the selection rule stays shared with
      // MainExecutor. The contiguous-prefix `break` on `pending` handles
      // cancellation and the stop-with-gap case (where lower indices were
      // never dispatched because stopRequested fired before the slots
      // reached them).
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

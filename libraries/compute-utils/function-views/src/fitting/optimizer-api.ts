import * as DG from 'datagrok-api/dg';
import {EarlyStoppingSettings, LOSS, ReproSettings} from './constants';
import {makeConstFunction} from './cost-functions';
import {performNelderMeadOptimization} from './optimizer';
import {OptimizationResult, OptimizerInputsConfig, OptimizerOutputsConfig} from './optimizer-misc';
import {nelderMeadSettingsOpts} from './optimizer-nelder-mead';
import {defaultEarlyStoppingSettings, defaultRandomSeedSettings} from './fitting-utils';
import {finalizeOptimizationResult, FinalizedFitting} from './finalize';
import type {ExecutorChoice} from './worker/executor';


// Public API for Compute2 to expose as a platform function

export type OptimizerParams = {
  lossType: LOSS;
  func: DG.Func;
  inputBounds: OptimizerInputsConfig;
  outputTargets: OptimizerOutputsConfig;
  samplesCount?: number;
  similarity?: number;
  settings?: Map<string, number>;
  reproSettings?: Partial<ReproSettings>;
  earlyStoppingSettings?: Partial<EarlyStoppingSettings>;
  /**
   * Default 'auto' — canHandle routes JS-language scripts annotated with
   * `//meta.workerSafe: true` to the worker arm; everything else falls back
   * to MainExecutor. Pass 'main' or 'worker' to force a specific arm.
   */
  executor?: ExecutorChoice;
};

type ResolvedParams = {
  lossType: LOSS;
  func: DG.Func;
  inputBounds: OptimizerInputsConfig;
  outputTargets: OptimizerOutputsConfig;
  samplesCount: number;
  similarity: number;
  settings: Map<string, number>;
  reproSettings: ReproSettings;
  earlyStoppingSettings: EarlyStoppingSettings;
  executor: ExecutorChoice;
};

function resolveParams(params: OptimizerParams): ResolvedParams {
  const {func, inputBounds, outputTargets, executor = 'auto', reproSettings, earlyStoppingSettings} = params;
  const lossType = params.lossType ?? LOSS.RMSE;

  const defaultsOverrides = JSON.parse(func.options['fittingSettings'] || '{}');

  const settings = params.settings ?? new Map<string, number>();
  nelderMeadSettingsOpts.forEach((opts, key) => {
    if (settings.get(key) == null)
      settings.set(key, defaultsOverrides[key] ?? opts.default);
  });

  const samplesCount = params.samplesCount ?? defaultsOverrides.samplesCount ?? 1;
  const similarity = params.similarity ?? defaultsOverrides.similarity ?? 10;

  const reproSettingsFull: ReproSettings = {...defaultRandomSeedSettings, ...defaultsOverrides, ...(reproSettings ?? {})};
  const earlyStoppingSettingsFull: EarlyStoppingSettings = {...defaultEarlyStoppingSettings,
    ...defaultsOverrides, ...(earlyStoppingSettings ?? {})};

  return {
    lossType, func, inputBounds, outputTargets,
    samplesCount, similarity, settings,
    reproSettings: reproSettingsFull,
    earlyStoppingSettings: earlyStoppingSettingsFull,
    executor,
  };
}

async function runRawOptimizer(p: ResolvedParams): Promise<OptimizationResult> {
  const objectiveFunc = makeConstFunction(p.lossType, p.func, p.inputBounds, p.outputTargets);
  return performNelderMeadOptimization({
    objectiveFunc,
    inputsBounds: p.inputBounds,
    settings: p.settings,
    samplesCount: p.samplesCount,
    reproSettings: p.reproSettings,
    earlyStoppingSettings: p.earlyStoppingSettings,
    executor: p.executor,
    func: p.func,
    outputTargets: p.outputTargets,
    lossType: p.lossType,
  });
}

/**
 * New finalized API. Returns sorted-and-similarity-filtered extrema together
 * with the unfiltered list and the materialized FuncCalls. Single source of
 * truth for post-processing — callers should not re-sort or re-run getNonSimilar.
 */
export async function runOptimizerFinalized(params: OptimizerParams): Promise<FinalizedFitting> {
  const resolved = resolveParams(params);
  const optResult = await runRawOptimizer(resolved);
  return finalizeOptimizationResult(optResult, {
    func: resolved.func,
    inputBounds: resolved.inputBounds,
    outputTargets: resolved.outputTargets,
    similarity: resolved.similarity,
  });
}

/**
 * @deprecated Prefer {@link runOptimizerFinalized}. Kept for backward
 * compatibility with deep-import callers; reimplemented as a thin shim over
 * the finalized API. The tuple shape and `OptimizationResult` shape are unchanged.
 */
export async function runOptimizer(params: OptimizerParams): Promise<[OptimizationResult, DG.FuncCall[]]> {
  const fin = await runOptimizerFinalized(params);
  return [{extremums: fin.allExtremums, fails: fin.fails}, fin.calls];
}

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EarlyStoppingSettings, LOSS, ReproSettings} from './constants';
import {makeConstFunction} from './cost-functions';
import {performNelderMeadOptimization} from './optimizer';
import {Extremum, OptimizationResult, OptimizerInputsConfig, OptimizerOutputsConfig,
  TargetTableOutput} from './optimizer-misc';
import {nelderMeadSettingsOpts} from './optimizer-nelder-mead';
import {defaultEarlyStoppingSettings, defaultRandomSeedSettings, getInputsData,
  makeGetCalledFuncCall} from './fitting-utils';
import {getNonSimilar} from './similarity-utils';


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
};

export async function runOptimizer(
  {
    lossType = LOSS.RMSE,
    func,
    inputBounds,
    outputTargets,
    samplesCount,
    similarity,
    settings,
    reproSettings,
    earlyStoppingSettings,
  }: OptimizerParams,
): Promise<[OptimizationResult, DG.FuncCall[]]> {
  const {variedInputNames, fixedInputs} = getInputsData(inputBounds);

  const objectiveFunc = makeConstFunction(lossType, func, inputBounds, outputTargets);

  const defaultsOverrides = JSON.parse(func.options['fittingSettings'] || '{}');

  settings ??= new Map();
  nelderMeadSettingsOpts.forEach((opts, key) => {
    if (settings.get(key) == null)
      settings.set(key, defaultsOverrides[key] ?? opts.default);
  });

  samplesCount = samplesCount ?? defaultsOverrides.samplesCount ?? 1;
  similarity = similarity ?? defaultsOverrides.similarity ?? 10;

  const reproSettingsFull = {...defaultRandomSeedSettings, ...defaultsOverrides, ...(reproSettings ?? {})};
  const earlyStoppingSettingsFull = {...defaultEarlyStoppingSettings,
    ...defaultsOverrides, ...(earlyStoppingSettings ?? {})};

  const optResult = await performNelderMeadOptimization({
    objectiveFunc,
    inputsBounds: inputBounds,
    settings,
    samplesCount: samplesCount!,
    reproSettings: reproSettingsFull,
    earlyStoppingSettings: earlyStoppingSettingsFull,
  });
  const extrema = optResult.extremums;
  extrema.sort((a: Extremum, b: Extremum) => a.cost - b.cost);

  const getCalledFuncCall = makeGetCalledFuncCall(func, fixedInputs, variedInputNames, true);
  const targetDfs: TargetTableOutput[] = outputTargets
    .filter((output) => output.type === DG.TYPE.DATA_FRAME)
    .map((output) => {
      return {name: output.propName, target: output.target as DG.DataFrame, argColName: output.argName};
    });

  const nonSimilarExtrema = await getNonSimilar(extrema, similarity!, getCalledFuncCall, targetDfs);

  const calls: DG.FuncCall[] = [];
  for (const extrema of nonSimilarExtrema) {
    const calledFuncCall = await getCalledFuncCall(extrema.point);
    calls.push(calledFuncCall);
  }

  return [optResult, calls];
}

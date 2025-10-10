import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EarlyStoppingSettings, LOSS, ReproSettings} from './constants';
import {makeConstFunction} from './cost-functions';
import {performNelderMeadOptimization} from './optimizer';
import {OptimizerInputsConfig, OptimizerOutputsConfig, ValueBoundsData} from './optimizer-misc';
import {nelderMeadSettingsOpts} from './optimizer-nelder-mead';
import {defaultEarlyStoppingSettings, defaultRandomSeedSettings} from './fitting-utils';


function getInputsData(inputBounds: Record<string, ValueBoundsData>) {
  const variedInputNames: string[] = [];
  const fixedInputs: Record<string, any> = {};
  for (const [name, bound] of Object.entries(inputBounds)) {
    if (bound.type === 'const')
      fixedInputs[name] = bound.value;
    else
      variedInputNames.push(name);
  }
  return {variedInputNames, fixedInputs};
}


// Public API for Compute2 to expose as a platform function

export async function runOptimizer(
  {
    lossType = LOSS.RMSE,
    func,
    inputsBounds,
    outputTargets,
    samplesCount = 1,
    settings,
    reproSettings,
    earlyStoppingSettings,
  }: {
    lossType: LOSS;
    func: DG.Func;
    inputsBounds: OptimizerInputsConfig;
    outputTargets: OptimizerOutputsConfig;
    samplesCount: number;
    settings?: Map<string, number>;
    reproSettings?: Partial<ReproSettings>;
    earlyStoppingSettings?: Partial<EarlyStoppingSettings>;
  },
) {
  const {variedInputNames, fixedInputs} = getInputsData(inputsBounds);

  const objectiveFunc = makeConstFunction(lossType, func, fixedInputs, variedInputNames, outputTargets);

  const defaultsOverrides = JSON.parse(func.options['fittingSettings'] || '{}');

  settings ??= new Map();
  nelderMeadSettingsOpts.forEach((opts, key) => {
    if (settings.get(key) == null)
      settings.set(key, defaultsOverrides[key] ?? opts.default);
  });

  const reproSettingsFull = {...defaultRandomSeedSettings, ...defaultsOverrides, ...(reproSettings ?? {})};
  const earlyStoppingSettingsFull = {...defaultEarlyStoppingSettings, ...defaultsOverrides, ...(earlyStoppingSettings ?? {})};

  const optResult = await performNelderMeadOptimization({
    objectiveFunc,
    inputsBounds,
    settings,
    samplesCount,
    reproSettings: reproSettingsFull,
    earlyStoppingSettings: earlyStoppingSettingsFull,
  });
  return optResult;
}

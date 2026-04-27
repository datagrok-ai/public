/** the Nelder-Mead optimizer */
import * as DG from 'datagrok-api/dg';

import {OptimizationResult, OutputTargetItem, ValueBoundsData} from './optimizer-misc';
import {LOSS, ReproSettings, EarlyStoppingSettings} from './constants';
import {ExecutorChoice, MainExecutor, canHandle, runWithEphemeralPool}
  from './worker/executor';

export async function performNelderMeadOptimization(
  {
    objectiveFunc,
    inputsBounds,
    samplesCount = 1,
    settings,
    reproSettings,
    earlyStoppingSettings,
    // Default 'main' — worker dispatch is opt-in. Flipping to 'auto' would
    // route any JS-language DG.Script through the worker by default, which
    // can break callers (e.g. DiffStudio's fallback path) whose script
    // bodies depend on platform globals not available in workers.
    executor = 'main',
    func,
    outputTargets,
    lossType,
    objectiveSource,
  }: {
    objectiveFunc: (x: Float64Array) => Promise<number|undefined>;
    inputsBounds: Record<string, ValueBoundsData>;
    samplesCount: number,
    settings: Map<string, number>;
    reproSettings: ReproSettings;
    earlyStoppingSettings: EarlyStoppingSettings;
    executor?: ExecutorChoice;
    func?: DG.Func;
    outputTargets?: OutputTargetItem[];
    lossType?: LOSS;
    objectiveSource?: string;
  },
): Promise<OptimizationResult> {
  const execArgs = {
    objectiveFunc,
    inputsBounds,
    samplesCount,
    settings,
    reproSettings,
    earlyStoppingSettings,
    func,
    outputTargets,
    lossType,
    objectiveSource,
  };

  const useWorker = executor === 'worker' ||
    (executor === 'auto' && canHandle(execArgs));

  if (useWorker)
    return runWithEphemeralPool(execArgs);

  return new MainExecutor().run(execArgs);
}

/** the Nelder-Mead optimizer */
import * as DG from 'datagrok-api/dg';

import {OptimizationResult, OutputTargetItem, ValueBoundsData} from './optimizer-misc';
import {LOSS, ReproSettings, EarlyStoppingSettings} from './constants';
import {ExecutorChoice, MainExecutor, canHandle, runWithSharedPool}
  from './worker/executor';

export async function performNelderMeadOptimization(
  {
    objectiveFunc,
    inputsBounds,
    samplesCount = 1,
    settings,
    reproSettings,
    earlyStoppingSettings,
    // Default 'auto' — canHandle decides per call. A script is routed to a
    // worker only when it carries `//meta.workerSafe: true` (the explicit
    // author opt-in), is JS-language, has outputTargets + lossType, and the
    // runtime has hardwareConcurrency ≥ 2. Anything else falls back to
    // MainExecutor — non-annotated scripts, non-Script Funcs, callers with
    // a closure-only objectiveFunc — so the UI and RunOptimizer paths get
    // worker speedup for free on annotated bodies without breaking anything
    // that isn't.
    executor = 'auto',
    func,
    outputTargets,
    lossType,
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
  };

  const useWorker = executor === 'worker' ||
    (executor === 'auto' && canHandle(execArgs));

  if (useWorker)
    return runWithSharedPool(execArgs);

  return new MainExecutor().run(execArgs);
}

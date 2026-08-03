// Dual-executor harness for fitting tests.
//
// Each fitting test is run once per executor in `ALL_EXECUTORS`.
// `runFitting('worker', args)` dispatches to the WorkerExecutor when the
// args carry a `DG.Func` body, and falls back to main otherwise so tests
// that only have a closure-based `objectiveFunc` keep working.

import * as DG from 'datagrok-api/dg';
import {performNelderMeadOptimization, LOSS} from './imports';
import type {OptimizationResult, ValueBoundsData,
  EarlyStoppingSettings, ReproSettings, OutputTargetItem} from './imports';

export type Executor = 'main' | 'worker';

export const ALL_EXECUTORS: Executor[] = ['main', 'worker'];

export type FittingArgs = {
  objectiveFunc: (x: Float64Array) => Promise<number | undefined>;
  inputsBounds: Record<string, ValueBoundsData>;
  samplesCount: number;
  settings: Map<string, number>;
  reproSettings: ReproSettings;
  earlyStoppingSettings: EarlyStoppingSettings;
  // Worker-arm payload — opt-in. When absent on the worker arm, the harness
  // dispatches to main so closure-only tests don't break.
  func?: DG.Func;
  outputTargets?: OutputTargetItem[];
  lossType?: LOSS;
};

export async function runFitting(executor: Executor, args: FittingArgs): Promise<OptimizationResult> {
  const wantsWorker = executor === 'worker' && args.func != null;
  return performNelderMeadOptimization({
    objectiveFunc: args.objectiveFunc,
    inputsBounds: args.inputsBounds,
    samplesCount: args.samplesCount,
    settings: args.settings,
    reproSettings: args.reproSettings,
    earlyStoppingSettings: args.earlyStoppingSettings,
    executor: wantsWorker ? 'worker' : 'main',
    func: args.func,
    outputTargets: args.outputTargets,
    lossType: args.lossType,
  });
}

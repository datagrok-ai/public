// Dual-executor harness for fitting tests.
//
// Each fitting test is run once per executor in `ALL_EXECUTORS`. While the
// worker path doesn't exist yet, `runFitting` always dispatches to the legacy
// main-thread implementation; once `WorkerExecutor` lands, swap the dispatch
// here and the same suite enforces the migration invariant.

import {performNelderMeadOptimization} from './imports';
import type {OptimizationResult, ValueBoundsData,
  EarlyStoppingSettings, ReproSettings} from './imports';

export type Executor = 'main' | 'worker';

export const ALL_EXECUTORS: Executor[] = ['main'];

export type FittingArgs = {
  objectiveFunc: (x: Float64Array) => Promise<number | undefined>;
  inputsBounds: Record<string, ValueBoundsData>;
  samplesCount: number;
  settings: Map<string, number>;
  reproSettings: ReproSettings;
  earlyStoppingSettings: EarlyStoppingSettings;
};

export async function runFitting(_executor: Executor, args: FittingArgs): Promise<OptimizationResult> {
  // TODO(workers-phase-3): branch on `_executor` once WorkerExecutor lands.
  return performNelderMeadOptimization(args);
}

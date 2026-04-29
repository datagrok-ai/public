// Single point of coupling to compute-utils' fitting internals. The library
// doesn't re-export these from its public index, so tests reach in via deep
// paths. When files move, only this barrel needs to change.

export {performNelderMeadOptimization} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer';
export {nelderMeadSettingsOpts} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-nelder-mead';
export {makeBoundsChecker} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-sampler';
export type {Extremum, OptimizationResult, OptimizerInputsConfig, OptimizerOutputsConfig,
  OutputTargetItem, ValueBoundsData} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-misc';
export type {EarlyStoppingSettings, ReproSettings} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/constants';
export {LOSS} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/constants';
export {makeConstFunction} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/cost-functions';
export {runOptimizer} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/optimizer-api';
export type {ExecutorChoice, ExecutorArgs} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/worker/executor';
export {canHandle, WorkerExecutor, runWithEphemeralPool} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/worker/executor';
export {buildSetup, buildRunSeed} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/worker/serialize';
export type {FitSessionSetup, FixedInputKind} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/worker/wire-types';
export {WorkerPool} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/worker/pool';
export type {WorkerPoolOptions, RunReply} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/worker/pool';

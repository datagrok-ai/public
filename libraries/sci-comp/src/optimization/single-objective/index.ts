// Types
export type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
  IterationState,
  OptimizationResult,
  IterationCallback,
  CommonSettings,
  Constraint,
  PenaltyMethod,
  PenaltyOptions,
} from './types';

// Base class (for people writing their own optimizers)
export {Optimizer} from './optimizer';

// Registry
export {registerOptimizer, getOptimizer, listOptimizers} from './registry';

// Penalty / constraints layer
export {applyPenalty, applyPenaltyAsync, boxConstraints} from './penalty';

// Built-in optimizers
export {NelderMead} from './optimizers/nelder-mead';
export type {NelderMeadSettings} from './optimizers/nelder-mead';

export {PSO} from './optimizers/pso';
export type {PSOSettings} from './optimizers/pso';

// --- Auto-register built-in optimizers ---
import {registerOptimizer} from './registry';
import {NelderMead} from './optimizers/nelder-mead';
import {PSO} from './optimizers/pso';

registerOptimizer('nelder-mead', () => new NelderMead());
registerOptimizer('pso', () => new PSO());

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

export {GradientDescent} from './optimizers/gradient-descent';
export type {GradientDescentSettings} from './optimizers/gradient-descent';

export {Adam} from './optimizers/adam';
export type {AdamSettings} from './optimizers/adam';

export {LBFGS} from './optimizers/lbfgs';
export type {LBFGSSettings} from './optimizers/lbfgs';

export {LBFGSB} from './optimizers/lbfgs-b';
export type {LBFGSBSettings, LBFGSBLineSearchSettings, LBFGSBBounds} from './optimizers/lbfgs-b';

// --- Auto-register built-in optimizers ---
import {registerOptimizer} from './registry';
import {NelderMead} from './optimizers/nelder-mead';
import {PSO} from './optimizers/pso';
import {GradientDescent} from './optimizers/gradient-descent';
import {Adam} from './optimizers/adam';
import {LBFGS} from './optimizers/lbfgs';
import {LBFGSB} from './optimizers/lbfgs-b';

registerOptimizer('nelder-mead', () => new NelderMead());
registerOptimizer('pso', () => new PSO());
registerOptimizer('gradient-descent', () => new GradientDescent());
registerOptimizer('adam', () => new Adam());
registerOptimizer('l-bfgs', () => new LBFGS());
registerOptimizer('l-bfgs-b', () => new LBFGSB());

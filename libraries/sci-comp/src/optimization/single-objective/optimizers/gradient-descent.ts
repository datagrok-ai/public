import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
  OptimizationResult,
  CommonSettings,
} from '../types';
import {runGradientDescentAsync, GradientDescentDeps} from './gradient-descent-driver';
import {runGradientDescentSync} from './gradient-descent-driver-sync';

/* ------------------------------------------------------------------ */
/*  Algorithm-specific settings                                        */
/* ------------------------------------------------------------------ */

export interface GradientDescentSettings extends CommonSettings {
  /** Step size for parameter updates (default: 0.01). */
  learningRate?: number;
  /** Momentum coefficient β ∈ [0, 1) — fraction of previous velocity retained (default: 0). */
  momentum?: number;
  /** Step size h for central finite-difference gradient approximation (default: 1e-7). */
  finiteDiffStep?: number;
  /** How many iterations without improvement before stopping (default: 50). */
  noImprovementLimit?: number;
  /** Maximum gradient norm — gradients are clipped to this magnitude (default: Infinity). */
  maxGradNorm?: number;
  /** Gradient norm convergence threshold — stop when ‖∇f‖ < gradTolerance (default: 1e-10). */
  gradTolerance?: number;
  /** Learning rate schedule strategy (default: 'constant'). */
  lrSchedule?: 'constant' | 'inverse' | 'exponential';
  /** Decay rate parameter for the learning rate schedule (default: 0). */
  lrDecayRate?: number;
}

/* ------------------------------------------------------------------ */
/*  Optimizer                                                          */
/* ------------------------------------------------------------------ */

/**
 * Gradient descent with momentum and finite-difference gradients. Public
 * methods inherited from `Optimizer<S>` delegate via `runInternal` /
 * `runInternalAsync` shims to the top-level functions in
 * `gradient-descent-driver.ts` (async, source of truth) and
 * `gradient-descent-driver-sync.ts` (codegen-generated).
 */
export class GradientDescent extends Optimizer<GradientDescentSettings> {
  constructor() {
    super('GradientDescent');
  }

  protected withDefaults(s: GradientDescentSettings): GradientDescentSettings {
    return {
      maxIterations: s.maxIterations ?? 10_000,
      tolerance: s.tolerance ?? 1e-8,
      learningRate: s.learningRate ?? 0.01,
      momentum: s.momentum ?? 0,
      finiteDiffStep: s.finiteDiffStep ?? 1e-7,
      noImprovementLimit: s.noImprovementLimit ?? 50,
      maxGradNorm: s.maxGradNorm ?? Infinity,
      gradTolerance: s.gradTolerance ?? 1e-10,
      lrSchedule: s.lrSchedule ?? 'constant',
      lrDecayRate: s.lrDecayRate ?? 0,
      onIteration: s.onIteration,
    };
  }

  protected runInternal(
    fn: ObjectiveFunction, x0: Float64Array, s: GradientDescentSettings,
  ): OptimizationResult {
    return runGradientDescentSync(fn, x0, s, this.makeDeps());
  }

  protected runInternalAsync(
    fn: AsyncObjectiveFunction, x0: Float64Array, s: GradientDescentSettings,
  ): Promise<OptimizationResult> {
    return runGradientDescentAsync(fn, x0, s, this.makeDeps());
  }

  private makeDeps(): GradientDescentDeps {
    return {notify: this.notify.bind(this)};
  }
}

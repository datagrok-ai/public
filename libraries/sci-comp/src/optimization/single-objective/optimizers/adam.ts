import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
  OptimizationResult,
  CommonSettings,
} from '../types';
import {runAdamAsync, AdamDeps} from './adam-driver';
import {runAdamSync} from './adam-driver-sync';

/* ------------------------------------------------------------------ */
/*  Algorithm-specific settings                                        */
/* ------------------------------------------------------------------ */

export interface AdamSettings extends CommonSettings {
  /** Base learning rate α (default: 0.001). */
  learningRate?: number;
  /** Exponential decay rate for first moment β₁ (default: 0.9). */
  beta1?: number;
  /** Exponential decay rate for second moment β₂ (default: 0.999). */
  beta2?: number;
  /** Small constant for numerical stability (default: 1e-8). */
  epsilon?: number;
  /** Step size h for central finite-difference gradient (default: 1e-7). */
  finiteDiffStep?: number;
  /** Gradient norm convergence threshold (default: 1e-10). */
  gradTolerance?: number;
  /** How many iterations without improvement before stopping (default: 50). */
  noImprovementLimit?: number;
  /** Maximum gradient norm — clip before update (default: Infinity). */
  maxGradNorm?: number;
}

/* ------------------------------------------------------------------ */
/*  Optimizer                                                          */
/* ------------------------------------------------------------------ */

/**
 * Adam optimizer. Public methods (`minimize`, `minimizeAsync`) inherited
 * from `Optimizer<S>` delegate through `runInternal`/`runInternalAsync`
 * shims to the top-level driver functions in `adam-driver.ts` (async,
 * source of truth) and `adam-driver-sync.ts` (codegen-generated).
 */
export class Adam extends Optimizer<AdamSettings> {
  constructor() {
    super('Adam');
  }

  protected withDefaults(s: AdamSettings): AdamSettings {
    return {
      maxIterations: s.maxIterations ?? 10_000,
      tolerance: s.tolerance ?? 1e-8,
      learningRate: s.learningRate ?? 0.001,
      beta1: s.beta1 ?? 0.9,
      beta2: s.beta2 ?? 0.999,
      epsilon: s.epsilon ?? 1e-8,
      finiteDiffStep: s.finiteDiffStep ?? 1e-7,
      gradTolerance: s.gradTolerance ?? 1e-10,
      noImprovementLimit: s.noImprovementLimit ?? 50,
      maxGradNorm: s.maxGradNorm ?? Infinity,
      onIteration: s.onIteration,
    };
  }

  protected runInternal(fn: ObjectiveFunction, x0: Float64Array, s: AdamSettings): OptimizationResult {
    return runAdamSync(fn, x0, s, this.makeDeps());
  }

  protected runInternalAsync(
    fn: AsyncObjectiveFunction, x0: Float64Array, s: AdamSettings,
  ): Promise<OptimizationResult> {
    return runAdamAsync(fn, x0, s, this.makeDeps());
  }

  private makeDeps(): AdamDeps {
    return {notify: this.notify.bind(this)};
  }
}

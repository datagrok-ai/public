import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
  OptimizationResult,
  CommonSettings,
} from '../types';
import {runNelderMeadAsync, NelderMeadDeps} from './nelder-mead-driver';
import {runNelderMeadSync} from './nelder-mead-driver-sync';

/* ------------------------------------------------------------------ */
/*  Algorithm-specific settings                                        */
/* ------------------------------------------------------------------ */

export interface NelderMeadSettings extends CommonSettings {
  /** Scale used when a component of x0 is zero (default: 0.00025). */
  nonZeroParam?: number;
  /** Multiplier for building the initial simplex (default: 0.05). */
  initialScale?: number;
  /** Reflection coefficient α (default: 1). */
  reflection?: number;
  /** Expansion coefficient γ (default: 2). */
  expansion?: number;
  /** Contraction coefficient ρ (default: 0.5). */
  contraction?: number;
  /** Shrink coefficient σ (default: 0.5). */
  shrink?: number;
  /** How many iterations without improvement before stopping (default: 2·dim). */
  noImprovementLimit?: number;
}

/* ------------------------------------------------------------------ */
/*  Optimizer                                                          */
/* ------------------------------------------------------------------ */

/**
 * Nelder-Mead simplex optimizer. Public methods inherited from
 * `Optimizer<S>` delegate via `runInternal` / `runInternalAsync` shims to
 * the top-level functions in `nelder-mead-driver.ts` (async, source of
 * truth) and `nelder-mead-driver-sync.ts` (codegen-generated).
 */
export class NelderMead extends Optimizer<NelderMeadSettings> {
  constructor() {
    super('NelderMead');
  }

  protected withDefaults(s: NelderMeadSettings): NelderMeadSettings {
    return {
      maxIterations: s.maxIterations ?? 1_000,
      tolerance: s.tolerance ?? 1e-8,
      nonZeroParam: s.nonZeroParam ?? 0.00025,
      initialScale: s.initialScale ?? 0.05,
      reflection: s.reflection ?? 1,
      expansion: s.expansion ?? 2,
      contraction: s.contraction ?? 0.5,
      shrink: s.shrink ?? 0.5,
      noImprovementLimit: s.noImprovementLimit,
      onIteration: s.onIteration,
    };
  }

  protected runInternal(
    fn: ObjectiveFunction, x0: Float64Array, s: NelderMeadSettings,
  ): OptimizationResult {
    return runNelderMeadSync(fn, x0, s, this.makeDeps());
  }

  protected runInternalAsync(
    fn: AsyncObjectiveFunction, x0: Float64Array, s: NelderMeadSettings,
  ): Promise<OptimizationResult> {
    return runNelderMeadAsync(fn, x0, s, this.makeDeps());
  }

  private makeDeps(): NelderMeadDeps {
    return {notify: this.notify.bind(this)};
  }
}

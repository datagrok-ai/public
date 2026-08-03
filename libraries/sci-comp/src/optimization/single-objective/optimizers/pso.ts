import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
  OptimizationResult,
  CommonSettings,
} from '../types';
import {runPSOAsync, PSODeps} from './pso-driver';
import {runPSOSync} from './pso-driver-sync';

/* ------------------------------------------------------------------ */
/*  Algorithm-specific settings                                        */
/* ------------------------------------------------------------------ */

export interface PSOSettings extends CommonSettings {
  /** Number of particles in the swarm (default: 30). */
  swarmSize?: number;
  /** Inertia weight ω — controls momentum (default: 0.7298). */
  inertia?: number;
  /** Cognitive coefficient c₁ — pull toward personal best (default: 1.4962). */
  cognitive?: number;
  /** Social coefficient c₂ — pull toward global best (default: 1.4962). */
  social?: number;
  /** Maximum velocity as a fraction of the search range (default: 0.2). */
  maxVelocityFraction?: number;
  /**
   * Search range for particle initialization.
   * NOT an optimization constraint — just tells PSO where to scatter
   * the initial swarm. Particles may leave this range during search.
   * Default: x0 ± 10·|x0| per component.
   */
  searchRange?: { lower: Float64Array; upper: Float64Array };
  /** How many iterations without improvement before stopping (default: 50). */
  noImprovementLimit?: number;
  /** Seed for the PRNG. If omitted, results are non-deterministic. */
  seed?: number;
}

/* ------------------------------------------------------------------ */
/*  Optimizer                                                          */
/* ------------------------------------------------------------------ */

/**
 * Particle Swarm Optimization. Public methods inherited from
 * `Optimizer<S>` delegate via `runInternal` / `runInternalAsync` shims to
 * the top-level functions in `pso-driver.ts` (async, source of truth) and
 * `pso-driver-sync.ts` (codegen-generated).
 */
export class PSO extends Optimizer<PSOSettings> {
  constructor() {
    super('PSO');
  }

  protected withDefaults(s: PSOSettings): PSOSettings {
    return {
      maxIterations: s.maxIterations ?? 1_000,
      tolerance: s.tolerance ?? 1e-8,
      swarmSize: s.swarmSize ?? 30,
      inertia: s.inertia ?? 0.7298,
      cognitive: s.cognitive ?? 1.4962,
      social: s.social ?? 1.4962,
      maxVelocityFraction: s.maxVelocityFraction ?? 0.2,
      noImprovementLimit: s.noImprovementLimit ?? 50,
      searchRange: s.searchRange,
      seed: s.seed,
      onIteration: s.onIteration,
    };
  }

  protected runInternal(
    fn: ObjectiveFunction, x0: Float64Array, s: PSOSettings,
  ): OptimizationResult {
    return runPSOSync(fn, x0, s, this.makeDeps());
  }

  protected runInternalAsync(
    fn: AsyncObjectiveFunction, x0: Float64Array, s: PSOSettings,
  ): Promise<OptimizationResult> {
    return runPSOAsync(fn, x0, s, this.makeDeps());
  }

  private makeDeps(): PSODeps {
    return {notify: this.notify.bind(this)};
  }
}

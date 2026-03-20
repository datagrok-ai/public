import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  OptimizationResult,
  CommonSettings,
} from '../types';

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
/*  Seeded PRNG (mulberry32)                                           */
/* ------------------------------------------------------------------ */

type RNG = () => number;

/** Returns a [0, 1) PRNG seeded with a 32-bit integer. */
function mulberry32(seed: number): RNG {
  let s = seed | 0;
  return () => {
    s = (s + 0x6D2B79F5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* ------------------------------------------------------------------ */
/*  Particle                                                           */
/* ------------------------------------------------------------------ */

interface Particle {
  position: Float64Array;
  velocity: Float64Array;
  value: number;
  bestPosition: Float64Array;
  bestValue: number;
}

/* ------------------------------------------------------------------ */
/*  Optimizer                                                          */
/* ------------------------------------------------------------------ */

export class PSO extends Optimizer<PSOSettings> {
  constructor() {
    super('PSO');
  }

  /* ----- defaults -------------------------------------------------- */

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

  /* ----- core algorithm -------------------------------------------- */

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: PSOSettings,
  ): OptimizationResult {
    const n = x0.length;
    const swarmSize = s.swarmSize!;
    const maxIter = s.maxIterations!;
    const tol = s.tolerance!;
    const w = s.inertia!;
    const c1 = s.cognitive!;
    const c2 = s.social!;
    const rng: RNG = s.seed != null ? mulberry32(s.seed) : Math.random;

    // Search range (only for initialization & velocity limits)
    const rangeLo = s.searchRange?.lower ?? this.defaultLower(x0);
    const rangeHi = s.searchRange?.upper ?? this.defaultUpper(x0);
    const vMax = new Float64Array(n);
    for (let i = 0; i < n; i++)
      vMax[i] = s.maxVelocityFraction! * (rangeHi[i] - rangeLo[i]);

    // --- Initialize swarm ---
    const swarm = this.initSwarm(fn, x0, swarmSize, rangeLo, rangeHi, vMax, rng);

    // Global best
    const gBest = new Float64Array(n);
    gBest.set(swarm[0].bestPosition);
    let gBestVal = swarm[0].bestValue;
    for (const p of swarm) {
      if (p.bestValue < gBestVal) {
        gBestVal = p.bestValue;
        gBest.set(p.bestPosition);
      }
    }

    const costHistory = new Float64Array(maxIter);
    let costLen = 0;
    let iteration = 0;
    let prevBest = Infinity;
    let noImprovement = 0;
    let converged = false;

    while (iteration < maxIter) {
      costHistory[costLen++] = gBestVal;

      // --- Convergence ---
      if (iteration > 0 && prevBest - gBestVal > tol)
        noImprovement = 0;
      else if (iteration > 0) {
        noImprovement++;
        if (noImprovement >= s.noImprovementLimit!) {converged = true; break;}
      }
      prevBest = gBestVal;

      // --- Callback ---
      if (this.notify(s.onIteration, {
        iteration,
        bestValue: gBestVal,
        bestPoint: gBest,
      }))
        break;


      // --- Update particles ---
      for (const p of swarm) {
        for (let i = 0; i < n; i++) {
          const r1 = rng();
          const r2 = rng();

          // Velocity update
          p.velocity[i] =
            w * p.velocity[i] +
            c1 * r1 * (p.bestPosition[i] - p.position[i]) +
            c2 * r2 * (gBest[i] - p.position[i]);

          // Clamp velocity (not position — particles roam freely)
          if (p.velocity[i] > vMax[i]) p.velocity[i] = vMax[i];
          else if (p.velocity[i] < -vMax[i]) p.velocity[i] = -vMax[i];

          // Position update — no clamping, unconstrained
          p.position[i] += p.velocity[i];
        }

        // Evaluate
        p.value = fn(p.position);

        // Update personal best
        if (p.value < p.bestValue) {
          p.bestValue = p.value;
          p.bestPosition.set(p.position);

          // Update global best
          if (p.value < gBestVal) {
            gBestVal = p.value;
            gBest.set(p.position);
          }
        }
      }

      iteration++;
    }

    return {
      point: gBest,
      value: gBestVal,
      iterations: iteration,
      converged,
      costHistory: costHistory.subarray(0, costLen),
    };
  }

  /* ----- private helpers ------------------------------------------- */

  private initSwarm(
    fn: ObjectiveFunction,
    x0: Float64Array,
    size: number,
    lo: Float64Array,
    hi: Float64Array,
    vMax: Float64Array,
    rng: RNG,
  ): Particle[] {
    const n = x0.length;
    const swarm: Particle[] = [];

    // First particle starts at x0
    const p0 = new Float64Array(n);
    p0.set(x0);
    const val0 = fn(p0);
    const bp0 = new Float64Array(n);
    bp0.set(p0);
    swarm.push({
      position: p0,
      velocity: new Float64Array(n),
      value: val0,
      bestPosition: bp0,
      bestValue: val0,
    });

    // Rest — random within search range
    for (let k = 1; k < size; k++) {
      const pos = new Float64Array(n);
      const vel = new Float64Array(n);

      for (let i = 0; i < n; i++) {
        pos[i] = lo[i] + rng() * (hi[i] - lo[i]);
        vel[i] = -vMax[i] + rng() * 2 * vMax[i];
      }

      const val = fn(pos);
      const bestPos = new Float64Array(n);
      bestPos.set(pos);

      swarm.push({
        position: pos,
        velocity: vel,
        value: val,
        bestPosition: bestPos,
        bestValue: val,
      });
    }

    return swarm;
  }

  /** Fallback lower range when no explicit searchRange given. */
  private defaultLower(x0: Float64Array): Float64Array {
    const lo = new Float64Array(x0.length);
    for (let i = 0; i < x0.length; i++)
      lo[i] = x0[i] !== 0 ? x0[i] - 10 * Math.abs(x0[i]) : -10;
    return lo;
  }

  /** Fallback upper range when no explicit searchRange given. */
  private defaultUpper(x0: Float64Array): Float64Array {
    const hi = new Float64Array(x0.length);
    for (let i = 0; i < x0.length; i++)
      hi[i] = x0[i] !== 0 ? x0[i] + 10 * Math.abs(x0[i]) : 10;
    return hi;
  }
}

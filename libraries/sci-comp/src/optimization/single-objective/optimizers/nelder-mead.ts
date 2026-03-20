import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  OptimizationResult,
  CommonSettings,
} from '../types';

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

export class NelderMead extends Optimizer<NelderMeadSettings> {
  constructor() {
    super('NelderMead');
  }

  /* ----- defaults -------------------------------------------------- */

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

  /* ----- core algorithm -------------------------------------------- */

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: NelderMeadSettings,
  ): OptimizationResult {
    const n = x0.length; // parameter dimension
    const nVerts = n + 1; // simplex has n+1 vertices
    const maxIter = s.maxIterations!;
    const tol = s.tolerance!;
    const noImpMax = s.noImprovementLimit ?? 2 * nVerts;

    // --- Build initial simplex ---
    const simplex = this.buildSimplex(fn, x0, s);

    // Index array sorted by objective value
    const idx = Array.from({length: nVerts}, (_, i) => i);
    const worst = () => idx[idx.length - 1];
    const best = () => idx[0];

    const centroid = new Float64Array(n);
    const trial = new Float64Array(n);
    const costHistory: number[] = [];

    let iteration = 0;
    let prevBest = Infinity;
    let noImprovement = 0;
    let converged = false;

    const sortIdx = () =>
      idx.sort((a, b) => simplex[a].value - simplex[b].value);

    sortIdx();

    while (iteration < maxIter) {
      sortIdx();

      const bestVal = simplex[best()].value;
      const secondWorstVal = simplex[idx[idx.length - 2]].value;
      const worstVal = simplex[worst()].value;
      costHistory.push(bestVal);

      // --- Convergence check ---
      if (iteration > 0 && prevBest - bestVal > tol)
        noImprovement = 0;
      else if (iteration > 0) {
        noImprovement++;
        if (noImprovement >= noImpMax) {converged = true; break;}
      }
      prevBest = bestVal;

      // --- Callback / early stop ---
      if (this.notify(s.onIteration, {
        iteration,
        bestValue: bestVal,
        bestPoint: simplex[best()].point,
      }))
        break;


      // --- Centroid (exclude worst vertex) ---
      this.computeCentroid(centroid, simplex, idx, worst(), n);

      // --- Reflection ---
      this.transformPoint(trial, centroid, simplex[worst()].point, s.reflection!, n);
      const reflVal = fn(trial);

      if (reflVal < bestVal) {
        // Reflected is best so far — try expansion
        const savedRefl = Float64Array.from(trial);
        this.transformPoint(trial, centroid, simplex[worst()].point, s.expansion!, n);
        const expVal = fn(trial);

        if (expVal < reflVal)
          this.accept(simplex, worst(), trial, expVal);
        else
          this.accept(simplex, worst(), savedRefl, reflVal);

        iteration++;
        continue;
      }

      if (reflVal < secondWorstVal) {
        // Reflected is better than second-worst — accept reflection
        this.accept(simplex, worst(), trial, reflVal);
        iteration++;
        continue;
      }

      // --- Contraction (reflected is >= second-worst) ---
      if (reflVal < worstVal) {
        // Outside contraction: reflected is between second-worst and worst
        this.transformPoint(trial, centroid, simplex[worst()].point, s.contraction!, n);
        const contrVal = fn(trial);

        if (contrVal <= reflVal) {
          this.accept(simplex, worst(), trial, contrVal);
          iteration++;
          continue;
        }
      } else {
        // Inside contraction: reflected is worse than or equal to worst
        this.transformPoint(trial, centroid, simplex[worst()].point, -s.contraction!, n);
        const contrVal = fn(trial);

        if (contrVal < worstVal) {
          this.accept(simplex, worst(), trial, contrVal);
          iteration++;
          continue;
        }
      }

      // --- Shrink ---
      this.shrinkSimplex(simplex, idx, best(), s.shrink!, fn);
      iteration++;
    }

    sortIdx();

    const finalVal = simplex[best()].value;
    while (costHistory.length < iteration) costHistory.push(finalVal);

    return {
      point: simplex[best()].point,
      value: finalVal,
      iterations: iteration,
      converged,
      costHistory,
    };
  }

  /* ----- private helpers ------------------------------------------- */

  private buildSimplex(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: NelderMeadSettings,
  ): { point: Float64Array; value: number }[] {
    const n = x0.length;
    const simplex: { point: Float64Array; value: number }[] = [];

    const p0 = Float64Array.from(x0);
    simplex.push({point: p0, value: fn(p0)});

    for (let i = 0; i < n; i++) {
      const p = Float64Array.from(x0);
      p[i] = p[i] === 0 ? s.nonZeroParam! : p[i] * (1 + s.initialScale!);
      simplex.push({point: p, value: fn(p)});
    }

    return simplex;
  }

  private computeCentroid(
    out: Float64Array,
    simplex: { point: Float64Array }[],
    idx: number[],
    worstIdx: number,
    n: number,
  ): void {
    out.fill(0);
    for (const j of idx) {
      if (j === worstIdx) continue;
      for (let i = 0; i < n; i++) out[i] += simplex[j].point[i];
    }
    for (let i = 0; i < n; i++) out[i] /= n;
  }

  /** centroid + scale * (centroid - worst) */
  private transformPoint(
    out: Float64Array,
    centroid: Float64Array,
    worstPoint: Float64Array,
    scale: number,
    n: number,
  ): void {
    for (let i = 0; i < n; i++)
      out[i] = centroid[i] + scale * (centroid[i] - worstPoint[i]);
  }

  private accept(
    simplex: { point: Float64Array; value: number }[],
    at: number,
    pt: Float64Array,
    val: number,
  ): void {
    simplex[at].point.set(pt);
    simplex[at].value = val;
  }

  private shrinkSimplex(
    simplex: { point: Float64Array; value: number }[],
    idx: number[],
    bestIdx: number,
    sigma: number,
    fn: ObjectiveFunction,
  ): void {
    const bestPt = simplex[bestIdx].point;
    for (const j of idx) {
      if (j === bestIdx) continue;
      const pt = simplex[j].point;
      for (let i = 0; i < pt.length; i++)
        pt[i] = bestPt[i] + sigma * (pt[i] - bestPt[i]);
      simplex[j].value = fn(pt);
    }
  }
}

import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
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

  /* ================================================================== */
  /*  Synchronous path                                                   */
  /* ================================================================== */

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: NelderMeadSettings,
  ): OptimizationResult {
    const n = x0.length;
    const nVerts = n + 1;
    const maxIter = s.maxIterations!;
    const tol = s.tolerance!;
    const noImpMax = s.noImprovementLimit ?? 2 * nVerts;

    const simplex = this.buildSimplex(fn, x0, s);

    const idx = new Uint32Array(nVerts);
    for (let i = 0; i < nVerts; i++) idx[i] = i;

    const centroid = new Float64Array(n);
    const trial = new Float64Array(n);
    const savedRefl = new Float64Array(n);
    const costHistory = new Float64Array(maxIter);
    let costLen = 0;

    let iteration = 0;
    let prevBest = Infinity;
    let noImprovement = 0;
    let converged = false;

    this.sortIdx(idx, simplex);

    while (iteration < maxIter) {
      const bestIdx = idx[0];
      const worstIdx = idx[nVerts - 1];
      const bestVal = simplex[bestIdx].value;
      const secondWorstVal = simplex[idx[nVerts - 2]].value;
      const worstVal = simplex[worstIdx].value;
      costHistory[costLen++] = bestVal;

      if (iteration > 0 && prevBest - bestVal > tol)
        noImprovement = 0;
      else if (iteration > 0) {
        noImprovement++;
        if (noImprovement >= noImpMax) {converged = true; break;}
      }
      prevBest = bestVal;

      if (this.notify(s.onIteration, {
        iteration,
        bestValue: bestVal,
        bestPoint: simplex[bestIdx].point,
      }))
        break;

      this.computeCentroid(centroid, simplex, idx, nVerts, n);

      // --- Reflection ---
      this.transformPoint(trial, centroid, simplex[worstIdx].point, s.reflection!, n);
      const reflVal = fn(trial);

      if (reflVal < bestVal) {
        savedRefl.set(trial);
        this.transformPoint(trial, centroid, simplex[worstIdx].point, s.expansion!, n);
        const expVal = fn(trial);

        if (expVal < reflVal)
          this.acceptAndInsert(simplex, idx, worstIdx, trial, expVal);
        else
          this.acceptAndInsert(simplex, idx, worstIdx, savedRefl, reflVal);

        iteration++;
        continue;
      }

      if (reflVal < secondWorstVal) {
        this.acceptAndInsert(simplex, idx, worstIdx, trial, reflVal);
        iteration++;
        continue;
      }

      // --- Contraction ---
      if (reflVal < worstVal) {
        this.transformPoint(trial, centroid, simplex[worstIdx].point, s.contraction!, n);
        const contrVal = fn(trial);

        if (contrVal <= reflVal) {
          this.acceptAndInsert(simplex, idx, worstIdx, trial, contrVal);
          iteration++;
          continue;
        }
      } else {
        this.transformPoint(trial, centroid, simplex[worstIdx].point, -s.contraction!, n);
        const contrVal = fn(trial);

        if (contrVal < worstVal) {
          this.acceptAndInsert(simplex, idx, worstIdx, trial, contrVal);
          iteration++;
          continue;
        }
      }

      // --- Shrink ---
      this.shrinkSimplex(simplex, idx, idx[0], s.shrink!, fn);
      this.sortIdx(idx, simplex);
      iteration++;
    }

    const finalBestIdx = idx[0];
    const finalVal = simplex[finalBestIdx].value;
    while (costLen < iteration) costHistory[costLen++] = finalVal;

    return {
      point: simplex[finalBestIdx].point,
      value: finalVal,
      iterations: iteration,
      converged,
      costHistory: costHistory.subarray(0, costLen),
    };
  }

  /* ================================================================== */
  /*  Asynchronous path                                                  */
  /* ================================================================== */

  protected async runInternalAsync(
    fn: AsyncObjectiveFunction,
    x0: Float64Array,
    s: NelderMeadSettings,
  ): Promise<OptimizationResult> {
    const n = x0.length;
    const nVerts = n + 1;
    const maxIter = s.maxIterations!;
    const tol = s.tolerance!;
    const noImpMax = s.noImprovementLimit ?? 2 * nVerts;

    const simplex = await this.buildSimplexAsync(fn, x0, s);

    const idx = new Uint32Array(nVerts);
    for (let i = 0; i < nVerts; i++) idx[i] = i;

    const centroid = new Float64Array(n);
    const trial = new Float64Array(n);
    const savedRefl = new Float64Array(n);
    const costHistory = new Float64Array(maxIter);
    let costLen = 0;

    let iteration = 0;
    let prevBest = Infinity;
    let noImprovement = 0;
    let converged = false;

    this.sortIdx(idx, simplex);

    while (iteration < maxIter) {
      const bestIdx = idx[0];
      const worstIdx = idx[nVerts - 1];
      const bestVal = simplex[bestIdx].value;
      const secondWorstVal = simplex[idx[nVerts - 2]].value;
      const worstVal = simplex[worstIdx].value;
      costHistory[costLen++] = bestVal;

      if (iteration > 0 && prevBest - bestVal > tol)
        noImprovement = 0;
      else if (iteration > 0) {
        noImprovement++;
        if (noImprovement >= noImpMax) {converged = true; break;}
      }
      prevBest = bestVal;

      if (this.notify(s.onIteration, {
        iteration,
        bestValue: bestVal,
        bestPoint: simplex[bestIdx].point,
      }))
        break;

      this.computeCentroid(centroid, simplex, idx, nVerts, n);

      // --- Reflection ---
      this.transformPoint(trial, centroid, simplex[worstIdx].point, s.reflection!, n);
      const reflVal = await fn(trial);

      if (reflVal < bestVal) {
        savedRefl.set(trial);
        this.transformPoint(trial, centroid, simplex[worstIdx].point, s.expansion!, n);
        const expVal = await fn(trial);

        if (expVal < reflVal)
          this.acceptAndInsert(simplex, idx, worstIdx, trial, expVal);
        else
          this.acceptAndInsert(simplex, idx, worstIdx, savedRefl, reflVal);

        iteration++;
        continue;
      }

      if (reflVal < secondWorstVal) {
        this.acceptAndInsert(simplex, idx, worstIdx, trial, reflVal);
        iteration++;
        continue;
      }

      // --- Contraction ---
      if (reflVal < worstVal) {
        this.transformPoint(trial, centroid, simplex[worstIdx].point, s.contraction!, n);
        const contrVal = await fn(trial);

        if (contrVal <= reflVal) {
          this.acceptAndInsert(simplex, idx, worstIdx, trial, contrVal);
          iteration++;
          continue;
        }
      } else {
        this.transformPoint(trial, centroid, simplex[worstIdx].point, -s.contraction!, n);
        const contrVal = await fn(trial);

        if (contrVal < worstVal) {
          this.acceptAndInsert(simplex, idx, worstIdx, trial, contrVal);
          iteration++;
          continue;
        }
      }

      // --- Shrink ---
      await this.shrinkSimplexAsync(simplex, idx, idx[0], s.shrink!, fn);
      this.sortIdx(idx, simplex);
      iteration++;
    }

    const finalBestIdx = idx[0];
    const finalVal = simplex[finalBestIdx].value;
    while (costLen < iteration) costHistory[costLen++] = finalVal;

    return {
      point: simplex[finalBestIdx].point,
      value: finalVal,
      iterations: iteration,
      converged,
      costHistory: costHistory.subarray(0, costLen),
    };
  }

  /* ================================================================== */
  /*  Shared helpers (no fn calls — used by both paths)                  */
  /* ================================================================== */

  /** Full sort of idx by simplex values. Used once at start and after shrink. */
  private sortIdx(
    idx: Uint32Array,
    simplex: { value: number }[],
  ): void {
    idx.sort((a, b) => simplex[a].value - simplex[b].value);
  }

  /** Centroid of all vertices except the last one in idx (worst). */
  private computeCentroid(
    out: Float64Array,
    simplex: { point: Float64Array }[],
    idx: Uint32Array,
    nVerts: number,
    n: number,
  ): void {
    out.fill(0);
    const count = nVerts - 1;
    for (let k = 0; k < count; k++) {
      const pt = simplex[idx[k]].point;
      for (let i = 0; i < n; i++) out[i] += pt[i];
    }
    for (let i = 0; i < n; i++) out[i] /= count;
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

  /**
   * Replace vertex `at` with new point/value, then re-insert into sorted idx.
   * O(n) insertion instead of O(n log n) full sort.
   */
  private acceptAndInsert(
    simplex: { point: Float64Array; value: number }[],
    idx: Uint32Array,
    at: number,
    pt: Float64Array,
    val: number,
  ): void {
    simplex[at].point.set(pt);
    simplex[at].value = val;

    const len = idx.length;
    let pos = len - 1;

    while (pos > 0 && simplex[idx[pos - 1]].value > val) {
      idx[pos] = idx[pos - 1];
      pos--;
    }
    idx[pos] = at;
  }

  /* ================================================================== */
  /*  Sync-only helpers (call fn synchronously)                          */
  /* ================================================================== */

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

  private shrinkSimplex(
    simplex: { point: Float64Array; value: number }[],
    idx: Uint32Array,
    bestIdx: number,
    sigma: number,
    fn: ObjectiveFunction,
  ): void {
    const bestPt = simplex[bestIdx].point;
    const len = idx.length;
    for (let k = 0; k < len; k++) {
      const j = idx[k];
      if (j === bestIdx) continue;
      const pt = simplex[j].point;
      for (let i = 0; i < pt.length; i++)
        pt[i] = bestPt[i] + sigma * (pt[i] - bestPt[i]);
      simplex[j].value = fn(pt);
    }
  }

  /* ================================================================== */
  /*  Async-only helpers (call fn with await)                            */
  /* ================================================================== */

  private async buildSimplexAsync(
    fn: AsyncObjectiveFunction,
    x0: Float64Array,
    s: NelderMeadSettings,
  ): Promise<{ point: Float64Array; value: number }[]> {
    const n = x0.length;
    const simplex: { point: Float64Array; value: number }[] = [];

    const p0 = Float64Array.from(x0);
    simplex.push({point: p0, value: await fn(p0)});

    for (let i = 0; i < n; i++) {
      const p = Float64Array.from(x0);
      p[i] = p[i] === 0 ? s.nonZeroParam! : p[i] * (1 + s.initialScale!);
      simplex.push({point: p, value: await fn(p)});
    }

    return simplex;
  }

  private async shrinkSimplexAsync(
    simplex: { point: Float64Array; value: number }[],
    idx: Uint32Array,
    bestIdx: number,
    sigma: number,
    fn: AsyncObjectiveFunction,
  ): Promise<void> {
    const bestPt = simplex[bestIdx].point;
    const len = idx.length;
    for (let k = 0; k < len; k++) {
      const j = idx[k];
      if (j === bestIdx) continue;
      const pt = simplex[j].point;
      for (let i = 0; i < pt.length; i++)
        pt[i] = bestPt[i] + sigma * (pt[i] - bestPt[i]);
      simplex[j].value = await fn(pt);
    }
  }
}

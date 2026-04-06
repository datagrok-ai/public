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

export class GradientDescent extends Optimizer<GradientDescentSettings> {
  constructor() {
    super('GradientDescent');
  }

  /* ----- defaults -------------------------------------------------- */

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

  /* ================================================================== */
  /*  Synchronous path                                                   */
  /* ================================================================== */

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: GradientDescentSettings,
  ): OptimizationResult {
    const n = x0.length;
    const maxIter = s.maxIterations!;
    const tol = s.tolerance!;
    const lr0 = s.learningRate!;
    const beta = s.momentum!;
    const h = s.finiteDiffStep!;
    const noImpMax = s.noImprovementLimit!;
    const maxGN = s.maxGradNorm!;
    const gradTol = s.gradTolerance!;
    const schedule = s.lrSchedule!;
    const rate = s.lrDecayRate!;

    const x = Float64Array.from(x0);
    const velocity = new Float64Array(n);
    const grad = new Float64Array(n);
    const costHistory = new Float64Array(maxIter);
    let costLen = 0;

    const bestPoint = new Float64Array(n);
    bestPoint.set(x);
    let bestValue = fn(x);
    let prevBest = Infinity;
    let noImprovement = 0;
    let converged = false;
    let iteration = 0;

    while (iteration < maxIter) {
      const fx = (iteration === 0) ? bestValue : fn(x);

      if (fx < bestValue) {
        bestValue = fx;
        bestPoint.set(x);
      }
      costHistory[costLen++] = bestValue;

      // Convergence check: no improvement streak
      if (iteration > 0 && prevBest - bestValue > tol)
        noImprovement = 0;
      else if (iteration > 0) {
        noImprovement++;
        if (noImprovement >= noImpMax) {
          converged = true;
          break;
        }
      }
      prevBest = bestValue;

      if (this.notify(s.onIteration, {
        iteration,
        bestValue,
        bestPoint,
      }))
        break;

      // Compute gradient via central finite differences
      this.computeGradient(fn, x, grad, h, n);

      // Gradient norm + clipping + convergence check
      let gradNorm = 0;
      for (let i = 0; i < n; i++) gradNorm += grad[i] * grad[i];
      gradNorm = Math.sqrt(gradNorm);
      if (gradNorm < gradTol) {
        converged = true;
        break;
      }
      if (gradNorm > maxGN) {
        const scale = maxGN / gradNorm;
        for (let i = 0; i < n; i++) grad[i] *= scale;
      }

      // Update with momentum: v = β·v + lr·∇f, x = x - v
      const lr = this.computeLR(lr0, schedule, rate, iteration);
      for (let i = 0; i < n; i++) {
        velocity[i] = beta * velocity[i] + lr * grad[i];
        x[i] -= velocity[i];
      }

      iteration++;
    }

    return {
      point: bestPoint,
      value: bestValue,
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
    s: GradientDescentSettings,
  ): Promise<OptimizationResult> {
    const n = x0.length;
    const maxIter = s.maxIterations!;
    const tol = s.tolerance!;
    const lr0 = s.learningRate!;
    const beta = s.momentum!;
    const h = s.finiteDiffStep!;
    const noImpMax = s.noImprovementLimit!;
    const maxGN = s.maxGradNorm!;
    const gradTol = s.gradTolerance!;
    const schedule = s.lrSchedule!;
    const rate = s.lrDecayRate!;

    const x = Float64Array.from(x0);
    const velocity = new Float64Array(n);
    const grad = new Float64Array(n);
    const costHistory = new Float64Array(maxIter);
    let costLen = 0;

    const bestPoint = new Float64Array(n);
    bestPoint.set(x);
    let bestValue = await fn(x);
    let prevBest = Infinity;
    let noImprovement = 0;
    let converged = false;
    let iteration = 0;

    while (iteration < maxIter) {
      const fx = (iteration === 0) ? bestValue : await fn(x);

      if (fx < bestValue) {
        bestValue = fx;
        bestPoint.set(x);
      }
      costHistory[costLen++] = bestValue;

      // Convergence check: no improvement streak
      if (iteration > 0 && prevBest - bestValue > tol)
        noImprovement = 0;
      else if (iteration > 0) {
        noImprovement++;
        if (noImprovement >= noImpMax) {
          converged = true;
          break;
        }
      }
      prevBest = bestValue;

      if (this.notify(s.onIteration, {
        iteration,
        bestValue,
        bestPoint,
      }))
        break;

      // Compute gradient via central finite differences
      await this.computeGradientAsync(fn, x, grad, h, n);

      // Gradient norm + clipping + convergence check
      let gradNorm = 0;
      for (let i = 0; i < n; i++) gradNorm += grad[i] * grad[i];
      gradNorm = Math.sqrt(gradNorm);
      if (gradNorm < gradTol) {
        converged = true;
        break;
      }
      if (gradNorm > maxGN) {
        const scale = maxGN / gradNorm;
        for (let i = 0; i < n; i++) grad[i] *= scale;
      }

      // Update with momentum: v = β·v + lr·∇f, x = x - v
      const lr = this.computeLR(lr0, schedule, rate, iteration);
      for (let i = 0; i < n; i++) {
        velocity[i] = beta * velocity[i] + lr * grad[i];
        x[i] -= velocity[i];
      }

      iteration++;
    }

    return {
      point: bestPoint,
      value: bestValue,
      iterations: iteration,
      converged,
      costHistory: costHistory.subarray(0, costLen),
    };
  }

  /* ================================================================== */
  /*  Shared helpers                                                     */
  /* ================================================================== */

  /** Compute effective learning rate for the given iteration. */
  private computeLR(lr0: number, schedule: string, rate: number, iteration: number): number {
    switch (schedule) {
    case 'inverse': return lr0 / (1 + rate * iteration);
    case 'exponential': return lr0 * (1 - rate) ** iteration;
    default: return lr0;
    }
  }

  /** Central finite-difference gradient: ∂f/∂xᵢ ≈ (f(x+hᵢ) − f(x−hᵢ)) / 2h */
  private computeGradient(
    fn: ObjectiveFunction,
    x: Float64Array,
    grad: Float64Array,
    h: number,
    n: number,
  ): void {
    for (let i = 0; i < n; i++) {
      const orig = x[i];
      x[i] = orig + h;
      const fPlus = fn(x);
      x[i] = orig - h;
      const fMinus = fn(x);
      x[i] = orig;
      grad[i] = (fPlus - fMinus) / (2 * h);
    }
  }

  /** Async variant of central finite-difference gradient. */
  private async computeGradientAsync(
    fn: AsyncObjectiveFunction,
    x: Float64Array,
    grad: Float64Array,
    h: number,
    n: number,
  ): Promise<void> {
    for (let i = 0; i < n; i++) {
      const orig = x[i];
      x[i] = orig + h;
      const fPlus = await fn(x);
      x[i] = orig - h;
      const fMinus = await fn(x);
      x[i] = orig;
      grad[i] = (fPlus - fMinus) / (2 * h);
    }
  }
}

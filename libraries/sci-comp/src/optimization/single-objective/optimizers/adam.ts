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

export class Adam extends Optimizer<AdamSettings> {
  constructor() {
    super('Adam');
  }

  /* ----- defaults -------------------------------------------------- */

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

  /* ================================================================== */
  /*  Synchronous path                                                   */
  /* ================================================================== */

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: AdamSettings,
  ): OptimizationResult {
    const n = x0.length;
    const maxIter = s.maxIterations!;
    const tol = s.tolerance!;
    const lr = s.learningRate!;
    const beta1 = s.beta1!;
    const beta2 = s.beta2!;
    const eps = s.epsilon!;
    const h = s.finiteDiffStep!;
    const noImpMax = s.noImprovementLimit!;
    const maxGN = s.maxGradNorm!;
    const gradTol = s.gradTolerance!;

    const x = Float64Array.from(x0);
    const m = new Float64Array(n);
    const v = new Float64Array(n);
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

    let beta1t = 1;
    let beta2t = 1;

    while (iteration < maxIter) {
      const fx = (iteration === 0) ? bestValue : fn(x);

      if (fx < bestValue) {
        bestValue = fx;
        bestPoint.set(x);
      }
      costHistory[costLen++] = bestValue;

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

      this.computeGradient(fn, x, grad, h, n);

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

      beta1t *= beta1;
      beta2t *= beta2;
      const mCorr = 1 / (1 - beta1t);
      const vCorr = 1 / (1 - beta2t);

      for (let i = 0; i < n; i++) {
        m[i] = beta1 * m[i] + (1 - beta1) * grad[i];
        v[i] = beta2 * v[i] + (1 - beta2) * grad[i] * grad[i];
        const mHat = m[i] * mCorr;
        const vHat = v[i] * vCorr;
        x[i] -= lr * mHat / (Math.sqrt(vHat) + eps);
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
    s: AdamSettings,
  ): Promise<OptimizationResult> {
    const n = x0.length;
    const maxIter = s.maxIterations!;
    const tol = s.tolerance!;
    const lr = s.learningRate!;
    const beta1 = s.beta1!;
    const beta2 = s.beta2!;
    const eps = s.epsilon!;
    const h = s.finiteDiffStep!;
    const noImpMax = s.noImprovementLimit!;
    const maxGN = s.maxGradNorm!;
    const gradTol = s.gradTolerance!;

    const x = Float64Array.from(x0);
    const m = new Float64Array(n);
    const v = new Float64Array(n);
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

    let beta1t = 1;
    let beta2t = 1;

    while (iteration < maxIter) {
      const fx = (iteration === 0) ? bestValue : await fn(x);

      if (fx < bestValue) {
        bestValue = fx;
        bestPoint.set(x);
      }
      costHistory[costLen++] = bestValue;

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

      await this.computeGradientAsync(fn, x, grad, h, n);

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

      beta1t *= beta1;
      beta2t *= beta2;
      const mCorr = 1 / (1 - beta1t);
      const vCorr = 1 / (1 - beta2t);

      for (let i = 0; i < n; i++) {
        m[i] = beta1 * m[i] + (1 - beta1) * grad[i];
        v[i] = beta2 * v[i] + (1 - beta2) * grad[i] * grad[i];
        const mHat = m[i] * mCorr;
        const vHat = v[i] * vCorr;
        x[i] -= lr * mHat / (Math.sqrt(vHat) + eps);
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

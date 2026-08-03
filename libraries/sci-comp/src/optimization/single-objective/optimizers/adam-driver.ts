// @async-source: adam-driver-sync.ts
// @codegen-rename: runAdamAsync=runAdamSync
// @codegen-rename: computeGradientAsync=computeGradient

// Adam outer driver. The async version is the source of truth; the sync
// sibling (`adam-driver-sync.ts`) is generated from this file by the
// async-to-sync codegen. The `Adam` class in `adam.ts` is a thin shim that
// delegates to these top-level functions.

import type {OptimizationResult, IterationCallback, IterationState} from '../types';
import type {AdamSettings} from './adam';

/** Stateless helpers the driver borrows from the `Adam` class. */
export interface AdamDeps {
  notify: (cb: IterationCallback | undefined, state: IterationState) => boolean;
}

/** Central finite-difference gradient: ∂f/∂xᵢ ≈ (f(x+hᵢ) − f(x−hᵢ)) / 2h */
async function computeGradientAsync(
  fn: (x: Float64Array) => Promise<number>,
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

export async function runAdamAsync(
  fn: (x: Float64Array) => Promise<number>,
  x0: Float64Array,
  s: AdamSettings,
  deps: AdamDeps,
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

    if (deps.notify(s.onIteration, {
      iteration,
      bestValue,
      bestPoint,
    }))
      break;

    await computeGradientAsync(fn, x, grad, h, n);

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

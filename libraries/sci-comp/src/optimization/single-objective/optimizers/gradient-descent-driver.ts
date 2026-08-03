// @async-source: gradient-descent-driver-sync.ts
// @codegen-rename: runGradientDescentAsync=runGradientDescentSync
// @codegen-rename: computeGradientAsync=computeGradient

// Gradient-descent outer driver. The async version is the source of truth;
// the sync sibling (`gradient-descent-driver-sync.ts`) is generated from
// this file by the async-to-sync codegen. The `GradientDescent` class in
// `gradient-descent.ts` is a thin shim that delegates to these top-level
// functions.

import type {OptimizationResult, IterationCallback, IterationState} from '../types';
import type {GradientDescentSettings} from './gradient-descent';

/** Stateless helpers the driver borrows from the `GradientDescent` class. */
export interface GradientDescentDeps {
  notify: (cb: IterationCallback | undefined, state: IterationState) => boolean;
}

/** Compute effective learning rate for the given iteration. */
export function computeLR(lr0: number, schedule: string, rate: number, iteration: number): number {
  switch (schedule) {
  case 'inverse': return lr0 / (1 + rate * iteration);
  case 'exponential': return lr0 * (1 - rate) ** iteration;
  default: return lr0;
  }
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

export async function runGradientDescentAsync(
  fn: (x: Float64Array) => Promise<number>,
  x0: Float64Array,
  s: GradientDescentSettings,
  deps: GradientDescentDeps,
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

    if (deps.notify(s.onIteration, {
      iteration,
      bestValue,
      bestPoint,
    }))
      break;

    // Compute gradient via central finite differences
    await computeGradientAsync(fn, x, grad, h, n);

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
    const lr = computeLR(lr0, schedule, rate, iteration);
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

/* eslint-disable */
// GENERATED — do not edit by hand.
// Run `npm run update-codegen` to regenerate.
// Source: ./gradient-descent-driver.ts
import {OptimizationResult} from '../types';
import {GradientDescentSettings} from './gradient-descent';
import {GradientDescentDeps, computeLR} from './gradient-descent-driver';

function computeGradient(
  fn: (x: Float64Array) => number,
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

export function runGradientDescentSync(
  fn: (x: Float64Array) => number,
  x0: Float64Array,
  s: GradientDescentSettings,
  deps: GradientDescentDeps,
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

    if (deps.notify(s.onIteration, {
      iteration,
      bestValue,
      bestPoint,
    }))
      break;

    // Compute gradient via central finite differences
    computeGradient(fn, x, grad, h, n);

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

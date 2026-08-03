/* eslint-disable */
// GENERATED — do not edit by hand.
// Run `npm run update-codegen` to regenerate.
// Source: ./nelder-mead-driver.ts
import {OptimizationResult} from '../types';
import {NelderMeadSettings} from './nelder-mead';
import {NelderMeadDeps, SimplexVertex, acceptAndInsert, computeCentroid, sortIdx, transformPoint} from './nelder-mead-driver';

function buildSimplex(
  fn: (x: Float64Array) => number,
  x0: Float64Array,
  s: NelderMeadSettings,
): SimplexVertex[] {
  const n = x0.length;
  const simplex: SimplexVertex[] = [];

  const p0 = Float64Array.from(x0);
  simplex.push({point: p0, value: fn(p0)});

  for (let i = 0; i < n; i++) {
    const p = Float64Array.from(x0);
    p[i] = p[i] === 0 ? s.nonZeroParam! : p[i] * (1 + s.initialScale!);
    simplex.push({point: p, value: fn(p)});
  }

  return simplex;
}

function shrinkSimplex(
  simplex: SimplexVertex[],
  idx: Uint32Array,
  bestIdx: number,
  sigma: number,
  fn: (x: Float64Array) => number,
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

export function runNelderMeadSync(
  fn: (x: Float64Array) => number,
  x0: Float64Array,
  s: NelderMeadSettings,
  deps: NelderMeadDeps,
): OptimizationResult {
  const n = x0.length;
  const nVerts = n + 1;
  const maxIter = s.maxIterations!;
  const tol = s.tolerance!;
  const noImpMax = s.noImprovementLimit ?? 2 * nVerts;

  const simplex = buildSimplex(fn, x0, s);

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

  sortIdx(idx, simplex);

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

    if (deps.notify(s.onIteration, {
      iteration,
      bestValue: bestVal,
      bestPoint: simplex[bestIdx].point,
    }))
      break;

    computeCentroid(centroid, simplex, idx, nVerts, n);

    // --- Reflection ---
    transformPoint(trial, centroid, simplex[worstIdx].point, s.reflection!, n);
    const reflVal = fn(trial);

    if (reflVal < bestVal) {
      savedRefl.set(trial);
      transformPoint(trial, centroid, simplex[worstIdx].point, s.expansion!, n);
      const expVal = fn(trial);

      if (expVal < reflVal)
        acceptAndInsert(simplex, idx, worstIdx, trial, expVal);
      else
        acceptAndInsert(simplex, idx, worstIdx, savedRefl, reflVal);

      iteration++;
      continue;
    }

    if (reflVal < secondWorstVal) {
      acceptAndInsert(simplex, idx, worstIdx, trial, reflVal);
      iteration++;
      continue;
    }

    // --- Contraction ---
    if (reflVal < worstVal) {
      transformPoint(trial, centroid, simplex[worstIdx].point, s.contraction!, n);
      const contrVal = fn(trial);

      if (contrVal <= reflVal) {
        acceptAndInsert(simplex, idx, worstIdx, trial, contrVal);
        iteration++;
        continue;
      }
    } else {
      transformPoint(trial, centroid, simplex[worstIdx].point, -s.contraction!, n);
      const contrVal = fn(trial);

      if (contrVal < worstVal) {
        acceptAndInsert(simplex, idx, worstIdx, trial, contrVal);
        iteration++;
        continue;
      }
    }

    // --- Shrink ---
    shrinkSimplex(simplex, idx, idx[0], s.shrink!, fn);
    sortIdx(idx, simplex);
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

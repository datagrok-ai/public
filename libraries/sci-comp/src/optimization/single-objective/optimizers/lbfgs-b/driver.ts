/**
 * Outer driver for L-BFGS-B.
 *
 * Implements §4 of the specification: Cauchy → subspace → Moré–Thuente
 * line search → curvature-gated memory update → convergence tests, with
 * a single memory-reset retry on line-search failure.
 *
 * Termination maps to `OptimizationResult.converged` as described in
 * §0.4 of the plan: only T1 (projected-gradient) and T2 (relative Δf)
 * set `converged = true`; all other exits (iteration/evaluation caps,
 * line-search abnormal termination, callback abort) return
 * `converged = false` with best-so-far point/value.
 */
import {BOUND_LOWER, BOUND_BOTH, BOUND_UPPER} from './types';
import type {LBFGSBSettings} from './types';
import type {ObjectiveFunction, AsyncObjectiveFunction, OptimizationResult} from '../../types';
import type {IterationCallback} from '../../types';
import {BFGSMat} from './bfgs-mat';
import {normalizeBounds, project, projectedGradient, maxFeasibleStep} from './bounds';
import {makeCauchyWorkspace, cauchyPoint} from './cauchy';
import {makeSubspaceWorkspace, subspaceMin} from './subspace';
import {runLineSearch, runLineSearchAsync} from './line-search';

/** Curvature-gate threshold: reject pair when `sᵀy ≤ CURVATURE_EPS · max(1, yᵀy)`. */
const CURVATURE_EPS = Number.EPSILON;

/* ================================================================== */
/*  Synchronous path                                                   */
/* ================================================================== */

export function runSync(
  fn: ObjectiveFunction,
  x0: Float64Array,
  s: Required<Pick<LBFGSBSettings,
    'maxIterations' | 'tolerance' | 'historySize' | 'gradTolerance' |
    'maxFunctionEvaluations' | 'finiteDiffStep' | 'lineSearch'>>
    & Pick<LBFGSBSettings, 'bounds' | 'gradFn' | 'onIteration'>,
): OptimizationResult {
  const n = x0.length;
  const m = s.historySize;
  const maxIter = s.maxIterations;
  const maxFEval = s.maxFunctionEvaluations;
  const gradTol = s.gradTolerance;
  const ftol = s.tolerance;
  const hFD = s.finiteDiffStep;
  const ls = s.lineSearch;
  const gradFn = s.gradFn;
  const onIter = s.onIteration;

  const {lower, upper, nbd, anyFinite} = normalizeBounds(s.bounds, n);

  // `boxed` mirrors the Fortran v3.0 lnsrlb.f flag: true iff every
  // variable has BOTH a finite lower and a finite upper bound. Used
  // to decide whether the unit-step heuristic for α₀ is appropriate.
  let boxed = anyFinite;
  if (boxed) {
    for (let i = 0; i < n; i++) {
      if (nbd[i] !== BOUND_BOTH) {boxed = false; break;}
    }
  }

  const x = Float64Array.from(x0);
  const g = new Float64Array(n);
  const xNew = new Float64Array(n);
  const gNew = new Float64Array(n);
  const pg = new Float64Array(n);
  const d = new Float64Array(n);
  const stepS = new Float64Array(n);
  const stepY = new Float64Array(n);
  const bestX = new Float64Array(n);

  project(x, lower, upper, nbd, x);

  const mat = new BFGSMat(n, m);
  const cauchy = makeCauchyWorkspace(n, m);
  const subspace = makeSubspaceWorkspace(n, m);

  let fEvalCount = 0;
  const evalF = (p: Float64Array): number => {
    fEvalCount++;
    return fn(p);
  };
  const evalGrad = (p: Float64Array, out: Float64Array): void => {
    if (gradFn)
      gradFn(p, out);
    else {
      // Central finite differences — 2n fn evaluations per gradient.
      for (let i = 0; i < n; i++) {
        const orig = p[i];
        p[i] = orig + hFD;
        const fPlus = evalF(p);
        p[i] = orig - hFD;
        const fMinus = evalF(p);
        p[i] = orig;
        out[i] = (fPlus - fMinus) / (2 * hFD);
      }
    }
  };

  let fCur = evalF(x);
  if (!Number.isFinite(fCur))
    throw new Error('L-BFGS-B: objective returned non-finite at x0');
  evalGrad(x, g);
  for (let i = 0; i < n; i++) {
    if (!Number.isFinite(g[i]))
      throw new Error('L-BFGS-B: gradient contains non-finite at x0');
  }

  let bestValue = fCur;
  bestX.set(x);

  const costHistory = new Float64Array(maxIter + 2);
  let costLen = 0;
  costHistory[costLen++] = bestValue;

  // T1 at start.
  let pgNorm = projectedGradient(x, g, lower, upper, nbd, pg);
  if (pgNorm <= gradTol)
    return finaliseResult(bestX, bestValue, 0, true, costHistory, costLen);

  let iteration = 0;
  let converged = false;
  let retriedMemReset = false;
  let lineSearchSteps = 0;
  let lastStepSize = 0;

  while (iteration < maxIter) {
    if (fEvalCount >= maxFEval) break;

    // ---- (a) Cauchy ------------------------------------------------
    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchy);
    if (!cr.ok) {
      if (mat.col === 0) break;
      mat.reset();
      continue;
    }
    const freeCount = cr.freeCount;

    // ---- (b, c) Subspace min ---------------------------------------
    let endpoint: Float64Array;
    if (freeCount > 0 && mat.col > 0) {
      subspaceMin(x, g, lower, upper, nbd,
        cauchy.xc, cauchy.c, cauchy.freeSet, freeCount, mat, subspace);
      endpoint = subspace.xHat;
    } else
      endpoint = cauchy.xc;


    // ---- (d) Search direction --------------------------------------
    for (let i = 0; i < n; i++) d[i] = endpoint[i] - x[i];
    let slope = dot(g, d);

    if (!(slope < 0)) {
      // Fall back to projected steepest descent.
      for (let i = 0; i < n; i++) {
        let di = -g[i];
        if ((nbd[i] === BOUND_LOWER || nbd[i] === BOUND_BOTH) && x[i] <= lower[i] && di < 0) di = 0;
        if ((nbd[i] === BOUND_UPPER || nbd[i] === BOUND_BOTH) && x[i] >= upper[i] && di > 0) di = 0;
        d[i] = di;
      }
      slope = dot(g, d);
      if (!(slope < 0)) {
        // Even projected gradient gives no descent — at KKT point but
        // pgNorm test didn't catch it. Mark converged and exit.
        converged = true;
        break;
      }
    }

    // ---- (e) Max feasible step -------------------------------------
    const stpMax = maxFeasibleStep(x, d, lower, upper, nbd);
    if (!(stpMax > 0)) {
      converged = true;
      break;
    }

    // ---- (f) Initial α ---------------------------------------------
    // Match Fortran lnsrlb.f: on the very first iteration of an
    // unboxed problem use stp = 1 / ‖d‖₂ (clamped to stpMax) to size
    // the step to the magnitude of the search direction; otherwise
    // start from the natural unit step. Once the BFGS approximation
    // is informative the unit step scales d correctly.
    let alpha0: number;
    if (iteration === 0 && !boxed) {
      const dNorm = Math.sqrt(dot(d, d));
      alpha0 = Math.min(1 / Math.max(dNorm, 1e-16), stpMax);
    } else
      alpha0 = Math.min(1, stpMax);

    // ---- (g) Line search -------------------------------------------
    const lsResult = runLineSearch(
      (alpha) => {
        for (let i = 0; i < n; i++) xNew[i] = x[i] + alpha * d[i];
        const f = evalF(xNew);
        if (!Number.isFinite(f)) return {phi: NaN, phiPrime: NaN};
        evalGrad(xNew, gNew);
        return {phi: f, phiPrime: dot(gNew, d)};
      },
      fCur, slope, alpha0,
      {ftol: ls.ftol!, gtol: ls.gtol!, xtol: ls.xtol!, maxSteps: ls.maxSteps!, stpMax},
    );

    if (!lsResult.ok) {
      if (retriedMemReset || mat.col === 0) break;
      mat.reset();
      retriedMemReset = true;
      continue;
    }
    retriedMemReset = false;

    // ---- (h) Accept step -------------------------------------------
    const stp = lsResult.stp;
    for (let i = 0; i < n; i++) xNew[i] = x[i] + stp * d[i];
    project(xNew, lower, upper, nbd, xNew);

    // Use (f, g) cached from the final line-search evaluation. The
    // projection may have nudged xNew by O(ε|x|); accepting this
    // inconsistency matches scipy / Fortran v3.0 behaviour.
    const fNew = lsResult.phi;

    // ---- (i) Relative Δf and (m) projected-gradient tests ----------
    // Both tests run BEFORE the BFGS update, matching Fortran v3.0
    // mainlb.f. If either passes we accept the step but skip the
    // (now redundant) memory update, keeping nskip / theta in sync
    // with scipy at the converged iteration.
    const denom = Math.max(Math.abs(fCur), Math.abs(fNew), 1);
    const fChange = Math.abs(fCur - fNew) / denom;
    pgNorm = projectedGradient(xNew, gNew, lower, upper, nbd, pg);

    if (fChange <= ftol || pgNorm <= gradTol) {
      x.set(xNew);
      fCur = fNew;
      g.set(gNew);
      if (fCur < bestValue) {
        bestValue = fCur;
        bestX.set(x);
      }
      lastStepSize = stp;
      lineSearchSteps = lsResult.nfev;
      iteration++;
      costHistory[costLen++] = bestValue;
      converged = true;
      break;
    }

    // ---- (j, k) Memory update (curvature-gated) --------------------
    for (let i = 0; i < n; i++) {
      stepS[i] = xNew[i] - x[i];
      stepY[i] = gNew[i] - g[i];
    }
    // negGTs = −gₖᵀ d · stp ≥ 0 along a descent direction. Matches
    // Fortran v3.0 ddum = -gdold * stp in mainlb.f.
    mat.update(stepS, stepY, CURVATURE_EPS, -slope * stp);

    // ---- (l) Advance state -----------------------------------------
    x.set(xNew);
    fCur = fNew;
    g.set(gNew);
    if (fCur < bestValue) {
      bestValue = fCur;
      bestX.set(x);
    }
    lastStepSize = stp;
    lineSearchSteps = lsResult.nfev;
    iteration++;
    costHistory[costLen++] = bestValue;

    // ---- (n) Evaluation limit --------------------------------------
    if (fEvalCount >= maxFEval) break;

    // ---- (o) Callback ----------------------------------------------
    if (fireCallback(onIter, iteration, bestValue, bestX, {
      projGradInfNorm: pgNorm,
      functionEvaluations: fEvalCount,
      stepSize: lastStepSize,
      lineSearchSteps,
      historyCount: mat.col,
      activeBounds: n - freeCount,
    })) break;
  }

  return finaliseResult(bestX, bestValue, iteration, converged, costHistory, costLen);
}

/* ================================================================== */
/*  Asynchronous path                                                  */
/* ================================================================== */

export async function runAsync(
  fn: AsyncObjectiveFunction,
  x0: Float64Array,
  s: Required<Pick<LBFGSBSettings,
    'maxIterations' | 'tolerance' | 'historySize' | 'gradTolerance' |
    'maxFunctionEvaluations' | 'finiteDiffStep' | 'lineSearch'>>
    & Pick<LBFGSBSettings, 'bounds' | 'gradFn' | 'gradFnAsync' | 'onIteration'>,
): Promise<OptimizationResult> {
  const n = x0.length;
  const m = s.historySize;
  const maxIter = s.maxIterations;
  const maxFEval = s.maxFunctionEvaluations;
  const gradTol = s.gradTolerance;
  const ftol = s.tolerance;
  const hFD = s.finiteDiffStep;
  const ls = s.lineSearch;
  const gradFnA = s.gradFnAsync;
  const gradFn = s.gradFn;
  const onIter = s.onIteration;

  const {lower, upper, nbd, anyFinite} = normalizeBounds(s.bounds, n);

  // See note in runSync: `boxed` controls the α₀ heuristic.
  let boxed = anyFinite;
  if (boxed) {
    for (let i = 0; i < n; i++) {
      if (nbd[i] !== BOUND_BOTH) {boxed = false; break;}
    }
  }

  const x = Float64Array.from(x0);
  const g = new Float64Array(n);
  const xNew = new Float64Array(n);
  const gNew = new Float64Array(n);
  const pg = new Float64Array(n);
  const d = new Float64Array(n);
  const stepS = new Float64Array(n);
  const stepY = new Float64Array(n);
  const bestX = new Float64Array(n);

  project(x, lower, upper, nbd, x);

  const mat = new BFGSMat(n, m);
  const cauchy = makeCauchyWorkspace(n, m);
  const subspace = makeSubspaceWorkspace(n, m);

  let fEvalCount = 0;
  const evalF = async (p: Float64Array): Promise<number> => {
    fEvalCount++;
    return fn(p);
  };
  const evalGrad = async (p: Float64Array, out: Float64Array): Promise<void> => {
    if (gradFnA)
      await gradFnA(p, out);
    else if (gradFn)
      gradFn(p, out);
    else {
      for (let i = 0; i < n; i++) {
        const orig = p[i];
        p[i] = orig + hFD;
        const fPlus = await evalF(p);
        p[i] = orig - hFD;
        const fMinus = await evalF(p);
        p[i] = orig;
        out[i] = (fPlus - fMinus) / (2 * hFD);
      }
    }
  };

  let fCur = await evalF(x);
  if (!Number.isFinite(fCur))
    throw new Error('L-BFGS-B: objective returned non-finite at x0');
  await evalGrad(x, g);
  for (let i = 0; i < n; i++) {
    if (!Number.isFinite(g[i]))
      throw new Error('L-BFGS-B: gradient contains non-finite at x0');
  }

  let bestValue = fCur;
  bestX.set(x);

  const costHistory = new Float64Array(maxIter + 2);
  let costLen = 0;
  costHistory[costLen++] = bestValue;

  let pgNorm = projectedGradient(x, g, lower, upper, nbd, pg);
  if (pgNorm <= gradTol)
    return finaliseResult(bestX, bestValue, 0, true, costHistory, costLen);

  let iteration = 0;
  let converged = false;
  let retriedMemReset = false;
  let lineSearchSteps = 0;
  let lastStepSize = 0;

  while (iteration < maxIter) {
    if (fEvalCount >= maxFEval) break;

    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchy);
    if (!cr.ok) {
      if (mat.col === 0) break;
      mat.reset();
      continue;
    }
    const freeCount = cr.freeCount;

    let endpoint: Float64Array;
    if (freeCount > 0 && mat.col > 0) {
      subspaceMin(x, g, lower, upper, nbd,
        cauchy.xc, cauchy.c, cauchy.freeSet, freeCount, mat, subspace);
      endpoint = subspace.xHat;
    } else
      endpoint = cauchy.xc;


    for (let i = 0; i < n; i++) d[i] = endpoint[i] - x[i];
    let slope = dot(g, d);
    if (!(slope < 0)) {
      for (let i = 0; i < n; i++) {
        let di = -g[i];
        if ((nbd[i] === BOUND_LOWER || nbd[i] === BOUND_BOTH) && x[i] <= lower[i] && di < 0) di = 0;
        if ((nbd[i] === BOUND_UPPER || nbd[i] === BOUND_BOTH) && x[i] >= upper[i] && di > 0) di = 0;
        d[i] = di;
      }
      slope = dot(g, d);
      if (!(slope < 0)) {
        converged = true;
        break;
      }
    }

    const stpMax = maxFeasibleStep(x, d, lower, upper, nbd);
    if (!(stpMax > 0)) {
      converged = true;
      break;
    }

    // See note in runSync's (f).
    let alpha0: number;
    if (iteration === 0 && !boxed) {
      const dNorm = Math.sqrt(dot(d, d));
      alpha0 = Math.min(1 / Math.max(dNorm, 1e-16), stpMax);
    } else
      alpha0 = Math.min(1, stpMax);

    const lsResult = await runLineSearchAsync(
      async (alpha) => {
        for (let i = 0; i < n; i++) xNew[i] = x[i] + alpha * d[i];
        const f = await evalF(xNew);
        if (!Number.isFinite(f)) return {phi: NaN, phiPrime: NaN};
        await evalGrad(xNew, gNew);
        return {phi: f, phiPrime: dot(gNew, d)};
      },
      fCur, slope, alpha0,
      {ftol: ls.ftol!, gtol: ls.gtol!, xtol: ls.xtol!, maxSteps: ls.maxSteps!, stpMax},
    );

    if (!lsResult.ok) {
      if (retriedMemReset || mat.col === 0) break;
      mat.reset();
      retriedMemReset = true;
      continue;
    }
    retriedMemReset = false;

    const stp = lsResult.stp;
    for (let i = 0; i < n; i++) xNew[i] = x[i] + stp * d[i];
    project(xNew, lower, upper, nbd, xNew);
    const fNew = lsResult.phi;

    // Convergence tests run before the BFGS update; see runSync's (i)/(m).
    const denom = Math.max(Math.abs(fCur), Math.abs(fNew), 1);
    const fChange = Math.abs(fCur - fNew) / denom;
    pgNorm = projectedGradient(xNew, gNew, lower, upper, nbd, pg);

    if (fChange <= ftol || pgNorm <= gradTol) {
      x.set(xNew);
      fCur = fNew;
      g.set(gNew);
      if (fCur < bestValue) {
        bestValue = fCur;
        bestX.set(x);
      }
      lastStepSize = stp;
      lineSearchSteps = lsResult.nfev;
      iteration++;
      costHistory[costLen++] = bestValue;
      converged = true;
      break;
    }

    for (let i = 0; i < n; i++) {
      stepS[i] = xNew[i] - x[i];
      stepY[i] = gNew[i] - g[i];
    }
    // negGTs = −gₖᵀ d · stp ≥ 0 along a descent direction. Matches
    // Fortran v3.0 ddum = -gdold * stp in mainlb.f.
    mat.update(stepS, stepY, CURVATURE_EPS, -slope * stp);

    x.set(xNew);
    fCur = fNew;
    g.set(gNew);
    if (fCur < bestValue) {
      bestValue = fCur;
      bestX.set(x);
    }
    lastStepSize = stp;
    lineSearchSteps = lsResult.nfev;
    iteration++;
    costHistory[costLen++] = bestValue;

    if (fEvalCount >= maxFEval) break;
    if (fireCallback(onIter, iteration, bestValue, bestX, {
      projGradInfNorm: pgNorm,
      functionEvaluations: fEvalCount,
      stepSize: lastStepSize,
      lineSearchSteps,
      historyCount: mat.col,
      activeBounds: n - freeCount,
    })) break;
  }

  return finaliseResult(bestX, bestValue, iteration, converged, costHistory, costLen);
}

/* ================================================================== */
/*  Helpers                                                            */
/* ================================================================== */

function dot(a: Float64Array, b: Float64Array): number {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}

function finaliseResult(
  bestX: Float64Array,
  bestValue: number,
  iterations: number,
  converged: boolean,
  costHistory: Float64Array,
  costLen: number,
): OptimizationResult {
  return {
    point: bestX,
    value: bestValue,
    iterations,
    converged,
    costHistory: costHistory.subarray(0, costLen),
  };
}

function fireCallback(
  cb: IterationCallback | undefined,
  iteration: number,
  bestValue: number,
  bestPoint: Float64Array,
  extra: Record<string, number>,
): boolean {
  return cb?.({iteration, bestValue, bestPoint, extra}) === true;
}

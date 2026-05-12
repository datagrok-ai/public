/* eslint-disable */
// GENERATED — do not edit by hand.
// Run `npm run update-codegen` to regenerate.
// Source: ./lbfgs-driver.ts
import {OptimizationResult} from '../types';
import {LBFGSSettings} from './lbfgs';
import {CURVATURE_EPS, LBFGSDeps, LineSearchResult, twoLoopRecursion} from './lbfgs-driver';

function lineSearch(
  fn: (x: Float64Array) => number,
  x: Float64Array,
  direction: Float64Array,
  fx: number,
  slope: number,
  c1: number,
  initStep: number,
  maxSteps: number,
  n: number,
  xNew: Float64Array,
  result: LineSearchResult,
): void {
  let alpha = initStep;
  let bestAlpha = 0;
  let accepted = false;
  let fNew = fx;
  let foundFinite = false;

  for (let step = 0; step < maxSteps; step++) {
    for (let i = 0; i < n; i++) xNew[i] = x[i] + alpha * direction[i];
    const f = fn(xNew);

    if (isFinite(f)) {
      foundFinite = true;
      bestAlpha = alpha;
      fNew = f;
      if (f <= fx + c1 * alpha * slope) {
        accepted = true;
        break;
      }
    }
    alpha *= 0.5;
  }

  if (!accepted) {
    if (foundFinite)
      for (let i = 0; i < n; i++) xNew[i] = x[i] + bestAlpha * direction[i];
    else {
      xNew.set(x);
      fNew = fx;
      bestAlpha = 0;
    }
  }

  result.stepSize = bestAlpha;
  result.accepted = accepted;
  result.fNew = fNew;
}

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

export function runLBFGSSync(
  fn: (x: Float64Array) => number,
  x0: Float64Array,
  s: LBFGSSettings,
  deps: LBFGSDeps,
): OptimizationResult {
  const n = x0.length;
  const maxIter = s.maxIterations!;
  const tol = s.tolerance!;
  const m = s.historySize!;
  const h = s.finiteDiffStep!;
  const gradTol = s.gradTolerance!;
  const c1 = s.c1!;
  const maxLS = s.maxLineSearchSteps!;
  const initStep = s.initialStepSize!;

  const x = Float64Array.from(x0);
  const grad = new Float64Array(n);
  const direction = new Float64Array(n);
  const xNew = new Float64Array(n);
  const gradNew = new Float64Array(n);

  const S: Float64Array[] = new Array(m);
  const Y: Float64Array[] = new Array(m);
  for (let i = 0; i < m; i++) {
    S[i] = new Float64Array(n);
    Y[i] = new Float64Array(n);
  }
  const rho = new Float64Array(m);
  const alphaScratch = new Float64Array(m);

  const costHistory = new Float64Array(maxIter + 1);
  let costLen = 0;

  const bestPoint = new Float64Array(n);
  bestPoint.set(x);

  let fCur = fn(x);
  computeGradient(fn, x, grad, h, n);

  let bestValue = fCur;
  let histCount = 0;
  let histHead = 0;
  let iteration = 0;
  let converged = false;
  const lsResult: LineSearchResult = {stepSize: 0, accepted: false, fNew: 0};

  while (iteration < maxIter) {
    // --- Gradient max-norm convergence (matches SciPy's pgtol) ---
    let gNormInf = 0;
    let gNaN = false;
    for (let i = 0; i < n; i++) {
      const g = Math.abs(grad[i]);
      if (g !== g) {gNaN = true; break;}
      if (g > gNormInf) gNormInf = g;
    }
    if (!gNaN && gNormInf < gradTol) {converged = true; break;}

    if (fCur < bestValue) {
      bestValue = fCur;
      bestPoint.set(x);
    }
    costHistory[costLen++] = bestValue;

    if (deps.notify(s.onIteration, {iteration, bestValue, bestPoint}))
      break;

    // --- Search direction via two-loop recursion (steepest descent if no history) ---
    if (histCount === 0)
      for (let i = 0; i < n; i++) direction[i] = -grad[i];
    else {
      twoLoopRecursion(
        grad, S, Y, rho, alphaScratch, histHead, histCount, m, n, direction,
      );
      for (let i = 0; i < n; i++) direction[i] = -direction[i];
    }

    // --- Descent check: if direction is not a descent, reset history ---
    let slope = 0;
    for (let i = 0; i < n; i++) slope += grad[i] * direction[i];
    if (!(slope < 0)) {
      histCount = 0;
      slope = 0;
      for (let i = 0; i < n; i++) {
        direction[i] = -grad[i];
        slope -= grad[i] * grad[i];
      }
    }

    // --- Line search ---
    lineSearch(fn, x, direction, fCur, slope, c1, initStep, maxLS, n, xNew, lsResult);

    // --- Gradient at new point ---
    computeGradient(fn, xNew, gradNew, h, n);

    // --- History update (only on accepted step, gated by curvature condition) ---
    if (lsResult.accepted) {
      const sBuf = S[histHead];
      const yBuf = Y[histHead];
      let ys = 0;
      let yy = 0;
      for (let i = 0; i < n; i++) {
        const si = xNew[i] - x[i];
        const yi = gradNew[i] - grad[i];
        sBuf[i] = si;
        yBuf[i] = yi;
        ys += yi * si;
        yy += yi * yi;
      }
      if (ys > CURVATURE_EPS * Math.max(1, yy)) {
        rho[histHead] = 1 / ys;
        histHead = (histHead + 1) % m;
        if (histCount < m) histCount++;
      }
    }

    // --- Function-value convergence: |Δf| / max(|f_k|, |f_{k+1}|, 1) < tolerance ---
    const fNew = lsResult.fNew;
    const denom = Math.max(Math.abs(fCur), Math.abs(fNew), 1);
    const fChange = Math.abs(fNew - fCur) / denom;

    // --- Advance state ---
    x.set(xNew);
    fCur = fNew;
    grad.set(gradNew);

    if (fCur < bestValue) {
      bestValue = fCur;
      bestPoint.set(x);
    }

    iteration++;

    if (fChange < tol) {converged = true; break;}
  }

  return {
    point: bestPoint,
    value: bestValue,
    iterations: iteration,
    converged,
    costHistory: costHistory.subarray(0, costLen),
  };
}

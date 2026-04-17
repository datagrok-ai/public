/**
 * L-BFGS — Limited-memory Broyden–Fletcher–Goldfarb–Shanno.
 *
 * Primary references:
 *  - J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage"
 *    (1980), Math. Comp., 35(151), 773–782.
 *    Source of the two-loop recursion (see {@link LBFGS.twoLoopRecursion}).
 *  - D. C. Liu, J. Nocedal, "On the limited memory BFGS method for large
 *    scale optimization" (1989), Math. Programming, 45, 503–528.
 *    Source of the initial Hessian scaling γ = (sᵀy)/(yᵀy).
 *  - R. H. Byrd, P. Lu, J. Nocedal, "A Limited Memory Algorithm for Bound
 *    Constrained Optimization" (1995), SIAM J. Sci. Comput., 16(5), 1190–1208.
 *    Basis of the SciPy/Fortran L-BFGS-B code; source of default constants and
 *    the max-norm gradient stopping criterion (pgtol/gtol).
 *
 * Reference implementation:
 *  - SciPy `scipy.optimize._lbfgsb_py._minimize_lbfgsb` (wrapper around
 *    Fortran L-BFGS-B v3.0 by Byrd, Lu, Nocedal & Morales). Defaults
 *    (historySize=10, gradTolerance=1e-5, maxLineSearchSteps=20) and the
 *    function-value relative-change criterion mirror SciPy.
 *
 * Standard-practice components (textbook, not from the primary papers):
 *  - Armijo backtracking line search — Nocedal & Wright, "Numerical
 *    Optimization" (2nd ed., 2006), Ch. 3. A future upgrade to strong-Wolfe
 *    Moré–Thuente (as used in SciPy) is planned.
 *  - Curvature guard `y·s > ε·max(1, y·y)` before history update —
 *    numerical-safety practice, Nocedal & Wright §7.2.
 *
 * Departures from the references:
 *  - No native bound constraints (L-BFGS-**B**): projected gradient / Cauchy
 *    point are not implemented. Constraints go through the penalty layer
 *    (`applyPenalty`).
 *  - No `LbfgsInvHessProduct`, `maxfun`, `factr`, or `iprint` from SciPy.
 */
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

/**
 * Settings for L-BFGS (limited-memory BFGS).
 *
 * L-BFGS is a quasi-Newton method that approximates the inverse Hessian
 * from a small history of past gradient/position changes. Unlike full BFGS
 * (which stores an n×n matrix), L-BFGS stores only `historySize` vector
 * pairs, making it practical for large `n`.
 *
 * Unlike `GradientDescent` / `Adam`, L-BFGS does **not** expose a
 * `noImprovementLimit` streak counter: quasi-Newton methods can
 * legitimately take a non-improving step and still be making progress,
 * so the relative function-value change (`tolerance`) is used instead.
 */
export interface LBFGSSettings extends CommonSettings {
  /** Number of corrections to store (history size m). Higher = better Hessian
   *  approximation but more memory. Typical range: 3–20 (default: 10). */
  historySize?: number;
  /** Step size h for central finite-difference gradient (default: 1e-7). */
  finiteDiffStep?: number;
  /** Gradient convergence threshold — stop when `max|∂f/∂xᵢ| < gradTolerance`
   *  (infinity norm, matches SciPy's `pgtol`/`gtol`). Default: 1e-5. */
  gradTolerance?: number;
  /** Armijo sufficient decrease parameter c₁ ∈ (0, 1) (default: 1e-4). */
  c1?: number;
  /** Maximum number of line search steps per iteration (default: 20). */
  maxLineSearchSteps?: number;
  /** Initial step size for line search (default: 1.0). */
  initialStepSize?: number;
}

/* ------------------------------------------------------------------ */
/*  Internal types                                                     */
/* ------------------------------------------------------------------ */

/** Reusable scratch struct for line search status (no `xNew` — caller owns that buffer). */
interface LineSearchResult {
  stepSize: number;
  accepted: boolean;
  fNew: number;
}

/** Curvature-condition threshold: skip history update when y·s ≤ CURVATURE_EPS·max(1, y·y). */
const CURVATURE_EPS = 1e-10;

/* ------------------------------------------------------------------ */
/*  Optimizer                                                          */
/* ------------------------------------------------------------------ */

/**
 * Limited-memory BFGS optimizer for unconstrained smooth optimization.
 *
 * Uses a fixed-size history of curvature pairs to approximate H·g
 * via the two-loop recursion, combined with Armijo backtracking line search.
 * Constraints are handled through the penalty layer (see `applyPenalty`);
 * native bound constraints (L-BFGS-B) are not implemented.
 */
export class LBFGS extends Optimizer<LBFGSSettings> {
  constructor() {
    super('L-BFGS');
  }

  /* ----- defaults -------------------------------------------------- */

  protected withDefaults(s: LBFGSSettings): LBFGSSettings {
    return {
      maxIterations: s.maxIterations ?? 1_000,
      tolerance: s.tolerance ?? 1e-8,
      historySize: s.historySize ?? 10,
      finiteDiffStep: s.finiteDiffStep ?? 1e-7,
      gradTolerance: s.gradTolerance ?? 1e-5,
      c1: s.c1 ?? 1e-4,
      maxLineSearchSteps: s.maxLineSearchSteps ?? 20,
      initialStepSize: s.initialStepSize ?? 1.0,
      onIteration: s.onIteration,
    };
  }

  /* ================================================================== */
  /*  Synchronous path                                                   */
  /* ================================================================== */

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: LBFGSSettings,
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

    // --- pre-allocated buffers (zero allocation in iteration body) ---
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

    // Evaluate f and ∇f once up-front; reuse across iterations.
    let fCur = fn(x);
    this.computeGradient(fn, x, grad, h, n);

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

      if (this.notify(s.onIteration, {iteration, bestValue, bestPoint}))
        break;

      // --- Search direction via two-loop recursion (steepest descent if no history) ---
      if (histCount === 0)
        for (let i = 0; i < n; i++) direction[i] = -grad[i];
      else {
        this.twoLoopRecursion(
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

      // --- Line search (writes xNew; fills lsResult) ---
      this.lineSearch(fn, x, direction, fCur, slope, c1, initStep, maxLS, n, xNew, lsResult);

      // --- Gradient at new point ---
      this.computeGradient(fn, xNew, gradNew, h, n);

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

  /* ================================================================== */
  /*  Asynchronous path                                                  */
  /* ================================================================== */

  protected async runInternalAsync(
    fn: AsyncObjectiveFunction,
    x0: Float64Array,
    s: LBFGSSettings,
  ): Promise<OptimizationResult> {
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

    let fCur = await fn(x);
    await this.computeGradientAsync(fn, x, grad, h, n);

    let bestValue = fCur;
    let histCount = 0;
    let histHead = 0;
    let iteration = 0;
    let converged = false;
    const lsResult: LineSearchResult = {stepSize: 0, accepted: false, fNew: 0};

    while (iteration < maxIter) {
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

      if (this.notify(s.onIteration, {iteration, bestValue, bestPoint}))
        break;

      if (histCount === 0)
        for (let i = 0; i < n; i++) direction[i] = -grad[i];
      else {
        this.twoLoopRecursion(
          grad, S, Y, rho, alphaScratch, histHead, histCount, m, n, direction,
        );
        for (let i = 0; i < n; i++) direction[i] = -direction[i];
      }

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

      await this.lineSearchAsync(fn, x, direction, fCur, slope, c1, initStep, maxLS, n, xNew, lsResult);

      await this.computeGradientAsync(fn, xNew, gradNew, h, n);

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

      const fNew = lsResult.fNew;
      const denom = Math.max(Math.abs(fCur), Math.abs(fNew), 1);
      const fChange = Math.abs(fNew - fCur) / denom;

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

  /* ================================================================== */
  /*  Private helpers                                                    */
  /* ================================================================== */

  /**
   * Two-loop recursion for L-BFGS inverse-Hessian–gradient product.
   *
   * Computes H·g in-place into `direction` without ever forming H.
   * History is stored in a ring buffer: `histHead` is the next write
   * position, and `histCount` pairs (oldest first) live at logical
   * positions 0..histCount-1, mapped to physical positions via
   * `(histHead − histCount + k + m) % m`.
   *
   * Uses the most recent (s, y) pair for initial Hessian scaling
   * γ = (s · y) / (y · y). Writes into caller-provided `direction`
   * buffer; no allocations.
   */
  private twoLoopRecursion(
    grad: Float64Array,
    S: Float64Array[],
    Y: Float64Array[],
    rho: Float64Array,
    alphaScratch: Float64Array,
    histHead: number,
    count: number,
    m: number,
    n: number,
    direction: Float64Array,
  ): void {
    // direction ← grad (acts as q)
    direction.set(grad);

    // First loop: logical indices count-1 → 0 (newest → oldest)
    for (let k = count - 1; k >= 0; k--) {
      const phys = (histHead - count + k + m) % m;
      const sVec = S[phys];
      const yVec = Y[phys];
      let sq = 0;
      for (let i = 0; i < n; i++) sq += sVec[i] * direction[i];
      const a = rho[phys] * sq;
      alphaScratch[k] = a;
      for (let i = 0; i < n; i++) direction[i] -= a * yVec[i];
    }

    // Scaling by γ using the most recent pair
    const newest = (histHead - 1 + m) % m;
    const sNewest = S[newest];
    const yNewest = Y[newest];
    let sy = 0;
    let yy = 0;
    for (let i = 0; i < n; i++) {
      sy += sNewest[i] * yNewest[i];
      yy += yNewest[i] * yNewest[i];
    }
    const gamma = (yy > 0) ? sy / yy : 1;
    for (let i = 0; i < n; i++) direction[i] *= gamma;

    // Second loop: logical indices 0 → count-1 (oldest → newest)
    for (let k = 0; k < count; k++) {
      const phys = (histHead - count + k + m) % m;
      const sVec = S[phys];
      const yVec = Y[phys];
      let yr = 0;
      for (let i = 0; i < n; i++) yr += yVec[i] * direction[i];
      const beta = rho[phys] * yr;
      const coef = alphaScratch[k] - beta;
      for (let i = 0; i < n; i++) direction[i] += coef * sVec[i];
    }
  }

  /**
   * Armijo backtracking line search.
   *
   * Finds a step size α along `direction` satisfying the Armijo
   * sufficient-decrease condition `f(x + α·d) ≤ f(x) + c₁·α·(∇f·d)`.
   * Starts at `initStep` and halves on each rejection or NaN/Inf eval.
   *
   * Writes the trial point into caller-provided `xNew`, and the status
   * (step, accepted, fNew) into caller-provided `result`. If no finite
   * eval is ever produced, `xNew` is restored to `x` and `stepSize=0`.
   */
  private lineSearch(
    fn: ObjectiveFunction,
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

  /** Async variant of {@link lineSearch}. */
  private async lineSearchAsync(
    fn: AsyncObjectiveFunction,
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
  ): Promise<void> {
    let alpha = initStep;
    let bestAlpha = 0;
    let accepted = false;
    let fNew = fx;
    let foundFinite = false;

    for (let step = 0; step < maxSteps; step++) {
      for (let i = 0; i < n; i++) xNew[i] = x[i] + alpha * direction[i];
      const f = await fn(xNew);

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

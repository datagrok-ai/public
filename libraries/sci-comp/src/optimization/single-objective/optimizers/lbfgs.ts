/**
 * L-BFGS — Limited-memory Broyden–Fletcher–Goldfarb–Shanno.
 *
 * Primary references:
 *  - J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage"
 *    (1980), Math. Comp., 35(151), 773–782.
 *  - D. C. Liu, J. Nocedal, "On the limited memory BFGS method for large
 *    scale optimization" (1989), Math. Programming, 45, 503–528.
 *  - R. H. Byrd, P. Lu, J. Nocedal, "A Limited Memory Algorithm for Bound
 *    Constrained Optimization" (1995), SIAM J. Sci. Comput., 16(5), 1190–1208.
 *
 * Reference implementation:
 *  - SciPy `scipy.optimize._lbfgsb_py._minimize_lbfgsb`.
 *
 * Standard-practice components (textbook, not from the primary papers):
 *  - Armijo backtracking line search — Nocedal & Wright, "Numerical
 *    Optimization" (2nd ed., 2006), Ch. 3.
 *  - Curvature guard `y·s > ε·max(1, y·y)` before history update.
 *
 * Departures from the references:
 *  - No native bound constraints (L-BFGS-**B**): use the penalty layer
 *    (`applyPenalty`) for constraints.
 *  - No `LbfgsInvHessProduct`, `maxfun`, `factr`, or `iprint` from SciPy.
 */
import {Optimizer} from '../optimizer';
import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
  OptimizationResult,
  CommonSettings,
} from '../types';
import {runLBFGSAsync, LBFGSDeps} from './lbfgs-driver';
import {runLBFGSSync} from './lbfgs-driver-sync';

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
/*  Optimizer                                                          */
/* ------------------------------------------------------------------ */

/**
 * Limited-memory BFGS optimizer for unconstrained smooth optimization.
 *
 * Public methods inherited from `Optimizer<S>` delegate via `runInternal` /
 * `runInternalAsync` shims to the top-level functions in `lbfgs-driver.ts`
 * (async, source of truth) and `lbfgs-driver-sync.ts` (codegen-generated).
 */
export class LBFGS extends Optimizer<LBFGSSettings> {
  constructor() {
    super('L-BFGS');
  }

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

  protected runInternal(
    fn: ObjectiveFunction, x0: Float64Array, s: LBFGSSettings,
  ): OptimizationResult {
    return runLBFGSSync(fn, x0, s, this.makeDeps());
  }

  protected runInternalAsync(
    fn: AsyncObjectiveFunction, x0: Float64Array, s: LBFGSSettings,
  ): Promise<OptimizationResult> {
    return runLBFGSAsync(fn, x0, s, this.makeDeps());
  }

  private makeDeps(): LBFGSDeps {
    return {notify: this.notify.bind(this)};
  }
}

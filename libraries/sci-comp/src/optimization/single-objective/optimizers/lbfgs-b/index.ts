/**
 * L-BFGS-B — Limited-memory BFGS with box constraints.
 *
 * Minimises a smooth nonlinear objective subject to simple bounds
 * `l ≤ x ≤ u` (any side may be ±∞). Combines the compact limited-memory
 * BFGS representation of Byrd, Nocedal & Schnabel (1994), the generalized
 * Cauchy point and subspace minimisation of Byrd, Lu, Nocedal & Zhu
 * (1995), the Morales–Nocedal (2011) subspace refinement, and the
 * Moré–Thuente line search (1994).
 *
 * Public API:
 *  - {@link LBFGSB} — the optimizer class (extends `Optimizer`).
 *  - {@link LBFGSBSettings} — settings interface.
 *  - {@link LBFGSBLineSearchSettings} — line-search tuning knobs.
 *  - {@link LBFGSBBounds} — box-constraint input shape.
 *
 * Implementation is split across sibling modules inside this folder
 * (line search, compact representation, Cauchy point, subspace
 * minimisation, bounds helpers, driver).
 *
 * References:
 *  - R. H. Byrd, P. Lu, J. Nocedal, C. Zhu (1995), *A Limited Memory
 *    Algorithm for Bound Constrained Optimization*,
 *    SIAM J. Sci. Comput. 16(5):1190–1208.
 *  - C. Zhu, R. H. Byrd, P. Lu, J. Nocedal (1997), *Algorithm 778:
 *    L-BFGS-B*, ACM TOMS 23(4):550–560.
 *  - J. L. Morales, J. Nocedal (2011), *Remark on Algorithm 778: L-BFGS-B*,
 *    ACM TOMS 38(1):7.
 *  - R. H. Byrd, J. Nocedal, R. B. Schnabel (1994), *Representations of
 *    Quasi-Newton Matrices*, Math. Programming 63:129–156.
 *  - J. J. Moré, D. J. Thuente (1994), *Line Search Algorithms with
 *    Guaranteed Sufficient Decrease*, ACM TOMS 20(3):286–307.
 */
import {Optimizer} from '../../optimizer';
import type {
  ObjectiveFunction,
  AsyncObjectiveFunction,
  OptimizationResult,
} from '../../types';
import type {LBFGSBSettings, LBFGSBLineSearchSettings} from './types';
import {runSync, runAsync} from './driver';

export type {
  LBFGSBSettings,
  LBFGSBLineSearchSettings,
  LBFGSBBounds,
} from './types';

/* ------------------------------------------------------------------ */
/*  Optimizer                                                          */
/* ------------------------------------------------------------------ */

/**
 * L-BFGS-B optimizer — limited-memory BFGS with native box constraints.
 *
 * Usage:
 * ```ts
 * const opt = new LBFGSB();
 * const result = opt.minimize(fn, x0, {
 *   bounds: { lower: 0, upper: [1, 5, 10] },
 *   gradFn: (x, gOut) => { ... },        // optional; FD fallback otherwise
 *   historySize: 10,
 *   gradTolerance: 1e-5,
 * });
 * ```
 */
export class LBFGSB extends Optimizer<LBFGSBSettings> {
  constructor() {
    super('L-BFGS-B');
  }

  /* ----- defaults + validation ------------------------------------- */

  protected withDefaults(s: LBFGSBSettings): LBFGSBSettings {
    const ls: LBFGSBLineSearchSettings = s.lineSearch ?? {};
    const ftol = ls.ftol ?? 1e-3;
    const gtol = ls.gtol ?? 0.9;
    const xtol = ls.xtol ?? 0.1;
    const maxSteps = ls.maxSteps ?? 20;
    const historySize = s.historySize ?? 10;
    const finiteDiffStep = s.finiteDiffStep ?? 1e-7;
    const hasAnalyticGrad = s.gradFn != null || s.gradFnAsync != null;

    if (!Number.isInteger(historySize) || historySize < 1)
      throw new Error('L-BFGS-B: historySize must be an integer ≥ 1');
    if (!(ftol > 0 && ftol < 1))
      throw new Error('L-BFGS-B: lineSearch.ftol must be in (0, 1)');
    if (!(gtol > 0 && gtol < 1))
      throw new Error('L-BFGS-B: lineSearch.gtol must be in (0, 1)');
    if (ftol >= gtol)
      throw new Error('L-BFGS-B: lineSearch.ftol must be < lineSearch.gtol');
    if (!(xtol > 0))
      throw new Error('L-BFGS-B: lineSearch.xtol must be > 0');
    if (!Number.isInteger(maxSteps) || maxSteps < 1)
      throw new Error('L-BFGS-B: lineSearch.maxSteps must be an integer ≥ 1');
    if (!hasAnalyticGrad && !(finiteDiffStep > 0))
      throw new Error('L-BFGS-B: finiteDiffStep must be > 0 when no analytic gradient is provided');

    return {
      maxIterations: s.maxIterations ?? 1_000,
      tolerance: s.tolerance ?? 1e-8,
      historySize,
      gradTolerance: s.gradTolerance ?? 1e-5,
      maxFunctionEvaluations: s.maxFunctionEvaluations ?? 15_000,
      finiteDiffStep,
      bounds: s.bounds,
      gradFn: s.gradFn,
      gradFnAsync: s.gradFnAsync,
      lineSearch: {ftol, gtol, xtol, maxSteps},
      onIteration: s.onIteration,
      constraints: s.constraints,
      penaltyOptions: s.penaltyOptions,
    };
  }

  /* ----- algorithm entry points (filled in by commits 2–7) --------- */

  protected runInternal(
    fn: ObjectiveFunction,
    x0: Float64Array,
    s: LBFGSBSettings,
  ): OptimizationResult {
    return runSync(fn, x0, s as Parameters<typeof runSync>[2]);
  }

  protected async runInternalAsync(
    fn: AsyncObjectiveFunction,
    x0: Float64Array,
    s: LBFGSBSettings,
  ): Promise<OptimizationResult> {
    return runAsync(fn, x0, s as Parameters<typeof runAsync>[2]);
  }
}

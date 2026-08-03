/**
 * Public and internal types for the L-BFGS-B optimizer.
 *
 * `LBFGSBSettings` is the public settings interface. The bound-type
 * constants (`BOUND_FREE` / `BOUND_LOWER` / `BOUND_BOTH` / `BOUND_UPPER`)
 * mirror Nocedal's `nbd[i]` encoding and are used only inside the
 * `lbfgs-b/` folder.
 */
import type {CommonSettings} from '../../types';

/* ------------------------------------------------------------------ */
/*  Internal: bound type encoding (Nocedal's nbd[i])                   */
/* ------------------------------------------------------------------ */

/** Variable bound type. 0=free, 1=lower only, 2=both sides, 3=upper only. */
export type BoundType = 0 | 1 | 2 | 3;

export const BOUND_FREE: BoundType = 0;
export const BOUND_LOWER: BoundType = 1;
export const BOUND_BOTH: BoundType = 2;
export const BOUND_UPPER: BoundType = 3;

/* ------------------------------------------------------------------ */
/*  Public: line search settings                                       */
/* ------------------------------------------------------------------ */

/**
 * Moré–Thuente line-search tuning knobs. Defaults match L-BFGS-B v3.0
 * (the values used in SciPy's `lnsrlb`).
 */
export interface LBFGSBLineSearchSettings {
  /** Armijo sufficient-decrease parameter c₁ ∈ (0, 1). Default 1e-3. */
  ftol?: number;
  /** Curvature parameter c₂ ∈ (c₁, 1). Default 0.9. */
  gtol?: number;
  /** Interval relative tolerance. Default 0.1. */
  xtol?: number;
  /** Maximum line-search iterations per outer step. Default 20. */
  maxSteps?: number;
}

/* ------------------------------------------------------------------ */
/*  Public: box-bounds input                                           */
/* ------------------------------------------------------------------ */

/**
 * Box constraints. Each side accepts a scalar (uniform across all
 * dimensions) or a per-dimension `ArrayLike<number>`. An omitted side
 * defaults to ±∞; individual entries may be `±Infinity` to mark a side
 * unbounded for that variable.
 */
export interface LBFGSBBounds {
  lower?: number | ArrayLike<number>;
  upper?: number | ArrayLike<number>;
}

/* ------------------------------------------------------------------ */
/*  Public: full settings                                              */
/* ------------------------------------------------------------------ */

/**
 * Settings for the L-BFGS-B optimizer (limited-memory BFGS with box
 * constraints).
 *
 * Inherits `maxIterations`, `tolerance`, `onIteration`, `constraints`,
 * and `penaltyOptions` from {@link CommonSettings}. `bounds` (native box
 * constraints) and `constraints` (ineq/eq via penalty) are orthogonal
 * and may be combined.
 */
export interface LBFGSBSettings extends CommonSettings {
  /**
   * Native box constraints. Handled geometrically via the generalized
   * Cauchy point and subspace minimisation — **not** via the penalty
   * layer. Omitted → all variables are (-∞, +∞).
   */
  bounds?: LBFGSBBounds;

  /**
   * Analytic gradient. Writes `∇f(x)` into `gOut`. If omitted, the
   * gradient is approximated by central finite differences with
   * step `finiteDiffStep`.
   */
  gradFn?: (x: Float64Array, gOut: Float64Array) => void;

  /**
   * Async analytic gradient — preferred by `minimizeAsync`/`maximizeAsync`
   * over `gradFn` when both are present. If neither is provided, async
   * central finite differences are used (2n awaits of the objective).
   */
  gradFnAsync?: (x: Float64Array, gOut: Float64Array) => Promise<void>;

  /**
   * Number of correction pairs stored (m in the literature). Higher = more
   * accurate Hessian approximation, more memory/work per iteration.
   * Typical range 3–20. Default 10.
   */
  historySize?: number;

  /**
   * Projected-gradient ∞-norm tolerance — primary convergence test:
   * `‖P(x − g, l, u) − x‖∞ ≤ gradTolerance`. Default 1e-5 (matches
   * SciPy's `pgtol`).
   */
  gradTolerance?: number;

  /**
   * Cap on total f/g evaluations across all line searches. Distinct
   * from `maxIterations` because one outer step can perform up to
   * `lineSearch.maxSteps` evaluations. Default 15000.
   */
  maxFunctionEvaluations?: number;

  /**
   * Step h for the central finite-difference gradient fallback used
   * when `gradFn` / `gradFnAsync` is absent. Default 1e-7.
   */
  finiteDiffStep?: number;

  /** Moré–Thuente line-search parameters. All optional. */
  lineSearch?: LBFGSBLineSearchSettings;
}

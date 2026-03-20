/**
 * Core types for the optimization library.
 */

/** Objective function signature. */
export type ObjectiveFunction = (x: Float64Array) => number;

/** Snapshot of solver state after each iteration — fed to callbacks. */
export interface IterationState {
  iteration: number;
  bestValue: number;
  bestPoint: Float64Array;
  /** Algorithm may attach extra info (e.g. simplex vertices for NM). */
  extra?: Record<string, unknown>;
}

/** Return value of every solver. */
export interface OptimizationResult {
  /** Best point found. */
  point: Float64Array;
  /** Objective value at that point. */
  value: number;
  /** Number of iterations performed. */
  iterations: number;
  /** True if a convergence criterion was met (not just maxIter). */
  converged: boolean;
  /** Per-iteration best value history. */
  costHistory: Float64Array;
}

/** Callback invoked after every iteration. Return `true` to stop early. */
export type IterationCallback = (state: IterationState) => boolean | void;

/** Settings shared by all solvers. */
export interface CommonSettings {
  /** Maximum number of iterations (default: 1000). */
  maxIterations?: number;
  /** Convergence tolerance (default: 1e-8). */
  tolerance?: number;
  /** Called after each iteration; can request early stop. */
  onIteration?: IterationCallback;
  /** Constraints handled internally — penalty is applied with correct sign for both minimize and maximize. */
  constraints?: Constraint[];
  /** Options for the penalty method used with `constraints`. */
  penaltyOptions?: PenaltyOptions;
}

/* ================================================================== */
/*  Constraints & penalty                                              */
/* ================================================================== */

/**
 * Single constraint of the form g(x) <= 0 (inequality)
 * or h(x) = 0 (equality).
 */
export interface Constraint {
  type: 'ineq' | 'eq';
  fn: (x: Float64Array) => number;
}

export type PenaltyMethod = 'quadratic' | 'barrier';

export interface PenaltyOptions {
  /** Penalty method (default: 'quadratic'). */
  method?: PenaltyMethod;
  /** Penalty coefficient μ (default: 1000). */
  mu?: number;
}

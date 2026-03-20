import type {ObjectiveFunction, Constraint, PenaltyOptions} from './types';

/**
 * Wraps an objective function with a penalty term that encodes constraints.
 *
 * The returned function can be passed to any unconstrained solver.
 * The solver minimises f(x) + P(x) without knowing about constraints.
 *
 * Two methods:
 *
 *   'quadratic' (default, exterior penalty):
 *       P(x) = μ · Σ max(0, gᵢ(x))²        for ineq constraints (g ≤ 0)
 *            + μ · Σ hᵢ(x)²                  for eq   constraints (h = 0)
 *     • Works even if x0 is infeasible.
 *     • Larger μ → more accurate but harder landscape.
 *
 *   'barrier' (interior / log-barrier):
 *       P(x) = −μ · Σ ln(−gᵢ(x))            for ineq constraints only
 *     • Requires x0 strictly inside the feasible region.
 *     • Equality constraints are not supported (throws).
 *     • Smaller μ → tighter approximation of the boundary.
 */
export function applyPenalty(
  fn: ObjectiveFunction,
  constraints: Constraint[],
  options: PenaltyOptions = {},
): ObjectiveFunction {
  const method = options.method ?? 'quadratic';
  const mu = options.mu ?? 1000;

  if (method === 'barrier') {
    if (constraints.some((c) => c.type === 'eq'))
      throw new Error('Barrier method does not support equality constraints');

    return (x: Float64Array): number => {
      let penalty = 0;

      for (const c of constraints) {
        const g = c.fn(x);
        if (g >= 0) return Infinity; // outside feasible region
        penalty -= Math.log(-g);
      }

      return fn(x) + mu * penalty;
    };
  }

  // quadratic (exterior penalty)
  return (x: Float64Array): number => {
    let penalty = 0;

    for (const c of constraints) {
      const val = c.fn(x);

      if (c.type === 'ineq') {
        // g(x) <= 0 → penalise positive values
        if (val > 0) penalty += val * val;
      } else {
        // h(x) = 0 → penalise any deviation
        penalty += val * val;
      }
    }

    return fn(x) + mu * penalty;
  };
}

/* ------------------------------------------------------------------ */
/*  Convenience: box constraints → Constraint[]                        */
/* ------------------------------------------------------------------ */

/**
 * Converts box constraints (lower ≤ x ≤ upper) into an array
 * of inequality Constraints suitable for `applyPenalty`.
 *
 * For each dimension i:
 *   x[i] ≥ lower[i]  →  lower[i] - x[i] ≤ 0
 *   x[i] ≤ upper[i]  →  x[i] - upper[i] ≤ 0
 */
export function boxConstraints(
  lower: Float64Array,
  upper: Float64Array,
): Constraint[] {
  const n = lower.length;
  const constraints: Constraint[] = [];

  for (let i = 0; i < n; i++) {
    const lo = lower[i];
    const hi = upper[i];

    if (lo > -Infinity) {
      const idx = i;
      constraints.push({
        type: 'ineq',
        fn: (x) => lo - x[idx], // lo - x ≤ 0  ⟺  x ≥ lo
      });
    }

    if (hi < Infinity) {
      const idx = i;
      constraints.push({
        type: 'ineq',
        fn: (x) => x[idx] - hi, // x - hi ≤ 0  ⟺  x ≤ hi
      });
    }
  }

  return constraints;
}

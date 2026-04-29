/**
 * Box-bounds helper tests.
 *
 * Covers:
 *   - normalizeBounds: scalar / array / missing sides, length mismatch,
 *     NaN, l > u;
 *   - classifyBounds: each of the four nbd codes;
 *   - project: feasible identity, per-side clamping, in-place;
 *   - projectedGradient: zero at optimum / active bound against gradient,
 *     correct ∞-norm, scalar return equals max(|out|);
 *   - maxFeasibleStep: unbounded → ∞, exact boundary distance, rounding
 *     clamp to 0 when x has drifted past bounds.
 */
import {
  normalizeBounds,
  classifyBounds,
  project,
  projectedGradient,
  maxFeasibleStep,
} from '../optimizers/lbfgs-b/bounds';
import {
  BOUND_FREE, BOUND_LOWER, BOUND_BOTH, BOUND_UPPER,
} from '../optimizers/lbfgs-b/types';

/* ================================================================== */
/*  normalizeBounds                                                    */
/* ================================================================== */

describe('normalizeBounds', () => {
  it('undefined input → all free, ±∞', () => {
    const {lower, upper, nbd, anyFinite} = normalizeBounds(undefined, 3);
    expect(Array.from(lower)).toEqual([-Infinity, -Infinity, -Infinity]);
    expect(Array.from(upper)).toEqual([Infinity, Infinity, Infinity]);
    expect(Array.from(nbd)).toEqual([BOUND_FREE, BOUND_FREE, BOUND_FREE]);
    expect(anyFinite).toBe(false);
  });

  it('scalar lower', () => {
    const {lower, upper, nbd, anyFinite} = normalizeBounds({lower: 0}, 3);
    expect(Array.from(lower)).toEqual([0, 0, 0]);
    expect(Array.from(upper)).toEqual([Infinity, Infinity, Infinity]);
    expect(Array.from(nbd)).toEqual([BOUND_LOWER, BOUND_LOWER, BOUND_LOWER]);
    expect(anyFinite).toBe(true);
  });

  it('both sides scalar', () => {
    const {lower, upper, nbd} = normalizeBounds({lower: -1, upper: 1}, 2);
    expect(Array.from(lower)).toEqual([-1, -1]);
    expect(Array.from(upper)).toEqual([1, 1]);
    expect(Array.from(nbd)).toEqual([BOUND_BOTH, BOUND_BOTH]);
  });

  it('per-dim arrays', () => {
    const {lower, upper, nbd} = normalizeBounds(
      {lower: [0, -Infinity, 0], upper: [5, 10, 1]}, 3);
    expect(Array.from(lower)).toEqual([0, -Infinity, 0]);
    expect(Array.from(upper)).toEqual([5, 10, 1]);
    expect(Array.from(nbd)).toEqual([BOUND_BOTH, BOUND_UPPER, BOUND_BOTH]);
  });

  it('mixed scalar + array', () => {
    const {nbd} = normalizeBounds({lower: 0, upper: [1, 2, Infinity]}, 3);
    expect(Array.from(nbd)).toEqual([BOUND_BOTH, BOUND_BOTH, BOUND_LOWER]);
  });

  it('throws on lower array length mismatch', () => {
    expect(() => normalizeBounds({lower: [0, 1]}, 3))
      .toThrow('bounds.lower.length (2)');
  });

  it('throws on upper array length mismatch', () => {
    expect(() => normalizeBounds({upper: [1]}, 3))
      .toThrow('bounds.upper.length (1)');
  });

  it('throws on NaN in scalar', () => {
    expect(() => normalizeBounds({lower: NaN}, 2)).toThrow('NaN');
  });

  it('throws on NaN in array', () => {
    expect(() => normalizeBounds({upper: [1, NaN]}, 2)).toThrow('NaN');
  });

  it('throws on l > u', () => {
    expect(() => normalizeBounds({lower: 1, upper: 0}, 2))
      .toThrow('bounds[0]: lower (1) > upper (0)');
  });

  it('allows l == u (fixed variable)', () => {
    const {nbd} = normalizeBounds({lower: 0.5, upper: 0.5}, 1);
    expect(Array.from(nbd)).toEqual([BOUND_BOTH]);
  });
});

/* ================================================================== */
/*  classifyBounds                                                     */
/* ================================================================== */

describe('classifyBounds', () => {
  it('classifies each nbd code', () => {
    const lower = new Float64Array([-Infinity, 0, -2, -Infinity]);
    const upper = new Float64Array([Infinity, Infinity, 3, 5]);
    const nbd = classifyBounds(lower, upper);
    expect(Array.from(nbd)).toEqual([
      BOUND_FREE, BOUND_LOWER, BOUND_BOTH, BOUND_UPPER,
    ]);
  });
});

/* ================================================================== */
/*  project                                                            */
/* ================================================================== */

describe('project', () => {
  it('identity for already-feasible points', () => {
    const lower = new Float64Array([-1, 0, -5]);
    const upper = new Float64Array([1, 10, 5]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([0.5, 3, 2]);
    const out = new Float64Array(3);
    project(x, lower, upper, nbd, out);
    expect(Array.from(out)).toEqual([0.5, 3, 2]);
  });

  it('clamps to nearest face', () => {
    const lower = new Float64Array([0, 0, -Infinity]);
    const upper = new Float64Array([Infinity, 5, 3]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([-1, 10, 7]);
    const out = new Float64Array(3);
    project(x, lower, upper, nbd, out);
    expect(Array.from(out)).toEqual([0, 5, 3]);
  });

  it('free variables pass through', () => {
    const lower = new Float64Array([-Infinity, -Infinity]);
    const upper = new Float64Array([Infinity, Infinity]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([-1e10, 1e10]);
    const out = new Float64Array(2);
    project(x, lower, upper, nbd, out);
    expect(Array.from(out)).toEqual([-1e10, 1e10]);
  });

  it('in-place (x === out) works', () => {
    const lower = new Float64Array([0, 0]);
    const upper = new Float64Array([1, 1]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([-0.5, 1.5]);
    project(x, lower, upper, nbd, x);
    expect(Array.from(x)).toEqual([0, 1]);
  });

  it('idempotent: project(project(x)) === project(x)', () => {
    const lower = new Float64Array([0, -1, -Infinity]);
    const upper = new Float64Array([1, 1, 5]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([-0.3, 2, 10]);
    const once = new Float64Array(3);
    const twice = new Float64Array(3);
    project(x, lower, upper, nbd, once);
    project(once, lower, upper, nbd, twice);
    expect(Array.from(twice)).toEqual(Array.from(once));
  });
});

/* ================================================================== */
/*  projectedGradient                                                  */
/* ================================================================== */

describe('projectedGradient', () => {
  it('zero at unconstrained minimum (g = 0)', () => {
    const lower = new Float64Array([0, 0]);
    const upper = new Float64Array([1, 1]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([0.5, 0.5]);
    const g = new Float64Array([0, 0]);
    const out = new Float64Array(2);
    const norm = projectedGradient(x, g, lower, upper, nbd, out);
    expect(norm).toBe(0);
    expect(Array.from(out)).toEqual([0, 0]);
  });

  it('zero at active bound with gradient pointing into infeasible region', () => {
    const lower = new Float64Array([0]);
    const upper = new Float64Array([Infinity]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([0]); // on lower bound
    const g = new Float64Array([1]); // g > 0 → want to decrease x, blocked
    const out = new Float64Array(1);
    const norm = projectedGradient(x, g, lower, upper, nbd, out);
    expect(norm).toBe(0);
    expect(out[0]).toBe(0);
  });

  it('equals -g for free variables', () => {
    const lower = new Float64Array([-Infinity, -Infinity]);
    const upper = new Float64Array([Infinity, Infinity]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([3, -2]);
    const g = new Float64Array([1, -4]);
    const out = new Float64Array(2);
    const norm = projectedGradient(x, g, lower, upper, nbd, out);
    expect(Array.from(out)).toEqual([-1, 4]);
    expect(norm).toBe(4);
  });

  it('∞-norm matches max(|out|)', () => {
    const lower = new Float64Array([0, -10, -Infinity]);
    const upper = new Float64Array([5, Infinity, 10]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([2, 0, 5]);
    const g = new Float64Array([-3, 7, -2]);
    const out = new Float64Array(3);
    const norm = projectedGradient(x, g, lower, upper, nbd, out);
    let maxAbs = 0;
    for (let i = 0; i < 3; i++) maxAbs = Math.max(maxAbs, Math.abs(out[i]));
    expect(norm).toBe(maxAbs);
  });
});

/* ================================================================== */
/*  maxFeasibleStep                                                    */
/* ================================================================== */

describe('maxFeasibleStep', () => {
  it('+∞ when no bounds are in play', () => {
    const lower = new Float64Array([-Infinity, -Infinity]);
    const upper = new Float64Array([Infinity, Infinity]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([0, 0]);
    const d = new Float64Array([1, -1]);
    expect(maxFeasibleStep(x, d, lower, upper, nbd)).toBe(Infinity);
  });

  it('exact distance to upper bound', () => {
    const lower = new Float64Array([-Infinity]);
    const upper = new Float64Array([10]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([2]);
    const d = new Float64Array([4]);
    expect(maxFeasibleStep(x, d, lower, upper, nbd)).toBe(2);
  });

  it('exact distance to lower bound', () => {
    const lower = new Float64Array([-3]);
    const upper = new Float64Array([Infinity]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([1]);
    const d = new Float64Array([-2]);
    expect(maxFeasibleStep(x, d, lower, upper, nbd)).toBe(2);
  });

  it('tightest binding coordinate wins', () => {
    const lower = new Float64Array([0, -Infinity]);
    const upper = new Float64Array([5, 10]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([2, 2]);
    const d = new Float64Array([1, 4]);
    // coord 0: (5 − 2)/1 = 3; coord 1: (10 − 2)/4 = 2 → min = 2.
    expect(maxFeasibleStep(x, d, lower, upper, nbd)).toBe(2);
  });

  it('d[i] === 0 contributes +∞ (ignored)', () => {
    const lower = new Float64Array([0, 0]);
    const upper = new Float64Array([1, 1]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([0.5, 0.5]);
    const d = new Float64Array([0, 0.25]);
    // coord 0 ignored, coord 1: (1 − 0.5) / 0.25 = 2.
    expect(maxFeasibleStep(x, d, lower, upper, nbd)).toBe(2);
  });

  it('direction toward unbounded side is ignored', () => {
    const lower = new Float64Array([0]); // only lower bound
    const upper = new Float64Array([Infinity]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([1]);
    const d = new Float64Array([5]); // positive d, no upper → unbounded
    expect(maxFeasibleStep(x, d, lower, upper, nbd)).toBe(Infinity);
  });

  it('negative cap (x drifted past bound) clamps to 0', () => {
    const lower = new Float64Array([0]);
    const upper = new Float64Array([1]);
    const nbd = classifyBounds(lower, upper);
    const x = new Float64Array([1.0001]); // slightly outside upper
    const d = new Float64Array([1]);
    // cap = (1 − 1.0001) / 1 = −0.0001 → clamp to 0.
    expect(maxFeasibleStep(x, d, lower, upper, nbd)).toBe(0);
  });
});

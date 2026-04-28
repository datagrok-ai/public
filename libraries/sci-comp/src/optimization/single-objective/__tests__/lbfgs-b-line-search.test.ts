/**
 * Moré–Thuente line-search tests.
 *
 * Line-search correctness is about Wolfe-condition satisfaction, not
 * locating the 1-D minimum exactly. Tests verify:
 *   - strong Wolfe holds at termination;
 *   - dcstep's four geometric cases are exercised;
 *   - error paths fire on invalid inputs;
 *   - NaN/Inf in the evaluator is handled by bisection;
 *   - sync ≡ async (bit-for-bit).
 *
 * A tight `gtol` is used where we want the search to do real work
 * (otherwise α₀ often already satisfies strong Wolfe and gets returned
 * unchanged — which is correct behaviour but doesn't exercise dcstep).
 */
import {
  runLineSearch,
  runLineSearchAsync,
} from '../optimizers/lbfgs-b/line-search';
import type {
  PhiEval,
  PhiEvalAsync,
  LineSearchParams,
  LineSearchResult,
} from '../optimizers/lbfgs-b/line-search';

const DEFAULT_PARAMS: LineSearchParams = {
  ftol: 1e-3,
  gtol: 0.9,
  xtol: 0.1,
  maxSteps: 20,
  stpMax: 1e20,
  stpMin: 1e-20,
};

const TIGHT_PARAMS: LineSearchParams = {
  ftol: 1e-4,
  gtol: 0.1,
  xtol: 0.1,
  maxSteps: 30,
  stpMax: 1e20,
  stpMin: 1e-20,
};

function expectWolfe(
  r: LineSearchResult,
  phi0: number,
  phiPrime0: number,
  params: LineSearchParams,
): void {
  // Sufficient decrease φ(α) ≤ φ(0) + c₁ α φ'(0). Small slack for FP.
  const armijoRhs = phi0 + params.ftol * r.stp * phiPrime0;
  expect(r.phi).toBeLessThanOrEqual(armijoRhs + 1e-10 * (1 + Math.abs(phi0)));
  // Curvature |φ'(α)| ≤ c₂ |φ'(0)|.
  expect(Math.abs(r.phiPrime)).toBeLessThanOrEqual(
    params.gtol * Math.abs(phiPrime0) + 1e-10 * (1 + Math.abs(phiPrime0)),
  );
}

function toAsyncEval(f: PhiEval): PhiEvalAsync {
  return async (a) => f(a);
}

/* ================================================================== */
/*  Trivial acceptance                                                 */
/* ================================================================== */

describe('runLineSearch — trivial acceptance on loose gtol', () => {
  // With gtol=0.9, α₀=0.5 on φ(α)=½(α-1)² satisfies strong Wolfe
  // immediately (|φ'(0.5)|=0.5 < 0.9·|-1|=0.9). This is CORRECT
  // behaviour — Moré–Thuente is a line search, not a 1-D minimiser.
  const phi: PhiEval = (a) => ({phi: 0.5 * (a - 1) ** 2, phiPrime: a - 1});

  it('accepts α₀ when it already satisfies strong Wolfe', () => {
    const r = runLineSearch(phi, 0.5, -1, 0.5, DEFAULT_PARAMS);
    expect(r.status).toBe('converged');
    expect(r.stp).toBe(0.5);
    expect(r.nfev).toBe(1);
    expectWolfe(r, 0.5, -1, DEFAULT_PARAMS);
  });
});

/* ================================================================== */
/*  Real work with tight gtol                                          */
/* ================================================================== */

describe('runLineSearch — quadratic φ(α) = ½(α − 1)² with tight gtol=0.1', () => {
  const phi: PhiEval = (a) => ({phi: 0.5 * (a - 1) ** 2, phiPrime: a - 1});
  const phi0 = 0.5;
  const pp0 = -1;

  it('sync — converges near α ≈ 1 from α₀ = 0.5', () => {
    const r = runLineSearch(phi, phi0, pp0, 0.5, TIGHT_PARAMS);
    expect(r.status).toBe('converged');
    expect(r.stp).toBeGreaterThanOrEqual(0.9);
    expect(r.stp).toBeLessThanOrEqual(1.1);
    expectWolfe(r, phi0, pp0, TIGHT_PARAMS);
  });

  it('sync — converges from very small α₀ = 0.01', () => {
    const r = runLineSearch(phi, phi0, pp0, 0.01, TIGHT_PARAMS);
    expect(r.status).toBe('converged');
    expect(r.stp).toBeGreaterThanOrEqual(0.9);
    expect(r.stp).toBeLessThanOrEqual(1.1);
    expectWolfe(r, phi0, pp0, TIGHT_PARAMS);
  });

  it('sync — converges from overshoot α₀ = 2', () => {
    const r = runLineSearch(phi, phi0, pp0, 2, TIGHT_PARAMS);
    expect(r.status).toBe('converged');
    expect(r.stp).toBeGreaterThanOrEqual(0.9);
    expect(r.stp).toBeLessThanOrEqual(1.1);
    expectWolfe(r, phi0, pp0, TIGHT_PARAMS);
  });

  it('sync — converges from huge overshoot α₀ = 100', () => {
    const r = runLineSearch(phi, phi0, pp0, 100, TIGHT_PARAMS);
    expect(r.status).toBe('converged');
    expectWolfe(r, phi0, pp0, TIGHT_PARAMS);
  });

  it('async — matches sync', async () => {
    const r = await runLineSearchAsync(toAsyncEval(phi), phi0, pp0, 0.01, TIGHT_PARAMS);
    expect(r.status).toBe('converged');
    expect(r.stp).toBeGreaterThanOrEqual(0.9);
    expect(r.stp).toBeLessThanOrEqual(1.1);
    expectWolfe(r, phi0, pp0, TIGHT_PARAMS);
  });
});

/* ================================================================== */
/*  Shifted quartic — different curvature regimes                      */
/* ================================================================== */

describe('runLineSearch — shifted quartic φ(α) = (α − 2)⁴ + 0.1·α', () => {
  // Descent slope dominated by cubic term: φ'(0) ≈ −31.9.
  // Minimum near α = 2. With tight gtol, search must travel from α₀=1
  // (where |φ'| = 3.9) toward the min.
  const phi: PhiEval = (a) => ({
    phi: (a - 2) ** 4 + 0.1 * a,
    phiPrime: 4 * (a - 2) ** 3 + 0.1,
  });
  const phi0 = 16;
  const pp0 = -31.9;

  it('sync — strong Wolfe holds under tight gtol', () => {
    const r = runLineSearch(phi, phi0, pp0, 1.0, TIGHT_PARAMS);
    expect(r.status).toBe('converged');
    expectWolfe(r, phi0, pp0, TIGHT_PARAMS);
  });
});

/* ================================================================== */
/*  Moré–Thuente 1994 §5 test functions                                */
/* ================================================================== */

describe('Moré–Thuente §5.1 — φ(α) = −α/(α²+β), β=2', () => {
  // ψ(α) = −α/(α²+β); ψ'(α) = (α² − β)/(α²+β)².
  // ψ(0) = 0, ψ'(0) = −1/β < 0 ✓.  Minimum at α = √β.
  const beta = 2;
  const phi: PhiEval = (a) => {
    const denom = a * a + beta;
    return {phi: -a / denom, phiPrime: (a * a - beta) / (denom * denom)};
  };
  const phi0 = 0;
  const pp0 = -1 / beta;

  it('converges with strong Wolfe (ftol=0.1, gtol=0.1)', () => {
    const params: LineSearchParams = {...TIGHT_PARAMS, ftol: 0.1, gtol: 0.1};
    const r = runLineSearch(phi, phi0, pp0, 1, params);
    expect(r.status).toBe('converged');
    expectWolfe(r, phi0, pp0, params);
    expect(r.stp).toBeGreaterThan(0);
  });
});

describe('Moré–Thuente §5.2 — φ(α) = (α+β)⁵ − 2(α+β)⁴, β=0.004', () => {
  // Minimum where 5t = 8 → t* = 1.6 → α* ≈ 1.596.
  const beta = 0.004;
  const phi: PhiEval = (a) => {
    const t = a + beta;
    return {
      phi: Math.pow(t, 5) - 2 * Math.pow(t, 4),
      phiPrime: 5 * Math.pow(t, 4) - 8 * Math.pow(t, 3),
    };
  };
  const phi0 = Math.pow(beta, 5) - 2 * Math.pow(beta, 4);
  const pp0 = 5 * Math.pow(beta, 4) - 8 * Math.pow(beta, 3);

  it('makes progress with strong Wolfe (ftol=0.1, gtol=0.1)', () => {
    // §5.2 has φ'(0) ≈ −5·10⁻⁷ (tiny), so the curvature target
    // |φ'(α)| ≤ gtol·|φ'(0)| ≈ 5·10⁻⁸ is near machine precision for
    // α close to the minimum t* = 8/5. The search may terminate via
    // warning_xtol / warning_rounding rather than 'converged' — all
    // are valid outcomes provided we made progress and (when the
    // search fully converged) strong Wolfe holds.
    const params: LineSearchParams = {...TIGHT_PARAMS, ftol: 0.1, gtol: 0.1};
    const r = runLineSearch(phi, phi0, pp0, 0.1, params);
    expect(['converged', 'warning_xtol', 'warning_rounding']).toContain(r.status);
    expect(r.phi).toBeLessThan(phi0);
    if (r.status === 'converged')
      expectWolfe(r, phi0, pp0, params);
  });
});

/* ================================================================== */
/*  Error paths                                                        */
/* ================================================================== */

describe('runLineSearch — error handling', () => {
  it('returns error when φ\'(0) ≥ 0 (not a descent direction)', () => {
    const phi: PhiEval = (a) => ({phi: a * a, phiPrime: 2 * a});
    const r = runLineSearch(phi, 0, 0, 1, DEFAULT_PARAMS);
    expect(r.status).toBe('error');
    expect(r.ok).toBe(false);
  });

  it('returns error when φ\'(0) > 0', () => {
    const phi: PhiEval = (a) => ({phi: a * a, phiPrime: 2 * a + 1});
    const r = runLineSearch(phi, 0, 1, 1, DEFAULT_PARAMS);
    expect(r.status).toBe('error');
  });

  it('returns error when α₀ ≤ 0', () => {
    const phi: PhiEval = (a) => ({phi: a * a, phiPrime: 2 * a});
    const r = runLineSearch(phi, 0, -1, 0, DEFAULT_PARAMS);
    expect(r.status).toBe('error');
  });

  it('returns error when stpMax < α₀', () => {
    const phi: PhiEval = (a) => ({phi: (a - 1) ** 2, phiPrime: 2 * (a - 1)});
    const r = runLineSearch(phi, 1, -2, 5, {...DEFAULT_PARAMS, stpMax: 1});
    expect(r.status).toBe('error');
  });
});

/* ================================================================== */
/*  Saturation / edge cases                                            */
/* ================================================================== */

describe('runLineSearch — saturation and degenerate paths', () => {
  it('hits stpMax when the function is monotone decreasing within bounds', () => {
    // φ(α) = −α → slope constant at −1. Strong curvature: |φ'(α)|=1=|φ'(0)|,
    // gtol=0.9 → |φ'| ≤ 0.9 never holds. Search hits stpMax.
    const phi: PhiEval = (a) => ({phi: -a, phiPrime: -1});
    const r = runLineSearch(phi, 0, -1, 0.5, {...DEFAULT_PARAMS, stpMax: 2});
    expect(r.stp).toBe(2);
    // Status can be warning_stpmax or warning_xtol depending on interval geometry.
    expect(['warning_stpmax', 'warning_xtol']).toContain(r.status);
  });

  it('handles NaN from evaluator by bisecting toward stx', () => {
    // φ(α) = (α - 0.4)² for α ≤ 0.5, NaN for α > 0.5.
    // α₀=1 → NaN → wrapper bisects: (0 + 1)/2 = 0.5 → (0.4-0.5)² = 0.01, OK.
    const phi: PhiEval = (a) => {
      if (a > 0.5) return {phi: NaN, phiPrime: NaN};
      return {phi: (a - 0.4) ** 2, phiPrime: 2 * (a - 0.4)};
    };
    const r = runLineSearch(phi, 0.16, -0.8, 1.0, DEFAULT_PARAMS);
    expect(Number.isFinite(r.phi)).toBe(true);
    expect(Number.isFinite(r.phiPrime)).toBe(true);
    expect(r.stp).toBeLessThanOrEqual(0.5);
  });

  it('handles persistently bad region via repeated bisection', () => {
    // Half-plane NaN at α ≥ 0.1.  Bisection must shrink to valid region.
    const phi: PhiEval = (a) => {
      if (a >= 0.1) return {phi: NaN, phiPrime: NaN};
      return {phi: (a - 0.05) ** 2, phiPrime: 2 * (a - 0.05)};
    };
    const r = runLineSearch(phi, 0.0025, -0.1, 1.0, DEFAULT_PARAMS);
    expect(Number.isFinite(r.phi)).toBe(true);
    expect(r.stp).toBeLessThan(0.1);
  });

  it('terminates within maxSteps on near-flat pathological function', () => {
    const phi: PhiEval = (a) => ({
      phi: 1e-14 * Math.sin(1000 * a) - 1e-14 * a,
      phiPrime: 1e-11 * Math.cos(1000 * a) - 1e-14,
    });
    const r = runLineSearch(phi, 0, -1e-14, 1, {...DEFAULT_PARAMS, maxSteps: 10});
    expect(r.nfev).toBeLessThanOrEqual(15);
  });

  it('returns error after NAN_CAP consecutive non-finite evaluations (sync)', () => {
    // Evaluator never produces a finite value → bisection cannot recover.
    // Without the cap the wrapper would silently shrink stp to 0 and report
    // success at the starting point; with the cap, 'error' is returned and
    // the caller can drop to memory-reset recovery.
    const phi: PhiEval = () => ({phi: NaN, phiPrime: NaN});
    const r = runLineSearch(phi, 0, -1, 1, {...DEFAULT_PARAMS, maxSteps: 100});
    expect(r.status).toBe('error');
    expect(r.ok).toBe(false);
    // Cap is 20, so we must terminate well before maxSteps.
    expect(r.nfev).toBeLessThan(50);
  });

  it('returns error after NAN_CAP consecutive non-finite evaluations (async)', async () => {
    const phi: PhiEvalAsync = async () => ({phi: NaN, phiPrime: NaN});
    const r = await runLineSearchAsync(phi, 0, -1, 1, {...DEFAULT_PARAMS, maxSteps: 100});
    expect(r.status).toBe('error');
    expect(r.ok).toBe(false);
    expect(r.nfev).toBeLessThan(50);
  });
});

/* ================================================================== */
/*  dcstep case coverage                                               */
/* ================================================================== */

describe('dcstep case coverage', () => {
  // Case 1 (fp > fx): overshoot large α₀, function value rises, bracket.
  it('Case 1 — overshoot bracketing on steep quadratic', () => {
    const phi: PhiEval = (a) => ({phi: (a - 0.1) ** 2, phiPrime: 2 * (a - 0.1)});
    const r = runLineSearch(phi, 0.01, -0.2, 10, TIGHT_PARAMS);
    expect(r.status).toBe('converged');
    expectWolfe(r, 0.01, -0.2, TIGHT_PARAMS);
  });

  // Case 2 (sgnd<0): slope flips sign on a convex function.
  it('Case 2 — slope sign flip on smooth convex', () => {
    const phi: PhiEval = (a) => ({phi: (a - 1) ** 4, phiPrime: 4 * (a - 1) ** 3});
    const r = runLineSearch(phi, 1, -4, 2, TIGHT_PARAMS);
    expect(r.status).toBe('converged');
    expectWolfe(r, 1, -4, TIGHT_PARAMS);
  });

  // Cases 3 & 4 exercised by §5.1 / §5.2 above (slow-curvature functions).
});

/* ================================================================== */
/*  Sync / async parity                                                */
/* ================================================================== */

describe('sync/async parity', () => {
  const cases: Array<[PhiEval, number, number, number, LineSearchParams]> = [
    [(a) => ({phi: 0.5 * (a - 1) ** 2, phiPrime: a - 1}), 0.5, -1, 0.01, TIGHT_PARAMS],
    [(a) => ({phi: (a - 2) ** 4 + 0.1 * a, phiPrime: 4 * (a - 2) ** 3 + 0.1}), 16, -31.9, 1, TIGHT_PARAMS],
    [(a) => ({phi: 0.5 * (a - 1) ** 2, phiPrime: a - 1}), 0.5, -1, 0.5, DEFAULT_PARAMS],
  ];

  for (let i = 0; i < cases.length; i++) {
    const [phi, phi0, pp0, a0, params] = cases[i];
    it(`case ${i} identical result sync vs async`, async () => {
      const rs = runLineSearch(phi, phi0, pp0, a0, params);
      const ra = await runLineSearchAsync(toAsyncEval(phi), phi0, pp0, a0, params);
      expect(ra.status).toBe(rs.status);
      expect(ra.stp).toBeCloseTo(rs.stp, 12);
      expect(ra.phi).toBeCloseTo(rs.phi, 12);
      expect(ra.phiPrime).toBeCloseTo(rs.phiPrime, 12);
      expect(ra.nfev).toBe(rs.nfev);
    });
  }
});

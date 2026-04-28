/**
 * L-BFGS-B scaffold tests.
 *
 * Commits 2–7 will replace this with integration tests (minimize
 * Rosenbrock/Beale/bounded/... variants in the standard sync+async
 * pattern used by other optimizers in this suite). For now the scope
 * is limited to verifying that the public surface, defaults, and
 * validation are wired correctly.
 */
import {LBFGSB, getOptimizer} from '..';

describe('LBFGSB scaffold', () => {
  describe('public surface', () => {
    it('can be instantiated', () => {
      const opt = new LBFGSB();
      expect(opt.name).toBe('L-BFGS-B');
    });

    it('is registered under "l-bfgs-b"', () => {
      const opt = getOptimizer('l-bfgs-b');
      expect(opt).toBeInstanceOf(LBFGSB);
      expect(opt.name).toBe('L-BFGS-B');
    });

    it('registry lookup is case-insensitive', () => {
      expect(getOptimizer('L-BFGS-B')).toBeInstanceOf(LBFGSB);
      expect(getOptimizer('l-BFGS-b')).toBeInstanceOf(LBFGSB);
    });

    it('is distinct from "l-bfgs"', () => {
      expect(getOptimizer('l-bfgs')).not.toBeInstanceOf(LBFGSB);
    });
  });

  describe('settings validation', () => {
    const opt = new LBFGSB();
    const x0 = new Float64Array([0, 0]);
    const f = (): number => 0;

    it('throws on historySize < 1', () => {
      expect(() => opt.minimize(f, x0, {historySize: 0}))
        .toThrow('historySize must be an integer ≥ 1');
    });

    it('throws on non-integer historySize', () => {
      expect(() => opt.minimize(f, x0, {historySize: 3.5}))
        .toThrow('historySize must be an integer ≥ 1');
    });

    it('throws when lineSearch.ftol is out of (0, 1)', () => {
      expect(() => opt.minimize(f, x0, {lineSearch: {ftol: 0}}))
        .toThrow('ftol must be in (0, 1)');
      expect(() => opt.minimize(f, x0, {lineSearch: {ftol: 1}}))
        .toThrow('ftol must be in (0, 1)');
      expect(() => opt.minimize(f, x0, {lineSearch: {ftol: 1.5}}))
        .toThrow('ftol must be in (0, 1)');
    });

    it('throws when lineSearch.gtol is out of (0, 1)', () => {
      expect(() => opt.minimize(f, x0, {lineSearch: {gtol: 0}}))
        .toThrow('gtol must be in (0, 1)');
      expect(() => opt.minimize(f, x0, {lineSearch: {gtol: 1.0001}}))
        .toThrow('gtol must be in (0, 1)');
    });

    it('throws when lineSearch.ftol ≥ lineSearch.gtol', () => {
      expect(() => opt.minimize(f, x0, {lineSearch: {ftol: 0.9, gtol: 0.5}}))
        .toThrow('ftol must be < lineSearch.gtol');
      expect(() => opt.minimize(f, x0, {lineSearch: {ftol: 0.5, gtol: 0.5}}))
        .toThrow('ftol must be < lineSearch.gtol');
    });

    it('throws when lineSearch.xtol ≤ 0', () => {
      expect(() => opt.minimize(f, x0, {lineSearch: {xtol: 0}}))
        .toThrow('xtol must be > 0');
      expect(() => opt.minimize(f, x0, {lineSearch: {xtol: -0.1}}))
        .toThrow('xtol must be > 0');
    });

    it('throws when lineSearch.maxSteps < 1', () => {
      expect(() => opt.minimize(f, x0, {lineSearch: {maxSteps: 0}}))
        .toThrow('maxSteps must be an integer ≥ 1');
    });

    it('throws on non-positive finiteDiffStep without analytic gradient', () => {
      expect(() => opt.minimize(f, x0, {finiteDiffStep: 0}))
        .toThrow('finiteDiffStep must be > 0 when no analytic gradient');
      expect(() => opt.minimize(f, x0, {finiteDiffStep: -1}))
        .toThrow('finiteDiffStep must be > 0 when no analytic gradient');
    });

    it('accepts non-positive finiteDiffStep when gradFn is present', () => {
      // Validation passes; the objective is trivially 0 with gradient 0,
      // so the solver converges immediately at iter 0.
      const r = opt.minimize(f, x0, {
        finiteDiffStep: 0,
        gradFn: (_x, gOut) => {gOut.fill(0);},
      });
      expect(r.converged).toBe(true);
      expect(r.iterations).toBe(0);
    });
  });
});

import {linearFit} from '../tests/linear-fit';
import {expectClose} from './helpers';

const TOL = 1e-3;

describe('linearFit', () => {
  it('recovers an exact line (slope 3, intercept 2)', () => {
    const out = linearFit([1, 2, 3, 4], [5, 8, 11, 14]);
    expectClose(out.slope, 3, TOL, 'slope');
    expectClose(out.intercept, 2, TOL, 'intercept');
    expectClose(out.slopeSe, 0, TOL, 'slopeSe');
    expectClose(out.rSquared, 1, TOL, 'rSquared');
    expect(out.df).toBe(2);
    expect(out.n).toBe(4);
    // Zero SE → CI collapses to the point estimate.
    expectClose(out.slopeCI[0], 3, TOL, 'ci.lo');
    expectClose(out.slopeCI[1], 3, TOL, 'ci.hi');
  });

  it('matches a hand-checked noisy regression (90% slope CI)', () => {
    // x=[1..5], y=[2.1,3.9,6.1,7.9,10.2]; R lm() reference:
    //   slope 2.02, intercept −0.02, slopeSe 0.047610, R² 0.998336,
    //   df 3, t(0.95,3)=2.353363 → 90% CI [1.90795, 2.13205].
    const out = linearFit([1, 2, 3, 4, 5], [2.1, 3.9, 6.1, 7.9, 10.2], {ciLevel: 0.90});
    expectClose(out.slope, 2.02, TOL, 'slope');
    expectClose(out.intercept, -0.02, TOL, 'intercept');
    expectClose(out.slopeSe, 0.047610, TOL, 'slopeSe');
    expectClose(out.rSquared, 0.998336, TOL, 'rSquared');
    expect(out.df).toBe(3);
    expect(out.n).toBe(5);
    expectClose(out.slopeCI[0], 1.90795, TOL, 'ci.lo');
    expectClose(out.slopeCI[1], 2.13205, TOL, 'ci.hi');
  });

  it('widens the CI at a higher confidence level', () => {
    const x = [1, 2, 3, 4, 5];
    const y = [2.1, 3.9, 6.1, 7.9, 10.2];
    const ci90 = linearFit(x, y, {ciLevel: 0.90}).slopeCI;
    const ci99 = linearFit(x, y, {ciLevel: 0.99}).slopeCI;
    expect(ci99[1] - ci99[0]).toBeGreaterThan(ci90[1] - ci90[0]);
  });

  it('drops NaN pairs before fitting', () => {
    const out = linearFit([1, 2, NaN, 4, 5], [5, 8, 99, 14, 17]);
    // Remaining points (1,5),(2,8),(4,14),(5,17) all lie on slope 3 / intercept 2.
    expect(out.n).toBe(4);
    expectClose(out.slope, 3, TOL, 'slope');
    expectClose(out.intercept, 2, TOL, 'intercept');
  });

  it('returns a degenerate (NaN) slope for zero x-spread', () => {
    const out = linearFit([2, 2, 2], [1, 2, 3]);
    expect(Number.isNaN(out.slope)).toBe(true);
    expect(Number.isNaN(out.slopeCI[0])).toBe(true);
    expect(out.n).toBe(3);
    expect(out.df).toBe(1);
  });

  it('returns slope but no SE/CI at df = 0 (two points)', () => {
    const out = linearFit([1, 2], [3, 5]);
    expectClose(out.slope, 2, TOL, 'slope');
    expectClose(out.intercept, 1, TOL, 'intercept');
    expect(Number.isNaN(out.slopeSe)).toBe(true);
    expect(Number.isNaN(out.slopeCI[0])).toBe(true);
    expect(out.df).toBe(0);
    expect(out.n).toBe(2);
  });

  it('rejects a ciLevel outside (0, 1)', () => {
    expect(() => linearFit([1, 2, 3], [1, 2, 3], {ciLevel: 1.5})).toThrow();
    expect(() => linearFit([1, 2, 3], [1, 2, 3], {ciLevel: 0})).toThrow();
  });
});

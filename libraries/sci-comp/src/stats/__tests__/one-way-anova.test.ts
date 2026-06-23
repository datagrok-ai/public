import {oneWayAnova} from '../tests/one-way-anova';
import {expectClose} from './helpers';

const TOL = 1e-3;

describe('oneWayAnova', () => {
  it('matches a hand-checked three-group F-test', () => {
    // Groups 1/2/3 = [1,2,3]/[4,5,6]/[7,8,9]; grand mean 5.
    //   SSB = 3·(9+0+9) = 54, dfB = 2, MSB = 27.
    //   SSW = 2·3 = 6,       dfW = 6, MSW = 1  →  F = 27, df (2, 6).
    const values = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    const groups = [1, 1, 1, 2, 2, 2, 3, 3, 3];
    const out = oneWayAnova(values, groups)!;
    expect(out).not.toBeNull();
    expectClose(out.fStatistic, 27.0, TOL, 'F');
    expect(out.dfBetween).toBe(2);
    expect(out.dfWithin).toBe(6);
    expect(out.pValue).toBeGreaterThan(0);
    expect(out.pValue).toBeLessThan(0.01);
    expect(out.groups).toBe(3);
    expect(out.n).toBe(9);
  });

  it('returns a near-zero F when group means coincide', () => {
    const values = [2, 4, 2, 4, 2, 4];
    const groups = [1, 1, 2, 2, 3, 3];
    const out = oneWayAnova(values, groups)!;
    expect(out).not.toBeNull();
    expectClose(out.fStatistic, 0, TOL, 'F');
    expect(out.pValue).toBeGreaterThan(0.99);
  });

  it('drops NaN observations before the test', () => {
    const out = oneWayAnova([1, 2, 3, NaN, 7, 8, 9], [1, 1, 1, 2, 3, 3, 3]);
    // Group 2 loses its only observation → still 2 groups (1 and 3) survive.
    expect(out).not.toBeNull();
    expect(out!.groups).toBe(2);
    expect(out!.n).toBe(6);
  });

  it('returns null for a single group', () => {
    expect(oneWayAnova([1, 2, 3], [1, 1, 1])).toBeNull();
  });

  it('returns null when there is no within-group variation', () => {
    expect(oneWayAnova([5, 5, 7, 7], [1, 1, 2, 2])).toBeNull();
  });

  it('throws on a length mismatch', () => {
    expect(() => oneWayAnova([1, 2, 3], [1, 1])).toThrow();
  });
});

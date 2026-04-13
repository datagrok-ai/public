import {linearTrend} from '../linear-trend';

describe('linearTrend', () => {
  it('perfect line (ramp)', () => {
    const x = new Float64Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    const r = linearTrend(x, 10);
    expect(r.slope).toBeCloseTo(1.0, 10);
    expect(r.intercept).toBeCloseTo(0.0, 10);
    expect(r.rvalue).toBeCloseTo(1.0, 10);
    expect(r.stderr).toBeCloseTo(0, 10);
    // Perfect fit: slope is significant, pvalue should be 0 (stderr=0 with non-zero slope)
    expect(r.pvalue).toBe(0.0);
  });

  it('constant', () => {
    const x = new Float64Array([5, 5, 5, 5, 5]);
    const r = linearTrend(x, 5);
    expect(r.slope).toBe(0);
    expect(r.intercept).toBe(5.0);
    expect(r.rvalue).toBe(0);
    expect(r.pvalue).toBe(1.0);
    expect(r.stderr).toBe(0);
  });

  it('negative slope', () => {
    const x = new Float64Array([9, 8, 7, 6, 5, 4, 3, 2, 1, 0]);
    const r = linearTrend(x, 10);
    expect(r.slope).toBeCloseTo(-1.0, 10);
    expect(r.rvalue).toBeCloseTo(-1.0, 10);
  });

  it('n=2', () => {
    const x = new Float64Array([1, 3]);
    const r = linearTrend(x, 2);
    expect(r.slope).toBeCloseTo(2.0, 10);
    expect(r.intercept).toBeCloseTo(1.0, 10);
    expect(r.rvalue).toBeCloseTo(1.0, 10);
    expect(r.stderr).toBe(0);
    // Perfect fit with non-zero slope: pvalue should be 0
    expect(r.pvalue).toBe(0.0);
  });

  it('n=3 — first non-trivial pvalue', () => {
    const x = new Float64Array([1, 2, 4]);
    const r = linearTrend(x, 3);
    expect(r.slope).toBeCloseTo(1.5, 10);
    // intercept = meanX - slope*meanT = 7/3 - 1.5*1 = 5/6 ≈ 0.8333
    expect(r.intercept).toBeCloseTo(5 / 6, 10);
    // With df=1, there should be a finite stderr and pvalue
    expect(r.stderr).toBeGreaterThan(0);
    expect(r.pvalue).toBeLessThan(1.0);
    expect(r.pvalue).toBeGreaterThan(0);
  });

  it('n < 2 returns zeros', () => {
    const x = new Float64Array([42]);
    const r = linearTrend(x, 1);
    expect(r.slope).toBe(0);
    expect(r.intercept).toBe(0);
    expect(r.rvalue).toBe(0);
    expect(r.pvalue).toBe(1.0);
    expect(r.stderr).toBe(0);
  });

  it('step function — matches golden values', () => {
    const x = new Float64Array([0, 0, 0, 0, 0, 5, 5, 5, 5, 5]);
    const r = linearTrend(x, 10);
    expect(Math.abs(r.slope - 0.7575757575757576)).toBeLessThan(1e-7);
    expect(Math.abs(r.intercept - (-0.9090909090909092))).toBeLessThan(1e-7);
    expect(Math.abs(r.rvalue - 0.8703882797784891)).toBeLessThan(1e-7);
    expect(Math.abs(r.pvalue - 0.0010528257934054874)).toBeLessThan(1e-6);
    expect(Math.abs(r.stderr - 0.1515151515151515)).toBeLessThan(1e-7);
  });
});

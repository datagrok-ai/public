import {boschlooExact, incidenceExactBoth} from '../tests/boschloo-exact';
import {expectClose, loadFixture} from './helpers';

interface BoschlooCase {
  name: string;
  category: string;
  inputs?: {
    table: number[][];
    alternative: 'two-sided' | 'less' | 'greater';
  };
  expected?: {
    statistic: number | null;
    p_value: number | null;
  };
  // randomised property block
  tables?: Array<{
    table: number[][];
    expected_statistic: number;
    expected_p_value: number;
  }>;
  data_origin?: string;
}

const fixture = loadFixture<BoschlooCase>('boschloo-exact');
const tolP = fixture.metadata.tolerances.p_value;
const tolStat = fixture.metadata.tolerances.statistic;

describe('boschlooExact — fixture from scipy.stats.boschloo_exact', () => {
  for (const c of fixture.cases) {
    if (c.category === 'randomized_property') {
      it(c.name, () => {
        for (const t of c.tables!) {
          const out = boschlooExact(t.table, {alternative: 'two-sided'});
          expectClose(out.statistic, t.expected_statistic, tolStat,
            `${c.name} stat ${JSON.stringify(t.table)}`);
          expectClose(out.pValue, t.expected_p_value, tolP,
            `${c.name} p ${JSON.stringify(t.table)}`);
        }
      });
      continue;
    }
    it(c.name, () => {
      const out = boschlooExact(c.inputs!.table, {alternative: c.inputs!.alternative});
      const exp = c.expected!;
      // The fixture serialises scipy's NaN as JSON `null`; collapse that
      // back to NaN so `expectClose` can compare with `nanIsEqual: true`.
      const expStat = exp.statistic === null ? NaN : exp.statistic;
      const expP = exp.p_value === null ? NaN : exp.p_value;
      expectClose(out.statistic, expStat, tolStat, `${c.name} stat`, /*nanIsEqual=*/ true);
      expectClose(out.pValue, expP, tolP, `${c.name} p`, /*nanIsEqual=*/ true);
    });
  }
});

describe('boschlooExact — input validation', () => {
  it('throws on a non-2×2 table', () => {
    expect(() => boschlooExact([[1, 2, 3], [4, 5, 6]] as number[][])).toThrow(/2×2/);
    expect(() => boschlooExact([[1, 2]] as number[][])).toThrow(/2×2/);
  });

  it('throws on a negative cell', () => {
    expect(() => boschlooExact([[-1, 2], [3, 4]])).toThrow(/non-negative/);
  });

  it('throws on a non-finite cell', () => {
    expect(() => boschlooExact([[1, NaN], [3, 4]])).toThrow(/non-negative/);
  });

  it('throws on an unknown alternative', () => {
    expect(() => boschlooExact([[1, 2], [3, 4]], {alternative: 'foo' as 'less'}))
      .toThrow(/alternative/);
  });

  it('throws when nGrid < 4', () => {
    expect(() => boschlooExact([[1, 2], [3, 4]], {nGrid: 3})).toThrow(/nGrid/);
  });

  it('returns NaN/NaN when a column total is zero (matches scipy)', () => {
    const r = boschlooExact([[0, 3], [0, 3]]);
    expect(Number.isNaN(r.statistic)).toBe(true);
    expect(Number.isNaN(r.pValue)).toBe(true);
  });
});

describe('boschlooExact — algebraic invariants', () => {
  it('two-sided p ≈ 2 · min(less, greater) clipped to 1', () => {
    const tab = [[10, 5], [3, 12]];
    const less = boschlooExact(tab, {alternative: 'less'});
    const greater = boschlooExact(tab, {alternative: 'greater'});
    const ts = boschlooExact(tab, {alternative: 'two-sided'});
    const minOneSided = Math.min(less.pValue, greater.pValue);
    const expected = Math.min(1, 2 * minOneSided);
    expect(Math.abs(ts.pValue - expected)).toBeLessThan(1e-9);
  });

  it('two-sided statistic equals the smaller-p side\'s Fisher statistic', () => {
    const tab = [[10, 5], [3, 12]];
    const less = boschlooExact(tab, {alternative: 'less'});
    const greater = boschlooExact(tab, {alternative: 'greater'});
    const ts = boschlooExact(tab, {alternative: 'two-sided'});
    const minSide = less.pValue < greater.pValue ? less : greater;
    expect(Math.abs(ts.statistic - minSide.statistic)).toBeLessThan(1e-12);
  });

  it('two-sided p is invariant under row transposition', () => {
    // The two-sided test depends only on the unordered row-vs-column structure.
    const a = boschlooExact([[10, 5], [3, 12]], {alternative: 'two-sided'});
    const b = boschlooExact([[3, 12], [10, 5]], {alternative: 'two-sided'});
    expect(Math.abs(a.pValue - b.pValue)).toBeLessThan(1e-9);
  });

  it('two-sided p is invariant under column transposition', () => {
    const a = boschlooExact([[10, 5], [3, 12]], {alternative: 'two-sided'});
    const b = boschlooExact([[5, 10], [12, 3]], {alternative: 'two-sided'});
    expect(Math.abs(a.pValue - b.pValue)).toBeLessThan(1e-9);
  });

  it('one-sided p is symmetric under group swap (less ↔ greater)', () => {
    // Swapping the two columns should swap "less" with "greater".
    const a = boschlooExact([[10, 5], [3, 12]], {alternative: 'less'});
    const b = boschlooExact([[5, 10], [12, 3]], {alternative: 'greater'});
    expect(Math.abs(a.pValue - b.pValue)).toBeLessThan(1e-9);
    expect(Math.abs(a.statistic - b.statistic)).toBeLessThan(1e-12);
  });

  it('Boschloo p ≤ Fisher p on directional shifts (uniform power gain)', () => {
    // Boschloo is uniformly more powerful than Fisher's two-sided test. On a
    // table with a clear directional shift, the Boschloo p-value should not
    // exceed Fisher's p-value (allowing tiny numerical slack).
    const cases = [
      [[10, 5], [3, 12]],
      [[8, 2], [2, 8]],
      [[20, 5], [10, 15]],
    ];
    for (const tab of cases) {
      const b = boschlooExact(tab).pValue;
      // Inline mini-Fisher two-sided computation via incidenceExactBoth
      const both = incidenceExactBoth(tab);
      expect(b).toBeLessThanOrEqual(both.pValueFisher + 1e-9);
    }
  });
});

describe('incidenceExactBoth', () => {
  it('returns Boschloo + Fisher together with the sample odds ratio', () => {
    const r = incidenceExactBoth([[10, 5], [3, 12]]);
    expect(r.oddsRatio).toBeCloseTo((10 * 12) / (5 * 3), 12);
    // Boschloo and Fisher should both be small for this table
    expect(r.pValue).toBeLessThan(0.05);
    expect(r.pValueFisher).toBeLessThan(0.05);
  });

  it('NaN on a degenerate table is collapsed to p = 1.0 (SEND policy)', () => {
    // Zero column total → boschlooExact returns NaN; wrapper substitutes 1.0
    // to match the SEND `incidence_exact_both` semantics.
    const r = incidenceExactBoth([[0, 3], [0, 3]]);
    expect(r.pValue).toBe(1.0);
  });
});

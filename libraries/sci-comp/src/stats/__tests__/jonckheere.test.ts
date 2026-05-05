import {mulberry32} from '../internal/random';
import {jonckheere, JonckheereMethod} from '../tests/jonckheere';
import type {Alternative} from '../types';
import {expectClose, loadFixture} from './helpers';

interface JTCase {
  name: string;
  category: string;
  inputs: {
    x: number[];
    g: number[];
    alternative: Alternative;
    method: JonckheereMethod;
    continuity: boolean;
    nperm: number | null;
  };
  expected: {
    j_statistic: number;
    z_statistic: number | null;
    p_value: number;
    significant: boolean;
  };
}

/** Re-group `(x, g)` flat lists into `groups: number[][]` ordered by sorted unique label. */
function regroup(x: readonly number[], g: readonly number[]): number[][] {
  const labels = Array.from(new Set(g)).sort((a, b) => a - b);
  const buckets = new Map<number, number[]>();
  for (const l of labels) buckets.set(l, []);
  for (let i = 0; i < x.length; i++) buckets.get(g[i])!.push(x[i]);
  return labels.map((l) => buckets.get(l)!);
}

const PERM_SEED = 42;
const PERM_NPERM_DEFAULT = 2000;

function runFixture(fixtureName: string, label: string): void {
  const fixture = loadFixture<JTCase>(fixtureName);
  const tolP = fixture.metadata.tolerances.p_value;
  const tolJ = fixture.metadata.tolerances.j_statistic ?? 0.5;
  const tolZ = fixture.metadata.tolerances.z_statistic;

  describe(label, () => {
    for (const c of fixture.cases) {
      it(c.name, () => {
        const groups = regroup(c.inputs.x, c.inputs.g);
        const out = jonckheere(groups, {
          alternative: c.inputs.alternative,
          method: c.inputs.method,
          continuity: c.inputs.continuity,
          nperm: c.inputs.nperm ?? PERM_NPERM_DEFAULT,
          rng: c.inputs.method === 'permutation' ? mulberry32(PERM_SEED) : undefined,
        });
        expectClose(out.jStatistic, c.expected.j_statistic, tolJ, `${c.name} J`);
        expectClose(out.pValue, c.expected.p_value, tolP, `${c.name} p`);
        if (tolZ !== undefined && c.expected.z_statistic !== null && out.statistic !== null)
          expectClose(out.statistic, c.expected.z_statistic, tolZ, `${c.name} Z`);
      });
    }
  });
}

runFixture('jonckheere-rp-approximate', 'jonckheere — approximate vs regressionpack');
runFixture('jonckheere-rp-exact', 'jonckheere — exact vs regressionpack');
runFixture('jonckheere-clinfun', 'jonckheere — approximate vs R clinfun');
runFixture('jonckheere-pmcmr', 'jonckheere — permutation vs R PMCMRplus');

describe('jonckheere — input validation', () => {
  const G3 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];

  it('throws when fewer than 3 groups are supplied', () => {
    expect(() => jonckheere([[1, 2], [3, 4]])).toThrow(/at least 3/);
  });

  it('error message points two-group callers to mannWhitneyU', () => {
    expect(() => jonckheere([[1, 2], [3, 4]])).toThrow(/mannWhitneyU/);
  });

  it('throws when method = "permutation" and nperm is missing', () => {
    expect(() => jonckheere(G3, {method: 'permutation'})).toThrow(/nperm/);
  });

  it('throws when method = "permutation" and nperm is non-positive', () => {
    expect(() => jonckheere(G3, {method: 'permutation', nperm: 0})).toThrow(/nperm/);
  });

  it('returns null result when all groups are empty', () => {
    const r = jonckheere([[], [], []]);
    expect(r).toEqual({statistic: null, pValue: null, jStatistic: null});
  });

  it('returns null result when N is below MIN_N = 4', () => {
    // N = 3, three non-empty groups — pooled sample is too small.
    const r = jonckheere([[1], [2], [3]]);
    expect(r).toEqual({statistic: null, pValue: null, jStatistic: null});
  });

  it('returns null result when fewer than two groups are non-empty', () => {
    // N = 5 ≥ MIN_N but only one group has data — variance is degenerate.
    const r = jonckheere([[1, 2, 3, 4, 5], [], []]);
    expect(r).toEqual({statistic: null, pValue: null, jStatistic: null});
  });

  it('accepts N = 4 with two non-empty groups', () => {
    // Boundary case: exactly at MIN_N, exactly two non-empty groups.
    const r = jonckheere([[1, 2], [3, 4], []]);
    expect(r.jStatistic).not.toBeNull();
    expect(r.statistic).not.toBeNull();
    expect(r.pValue).not.toBeNull();
  });

  it('strips NaN before computing', () => {
    const a = jonckheere([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
    const b = jonckheere([[1, NaN, 2, 3], [4, 5, 6], [7, NaN, 8, 9]]);
    expect(b.jStatistic).toBe(a.jStatistic);
  });

  it('falls below MIN_N once NaNs are stripped', () => {
    // 6 raw values but only 3 finite → below MIN_N.
    const r = jonckheere([[1, NaN], [NaN, 2], [3, NaN]]);
    expect(r).toEqual({statistic: null, pValue: null, jStatistic: null});
  });
});

describe('jonckheere — alternative direction', () => {
  // Strongly increasing data — `increasing` should give a tiny p-value,
  // `decreasing` ≈ 1, two-sided ≈ 2 × increasing.
  const G = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]];

  it('increasing trend gives a small one-sided p-value', () => {
    const r = jonckheere(G, {alternative: 'increasing'});
    expect(r.pValue).not.toBeNull();
    expect(r.pValue!).toBeLessThan(0.001);
  });

  it('decreasing on increasing data gives a near-1 p-value', () => {
    const r = jonckheere(G, {alternative: 'decreasing'});
    expect(r.pValue).not.toBeNull();
    expect(r.pValue!).toBeGreaterThan(0.99);
  });

  it('two-sided p ≈ 2 × min(increasing, decreasing)', () => {
    const inc = jonckheere(G, {alternative: 'increasing'}).pValue!;
    const ts = jonckheere(G, {alternative: 'two-sided'}).pValue!;
    expect(Math.abs(ts - 2 * inc)).toBeLessThan(1e-12);
  });
});

describe('jonckheere — reproducibility & determinism', () => {
  const G = [[1, 2, 3, 4], [3, 4, 5, 6], [5, 6, 7, 8]];

  it('seeded permutation runs are reproducible', () => {
    const r1 = jonckheere(G, {method: 'permutation', nperm: 500, rng: mulberry32(123)});
    const r2 = jonckheere(G, {method: 'permutation', nperm: 500, rng: mulberry32(123)});
    expect(r1.statistic).toBe(r2.statistic);
    expect(r1.pValue).toBe(r2.pValue);
  });

  it('exact and approximate roughly agree for small N', () => {
    const ex = jonckheere(G, {method: 'exact'});
    const ap = jonckheere(G, {method: 'approximate'});
    expect(ex.jStatistic).toBe(ap.jStatistic);
    expect(Math.abs(ex.pValue! - ap.pValue!)).toBeLessThan(0.05);
  });
});

describe('jonckheere — tie handling (midrank convention)', () => {
  // 1. Untied data: midrank reduces to strict-< (no half-credit terms),
  //    so J is integer-valued and equal to the number of strictly-ordered
  //    cross-group pairs. With G1=[1,4,7], G2=[2,5,8], G3=[3,6,9] there are
  //    18 such pairs (6 per pair-of-groups, by hand).
  it('matches strict-< count when there are no ties', () => {
    const G = [[1, 4, 7], [2, 5, 8], [3, 6, 9]];
    const r = jonckheere(G);
    expect(r.jStatistic).toBe(18);
    expect(Number.isInteger(r.jStatistic!)).toBe(true);
    expect(r.statistic).not.toBeNull();
  });

  // 2. Cross-group ties: J accumulates half-credit. The expected values
  //    here come from the closed-form Lehmann §6.2 / clinfun formulas —
  //    they should match `clinfun::jonckheere.test(method="asymptotic")`
  //    on the same data within numerical precision (the tie-corrected
  //    variance formula is identical). The 144-case clinfun fixture
  //    provides the broader cross-validation against R.
  //
  //    G1=[1,2,3], G2=[2,3,4], G3=[3,4,5]
  //      Pooled multiplicities: t = [2 (val=2), 3 (val=3), 2 (val=4)]
  //      J = 7 + 8.5 + 7 = 22.5
  //      E[J] = (81 − 27)/4 = 13.5
  //      Var[J] = 1356/72 + 108/18144 + 180/576 ≈ 19.15179
  //      Z (continuity)    = (22.5 − 13.5 − 0.5) / √Var ≈ 1.94229
  //      p (continuity, two-sided) ≈ 0.05210
  it('handles cross-group ties with the midrank convention (continuity=true)', () => {
    const G = [[1, 2, 3], [2, 3, 4], [3, 4, 5]];
    const r = jonckheere(G, {continuity: true});
    expect(r.jStatistic).toBe(22.5);
    expect(Math.abs(r.statistic! - 1.9422909597931441)).toBeLessThan(1e-9);
    expect(Math.abs(r.pValue! - 0.052101886806283026)).toBeLessThan(1e-9);
  });

  // 2b. Default omits continuity correction (matches clinfun / PMCMRplus /
  //     SAS PROC FREQ — and the SEND `jonckheere_test` Python library).
  //     Without the −0.5 shift, |Z| is strictly larger and p is smaller.
  //     We assert the relationship algebraically so the test does not
  //     depend on hand-computed normal-CDF values.
  it('default omits continuity correction', () => {
    const G = [[1, 2, 3], [2, 3, 4], [3, 4, 5]];
    const rDefault = jonckheere(G);
    const rNoCC = jonckheere(G, {continuity: false});
    const rCC = jonckheere(G, {continuity: true});

    // Default behaves exactly like an explicit `continuity: false`
    expect(rDefault.statistic).toBe(rNoCC.statistic);
    expect(rDefault.pValue).toBe(rNoCC.pValue);

    // Without continuity, |Z| is larger than with continuity (J > E[J] here)
    expect(Math.abs(rNoCC.statistic!)).toBeGreaterThan(Math.abs(rCC.statistic!));

    // (Z_no_cc − Z_cc) · σ = 0.5 exactly, where σ = (J − E[J]) / Z_no_cc
    const sigma = (22.5 - 13.5) / rNoCC.statistic!;
    expect(Math.abs((rNoCC.statistic! - rCC.statistic!) * sigma - 0.5)).toBeLessThan(1e-12);
  });

  // 3. All values equal: J = E[J] (no information about ordering), and
  //    the tie-corrected variance vanishes → degenerate normal at 0.
  //    Per design choice: return the standard NULL_RESULT (variance ≤ 0
  //    short-circuit), preserving J for inspection.
  it('returns null statistic/pValue when all values are equal (variance = 0)', () => {
    const G = [[5, 5, 5], [5, 5, 5], [5, 5, 5]];
    const r = jonckheere(G);
    // J = ½ · Σ_{i<j} nᵢnⱼ = ½ · 27 = 13.5 = E[J]
    expect(r.jStatistic).toBe(13.5);
    expect(r.statistic).toBeNull();
    expect(r.pValue).toBeNull();
  });
});

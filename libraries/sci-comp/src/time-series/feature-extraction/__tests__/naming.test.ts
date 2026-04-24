import {extractFeatures} from '../extract';
import {makeDataFrame} from './helpers';

describe('feature naming', () => {
  const df = makeDataFrame(new Float64Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), 'value');
  const result = extractFeatures(df);
  const names = result.columns.map((c) => c.name);

  it('follows tsfresh naming convention', () => {
    const pattern = /^[a-zA-Z0-9_]+__[a-zA-Z0-9_]+(__[a-zA-Z0-9_.]+)?$/;
    for (const name of names)
      expect(name).toMatch(pattern);
  });

  it('includes parameterized feature names', () => {
    expect(names).toContain('value__cid_ce__normalize_true');
    expect(names).toContain('value__cid_ce__normalize_false');
    expect(names).toContain('value__symmetry_looking__r_0.05');
    expect(names).toContain('value__linear_trend__attr_slope');
    expect(names).toContain('value__ratio_beyond_r_sigma__r_1.5');
    expect(names).toContain('value__number_crossing_m__m_0');
  });

  it('multi-column names are prefixed correctly', () => {
    const df2 = {
      ids: new Uint32Array([1, 1, 1]),
      time: new Float64Array([0, 1, 2]),
      rowCount: 3,
      columns: [
        {name: 'pH', data: new Float64Array([7.0, 7.1, 7.2])},
        {name: 'conc', data: new Float64Array([1.0, 1.5, 2.0])},
      ],
    };
    const result2 = extractFeatures(df2);
    const names2 = result2.columns.map((c) => c.name);

    const phNames = names2.filter((n) => n.startsWith('pH__'));
    const concNames = names2.filter((n) => n.startsWith('conc__'));

    expect(phNames.length).toBeGreaterThan(0);
    expect(concNames.length).toBeGreaterThan(0);
    expect(phNames.length + concNames.length).toBe(names2.length);
  });
});

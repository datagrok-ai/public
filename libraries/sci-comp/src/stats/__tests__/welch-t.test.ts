import {welchTTest} from '../tests/welch-t';
import {arrayWithNaN, decodeNum, expectClose, loadFixture} from './helpers';

interface WelchCase {
  name: string;
  category: string;
  inputs: {x: (number | null)[]; y: (number | null)[]};
  expected:
    | {statistic: number | string | null; p_value: number | string | null};
}

const fixture = loadFixture<WelchCase>('welch-t');
const tolP = fixture.metadata.tolerances.p_value;
const tolS = fixture.metadata.tolerances.statistic;

describe('welchTTest', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      const x = arrayWithNaN(c.inputs.x);
      const y = arrayWithNaN(c.inputs.y);
      const out = welchTTest(x, y);

      const expS = decodeNum(c.expected.statistic);
      const expP = decodeNum(c.expected.p_value);

      // Tolerances are looser than scipy comparison for "published" cases —
      // the published values themselves are rounded to 4 decimals.
      const useTolS = c.category === 'published_reference' ? 1e-3 : tolS;
      const useTolP = c.category === 'published_reference' ? 1e-3 : tolP;

      // NaN parity for zero-variance edge case
      const nanIsEqual = c.category === 'edge_nan';
      expectClose(out.statistic, expS, useTolS, `${c.name} statistic`, nanIsEqual);
      expectClose(out.pValue, expP, useTolP, `${c.name} p_value`, nanIsEqual);
    });
  }
});


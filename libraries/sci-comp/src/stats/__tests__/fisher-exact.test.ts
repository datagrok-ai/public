import {fisherExact2x2} from '../tests/fisher-exact';
import {decodeNum, expectClose, loadFixture} from './helpers';

interface FisherCase {
  name: string;
  category: string;
  inputs?: {table: number[][]};
  expected?: {
    p_two_sided: number;
    p_greater: number;
    p_less: number;
    odds_ratio: number | string;
  };
  // randomized property cases
  tables?: Array<{table: number[][]; expected_p_two_sided: number}>;
  data_origin?: string;
}

const fixture = loadFixture<FisherCase>('fisher-exact');
const tolP = fixture.metadata.tolerances.p_value;

describe('fisherExact2x2', () => {
  for (const c of fixture.cases) {
    if (c.category === 'randomized_property') {
      it(c.name, () => {
        for (const t of c.tables!) {
          const out = fisherExact2x2(t.table);
          expectClose(out.pValue, t.expected_p_two_sided, 1e-10,
            `${c.name} ${JSON.stringify(t.table)}`);
        }
      });
      continue;
    }
    it(c.name, () => {
      const out = fisherExact2x2(c.inputs!.table);
      const exp = c.expected!;
      expectClose(out.pValue, exp.p_two_sided, tolP, `${c.name} p_two`);
      expectClose(out.pGreater, exp.p_greater, tolP, `${c.name} p_gt`);
      expectClose(out.pLess, exp.p_less, tolP, `${c.name} p_lt`);
      const expOR = decodeNum(exp.odds_ratio);
      if (expOR === Infinity) expect(out.oddsRatio).toBe(Infinity);
      else if (expOR === 0) expect(out.oddsRatio).toBe(0);
      else expectClose(out.oddsRatio, expOR, 1e-10, `${c.name} OR`);
    });
  }
});

import {mannWhitneyU} from '../tests/mann-whitney';
import {arrayWithNaN, expectClose, loadFixture} from './helpers';

interface MWUCase {
  name: string;
  category: string;
  inputs: {x: (number | null)[]; y: (number | null)[]; alternative?: string};
  expected: {statistic: number | null; p_value: number | null};
}

const fixture = loadFixture<MWUCase>('mann-whitney');
const tolP = fixture.metadata.tolerances.p_value;
const tolS = fixture.metadata.tolerances.statistic;

describe('mannWhitneyU', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      const x = arrayWithNaN(c.inputs.x);
      const y = arrayWithNaN(c.inputs.y);
      const out = mannWhitneyU(x, y, 'two-sided');
      expectClose(out.statistic, c.expected.statistic, tolS, `${c.name} statistic`);
      expectClose(out.pValue, c.expected.p_value, tolP, `${c.name} p_value`);
    });
  }
});

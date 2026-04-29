import {dunnettPairwise} from '../tests/dunnett';
import {expectClose, loadFixture} from './helpers';

interface DunnettCase {
  name: string;
  category: string;
  inputs: {
    control: number[];
    treated: {dose_level: number; values: number[]}[];
  };
  expected: {
    dose_level: number;
    p_value: number;
    statistic: number;
  }[];
}

const fixture = loadFixture<DunnettCase>('dunnett');
const tolP = fixture.metadata.tolerances.p_value;
const tolS = fixture.metadata.tolerances.statistic;

describe('dunnettPairwise', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      const treated = c.inputs.treated.map((t) => ({
        doseLevel: t.dose_level,
        values: t.values,
      }));
      const out = dunnettPairwise(c.inputs.control, treated);
      expect(out.length).toBe(c.expected.length);
      for (let i = 0; i < out.length; i++) {
        expect(out[i].doseLevel).toBe(c.expected[i].dose_level);
        expectClose(out[i].statistic, c.expected[i].statistic, tolS,
          `${c.name} dose=${c.expected[i].dose_level} statistic`);
        expectClose(out[i].pValue, c.expected[i].p_value, tolP,
          `${c.name} dose=${c.expected[i].dose_level} p`);
      }
    });
  }
});

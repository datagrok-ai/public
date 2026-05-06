import {welchPairwise} from '../tests/welch-pairwise';
import {arrayWithNaN, expectClose, loadFixture} from './helpers';

interface PairwiseCase {
  name: string;
  category: string;
  inputs: {
    control: (number | null)[];
    treated: {dose_level: number; values: (number | null)[]}[];
  };
  expected: {dose_level: number; p_value_welch: number | null}[];
}

const fixture = loadFixture<PairwiseCase>('welch-pairwise');
const tolP = fixture.metadata.tolerances.p_value;

describe('welchPairwise', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      const ctrl = arrayWithNaN(c.inputs.control);
      const treated = c.inputs.treated.map((t) => ({
        doseLevel: t.dose_level,
        values: arrayWithNaN(t.values),
      }));
      const out = welchPairwise(ctrl, treated);
      expect(out.length).toBe(c.expected.length);
      for (let i = 0; i < out.length; i++) {
        expect(out[i].doseLevel).toBe(c.expected[i].dose_level);
        expectClose(
          out[i].pValueWelch,
          c.expected[i].p_value_welch,
          tolP,
          `${c.name}[${i}] p`,
        );
      }
    });
  }
});

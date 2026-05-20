import {bonferroniCorrect} from '../multiple-comparison/bonferroni';
import {expectClose, loadFixture} from './helpers';

interface BonferroniCase {
  name: string;
  category: string;
  inputs: {p_values: (number | null)[]; n_tests?: number};
  expected: (number | null)[];
}

const fixture = loadFixture<BonferroniCase>('bonferroni');

describe('bonferroniCorrect', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      const out = bonferroniCorrect(c.inputs.p_values, c.inputs.n_tests);
      expect(out.length).toBe(c.expected.length);
      for (let i = 0; i < out.length; i++)
        expectClose(out[i], c.expected[i], 1e-10, `${c.name}[${i}]`);
    });
  }
});

import {jonckheere} from '../tests/jonckheere';
import {expectClose, loadFixture} from './helpers';

interface JTCase {
  name: string;
  category: string;
  inputs: {groups: number[][]};
  expected: {statistic: number | null; p_value: number | null};
}

const fixture = loadFixture<JTCase>('jonckheere');
const tolP = fixture.metadata.tolerances.p_value;
const tolS = fixture.metadata.tolerances.statistic;

describe('jonckheere', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      const out = jonckheere(c.inputs.groups);
      expectClose(out.statistic, c.expected.statistic, tolS, `${c.name} Z`);
      expectClose(out.pValue, c.expected.p_value, tolP, `${c.name} p`);
    });
  }
});

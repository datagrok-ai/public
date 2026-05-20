import {hedgesG} from '../tests/hedges-g';
import {expectClose, loadFixture} from './helpers';

interface HedgesCase {
  name: string;
  category: string;
  inputs: {y1: number[]; y2: number[]};
  expected: null | {
    cohens_d_uncorrected: number;
    j_approx: number;
    j_exact: number;
    hedges_g_with_correction: number;
    df: number;
  };
}

const fixture = loadFixture<HedgesCase>('hedges-g');

describe('hedgesG', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      const out = hedgesG(c.inputs.y1, c.inputs.y2);
      if (c.expected === null) {
        expect(out).toBeNull();
        return;
      }
      // We compute d × J_approx, which should match `hedges_g_with_correction`
      // to high precision (our implementation uses the same formula as Python).
      expectClose(out, c.expected.hedges_g_with_correction, 1e-10, c.name);
    });
  }
});

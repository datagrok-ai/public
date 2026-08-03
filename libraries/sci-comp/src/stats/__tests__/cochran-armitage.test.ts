import {cochranArmitage, cochranArmitageBasic, thresholdTest} from '../tests/cochran-armitage';
import {expectClose, loadFixture} from './helpers';

interface BasicCase {
  name: string;
  category: string;
  inputs: {counts: number[]; totals: number[]};
  expected: {statistic: number | null; p_value: number | null};
}

interface ModifiedCase {
  name: string;
  category: string;
  inputs: {
    counts: number[]; totals: number[];
    scores?: number[];
    alternative?: string;
    variance?: string;
    modified?: boolean;
  };
  expected?: {
    z_statistic: number;
    p_value: number;
    z_modified?: number;
    p_value_modified?: number;
  };
  expected_error?: string;
}

interface ThresholdCase {
  name: string;
  category: string;
  inputs: {
    counts: number[]; totals: number[];
    alpha?: number; adjust_alpha?: boolean;
  };
  expected: Array<{
    z: number | null; p: number | null;
    significant: boolean;
    effect_group?: number | null;
    noel_groups?: number[];
  }>;
  published_zs?: number[];
  published_effect_group?: number | null;
  published_noel_groups?: number[];
  expected_alpha_adj?: number;
}

describe('cochranArmitageBasic', () => {
  const fixture = loadFixture<BasicCase>('cochran-armitage');
  const tolP = fixture.metadata.tolerances.p_value;
  const tolS = fixture.metadata.tolerances.statistic;
  for (const c of fixture.cases) {
    it(c.name, () => {
      const out = cochranArmitageBasic(c.inputs.counts, c.inputs.totals);
      expectClose(out.statistic, c.expected.statistic, tolS, `${c.name} z`);
      expectClose(out.pValue, c.expected.p_value, tolP, `${c.name} p`);
    });
  }
});

describe('cochranArmitage (modified)', () => {
  const fixture = loadFixture<ModifiedCase>('cochran-armitage-modified');
  const tolP = fixture.metadata.tolerances.p_value;
  const tolS = fixture.metadata.tolerances.statistic;
  for (const c of fixture.cases) {
    it(c.name, () => {
      const settings = {
        scores: c.inputs.scores,
        alternative: c.inputs.alternative as any,
        variance: c.inputs.variance as any,
        modified: c.inputs.modified,
      };
      if (c.expected_error) {
        expect(() =>
          cochranArmitage(c.inputs.counts, c.inputs.totals, settings),
        ).toThrow();
        return;
      }
      const out = cochranArmitage(c.inputs.counts, c.inputs.totals, settings);
      const exp = c.expected!;
      // Tolerance for published Young 5b tables: ~0.20 (rounding artifacts in
      // the 1985 paper); apply only when explicitly stated.
      const useTolS = c.category === 'published_dose_response' ? 0.20 : tolS;
      const useTolP = c.category === 'published_dose_response' ? 0.20 : tolP;
      expectClose(out.zStatistic, exp.z_statistic, useTolS, `${c.name} z`);
      expectClose(out.pValue, exp.p_value, useTolP, `${c.name} p`);
      if (exp.z_modified !== undefined)
        expectClose(out.zModified ?? null, exp.z_modified, tolS, `${c.name} z_mod`);

      if (exp.p_value_modified !== undefined)
        expectClose(out.pValueModified ?? null, exp.p_value_modified, tolP, `${c.name} p_mod`);
    });
  }
});

describe('thresholdTest', () => {
  const fixture = loadFixture<ThresholdCase>('cochran-armitage-modified');
  const cases = fixture.threshold_cases ?? [];
  const tolP = fixture.metadata.tolerances.p_value;
  const tolS = fixture.metadata.tolerances.statistic;
  for (const c of cases) {
    it(c.name, () => {
      const out = thresholdTest(c.inputs.counts, c.inputs.totals, {
        alpha: c.inputs.alpha,
        adjustAlpha: c.inputs.adjust_alpha,
      });
      expect(out.length).toBe(c.expected.length);
      for (let i = 0; i < out.length; i++) {
        expectClose(out[i].z, c.expected[i].z, tolS, `${c.name}[${i}].z`);
        expectClose(out[i].p, c.expected[i].p, tolP, `${c.name}[${i}].p`);
        expect(out[i].significant).toBe(c.expected[i].significant);
        if (c.expected[i].effect_group !== undefined)
          expect(out[i].effectGroup ?? null).toBe(c.expected[i].effect_group);
        if (c.expected[i].noel_groups !== undefined)
          expect(out[i].noelGroups).toEqual(c.expected[i].noel_groups);
      }
      if (c.expected_alpha_adj !== undefined && out.length > 0)
        expectClose(out[0].alphaAdj, c.expected_alpha_adj, 1e-10, `${c.name} alpha_adj`);
    });
  }
});

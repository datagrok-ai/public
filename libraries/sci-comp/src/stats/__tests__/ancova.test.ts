import {AncovaResult, runAncova} from '../tests/ancova';
import {expectClose, loadFixture} from './helpers';

interface AncovaCase {
  name: string;
  category: string;
  inputs: {
    organ_values: number[];
    body_weights: number[];
    groups: number[];
    control_group?: number;
    use_organ_free_bw?: boolean;
    alpha?: number;
  };
  expected: AncovaJson | null;
  expected_covariate_mean?: number;
  expected_property?: string;
}

interface AdjustedMeanJson {
  group: number; raw_mean: number; adjusted_mean: number; n: number; se: number;
}
interface PairwiseJson {
  group: number; difference: number; se: number;
  t_statistic: number; p_value: number; significant: boolean;
}
interface EffectDecompJson {
  group: number; total_effect: number; direct_effect: number;
  indirect_effect: number; proportion_direct: number;
  direct_g: number; direct_p: number;
}
interface AncovaJson {
  adjusted_means: AdjustedMeanJson[];
  pairwise: PairwiseJson[];
  slope: {estimate: number; se: number; t_statistic: number; p_value: number};
  slope_homogeneity: {
    f_statistic: number | null; p_value: number | null; homogeneous: boolean;
  };
  effect_decomposition: EffectDecompJson[];
  model_r_squared: number;
  mse: number;
  use_organ_free_bw: boolean;
  covariate_mean: number;
}

const fixture = loadFixture<AncovaCase>('ancova');
// Use 1e-3 — the fixture rounds to 4 decimals, so any 4-decimal match is fine.
const TOL = 1e-3;

describe('runAncova', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      const out = runAncova(
        c.inputs.organ_values,
        c.inputs.body_weights,
        c.inputs.groups,
        {
          controlGroup: c.inputs.control_group,
          useOrganFreeBw: c.inputs.use_organ_free_bw,
          alpha: c.inputs.alpha,
        },
      );
      if (c.expected === null) {
        expect(out).toBeNull();
        return;
      }
      expect(out).not.toBeNull();
      compareAncova(out!, c.expected, c.name);
      if (c.expected_property === 'all_pairwise_p_above_0.05')
        for (const pw of out!.pairwise) expect(pw.pValue).toBeGreaterThan(0.05);
    });
  }
});

function compareAncova(out: AncovaResult, exp: AncovaJson, name: string) {
  expect(out.adjustedMeans.length).toBe(exp.adjusted_means.length);
  for (let i = 0; i < out.adjustedMeans.length; i++) {
    const a = out.adjustedMeans[i];
    const b = exp.adjusted_means[i];
    expect(a.group).toBe(b.group);
    expect(a.n).toBe(b.n);
    expectClose(a.rawMean, b.raw_mean, TOL, `${name} adjMeans[${i}].rawMean`);
    expectClose(a.adjustedMean, b.adjusted_mean, TOL, `${name} adjMeans[${i}].adjustedMean`);
    expectClose(a.se, b.se, TOL, `${name} adjMeans[${i}].se`);
  }

  expect(out.pairwise.length).toBe(exp.pairwise.length);
  for (let i = 0; i < out.pairwise.length; i++) {
    const a = out.pairwise[i];
    const b = exp.pairwise[i];
    expect(a.group).toBe(b.group);
    expectClose(a.difference, b.difference, TOL, `${name} pw[${i}].diff`);
    expectClose(a.se, b.se, TOL, `${name} pw[${i}].se`);
    expectClose(a.tStatistic, b.t_statistic, TOL, `${name} pw[${i}].t`);
    expectClose(a.pValue, b.p_value, TOL, `${name} pw[${i}].p`);
    expect(a.significant).toBe(b.significant);
  }

  expectClose(out.slope.estimate, exp.slope.estimate, TOL, `${name} slope.est`);
  expectClose(out.slope.se, exp.slope.se, TOL, `${name} slope.se`);
  expectClose(out.slope.tStatistic, exp.slope.t_statistic, TOL, `${name} slope.t`);
  expectClose(out.slope.pValue, exp.slope.p_value, TOL, `${name} slope.p`);

  expectClose(out.slopeHomogeneity.fStatistic, exp.slope_homogeneity.f_statistic, TOL, `${name} hom.F`);
  expectClose(out.slopeHomogeneity.pValue, exp.slope_homogeneity.p_value, TOL, `${name} hom.p`);
  expect(out.slopeHomogeneity.homogeneous).toBe(exp.slope_homogeneity.homogeneous);

  expectClose(out.modelRSquared, exp.model_r_squared, TOL, `${name} R²`);
  expectClose(out.mse, exp.mse, TOL, `${name} mse`);
  expect(out.useOrganFreeBw).toBe(exp.use_organ_free_bw);
  expectClose(out.covariateMean, exp.covariate_mean, TOL, `${name} covMean`);
}

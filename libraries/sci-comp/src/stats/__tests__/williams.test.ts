import {pavaDecreasing, pavaIncreasing, williamsTest} from '../tests/williams';
import {lookup1971, lookup1972} from '../tests/williams-tables';
import {expectClose, loadFixture} from './helpers';

interface WilliamsCase {
  name: string;
  category: string;
  inputs: any;
  expected: any;
  tolerance?: number;
}

const fixture = loadFixture<WilliamsCase>('williams');

function decodeV(v: number | string): number {
  if (v === 'inf') return Infinity;
  return v as number;
}

describe('williams', () => {
  for (const c of fixture.cases) {
    it(c.name, () => {
      switch (c.category) {
      case 'pava': {
        const dir = c.inputs.direction;
        const out = dir === 'increase' ?
          pavaIncreasing(c.inputs.values, c.inputs.weights) :
          pavaDecreasing(c.inputs.values, c.inputs.weights);
        const tol = c.tolerance ?? 1e-10;
        expect(out.length).toBe(c.expected.length);
        for (let i = 0; i < out.length; i++)
          expectClose(out[i], c.expected[i], tol, `${c.name}[${i}]`);
        return;
      }
      case 'lookup_1971': {
        const v = decodeV(c.inputs.v);
        const out = lookup1971(c.inputs.k, v, c.inputs.alpha);
        expectClose(out, c.expected, c.tolerance ?? 1e-10, c.name);
        return;
      }
      case 'lookup_1972': {
        const v = decodeV(c.inputs.v);
        const out = lookup1972(c.inputs.i, v, c.inputs.alpha, c.inputs.w);
        expectClose(out, c.expected, c.tolerance ?? 1e-10, c.name);
        return;
      }
      case 'integration':
      case 'integration_edge': {
        const inp = c.inputs;
        const out = williamsTest(
          inp.means, inp.sds, inp.ns, inp.dose_labels,
          {direction: inp.direction, alpha: inp.alpha},
        );
        const exp = c.expected;
        if (exp.minimum_effective_dose !== undefined)
          expect(out.minimumEffectiveDose).toBe(exp.minimum_effective_dose);
        if (exp.n_step_down_results !== undefined)
          expect(out.stepDownResults.length).toBe(exp.n_step_down_results);
        if (exp.first_three_significant)
          for (let i = 0; i < 3; i++) expect(out.stepDownResults[i].significant).toBe(true);

        if (exp.fourth_significant === false)
          expect(out.stepDownResults[3].significant).toBe(false);
        if (exp.t_dose6 !== undefined && out.stepDownResults[0])
          expectClose(out.stepDownResults[0].testStatistic, exp.t_dose6, exp.tol ?? 0.01, `${c.name} t6`);
        if (exp.t_dose3 !== undefined && out.stepDownResults[3])
          expectClose(out.stepDownResults[3].testStatistic, exp.t_dose3, exp.tol ?? 0.01, `${c.name} t3`);
        if (exp.constrained_means !== undefined) {
          const tol = exp.tol_means ?? 0.05;
          for (let i = 0; i < exp.constrained_means.length; i++)
            expectClose(out.constrainedMeans[i], exp.constrained_means[i], tol, `${c.name} cm[${i}]`);
        }
        if (exp.stats_by_dose_index !== undefined) {
          const tol = exp.tol_stats ?? 0.02;
          for (const r of out.stepDownResults) {
            const want = exp.stats_by_dose_index[String(r.doseIndex)];
            if (want !== undefined)
              expectClose(r.testStatistic, want, tol, `${c.name} t${r.doseIndex}`);
          }
        }
        if (exp.dose2_significant === false) {
          const dose2 = out.stepDownResults.find((r) => r.doseIndex === 2);
          if (dose2) expect(dose2.significant).toBe(false);
        }
        if (exp.critical_value_step1 !== undefined && out.stepDownResults.length > 0) {
          expectClose(out.stepDownResults[0].criticalValue, exp.critical_value_step1,
            exp.tol_cv ?? 0.001, `${c.name} cv1`);
        }
        if (exp.all_groups_tested !== undefined)
          expect(out.allGroupsTested).toBe(exp.all_groups_tested);
        if (exp.direction !== undefined)
          expect(out.direction).toBe(exp.direction);
        return;
      }
      }
    });
  }
});

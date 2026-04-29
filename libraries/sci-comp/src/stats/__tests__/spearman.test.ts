import {severityTrend, spearman} from '../tests/spearman';
import {arrayWithNaN, decodeNum, expectClose, loadFixture} from './helpers';

interface SpearmanCase {
  name: string;
  category: string;
  inputs: {x?: (number | null)[]; y?: (number | null)[]};
  expected: {rho: number | string | null; p_value: number | string | null};
}

interface SeverityCase {
  name: string;
  category: string;
  inputs: {dose_levels: (number | null)[]; avg_severities: (number | null)[]};
  expected: {rho: number | string | null; p_value: number | string | null};
}

describe('spearman', () => {
  const fixture = loadFixture<SpearmanCase>('spearman');
  const tolP = fixture.metadata.tolerances.p_value;
  const tolS = fixture.metadata.tolerances.statistic;
  for (const c of fixture.cases) {
    it(c.name, () => {
      const x = arrayWithNaN(c.inputs.x!);
      const y = arrayWithNaN(c.inputs.y!);
      const out = spearman(x, y);
      expectClose(out.rho, decodeNum(c.expected.rho), tolS, `${c.name} rho`, true);
      expectClose(out.pValue, decodeNum(c.expected.p_value), tolP, `${c.name} p_value`, true);
    });
  }
});

describe('severityTrend', () => {
  const fixture = loadFixture<SeverityCase>('severity-trend');
  const tolP = fixture.metadata.tolerances.p_value;
  const tolS = fixture.metadata.tolerances.statistic;
  for (const c of fixture.cases) {
    it(c.name, () => {
      const dl = arrayWithNaN(c.inputs.dose_levels);
      const sev = arrayWithNaN(c.inputs.avg_severities);
      const out = severityTrend(dl, sev);
      expectClose(out.rho, decodeNum(c.expected.rho), tolS, `${c.name} rho`, true);
      expectClose(out.pValue, decodeNum(c.expected.p_value), tolP, `${c.name} p_value`, true);
    });
  }
});

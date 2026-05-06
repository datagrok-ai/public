import {computeNca} from '../compute-nca';
import {ROUTE_IV_BOLUS, ROUTE_PO} from '../types';
import type {ProfileInputs, NcaRules} from '../types';

const DEFAULT_RULES: NcaRules = {
  aucMethod: 'linear-up-log-down',
  blq: {
    preFirstMeasurable: 'set-zero',
    embedded: 'set-zero',
    afterLast: 'set-zero',
    consecutiveAfterLast: 'set-zero',
  },
  lambdaZ: {
    mode: 'auto-best-fit',
    minPoints: 3,
    minRSquared: 0.85,
    excludeCmax: true,
    adjRSquaredFactor: 1e-4,
  },
  extrapWarnPct: 20,
  extrapErrorPct: 50,
  compensatedSummation: false,
};

function poInputs(time: number[], conc: number[], dose: number): ProfileInputs {
  return {
    time: Float64Array.from(time),
    conc: Float64Array.from(conc),
    blqMask: new Uint8Array(time.length),
    lloq: 0.01,
    dose,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    route: ROUTE_PO,
    infusionDuration: null,
    bodyWeight: null,
  };
}

function ivInputs(time: number[], conc: number[], dose: number): ProfileInputs {
  return {...poInputs(time, conc, dose), route: ROUTE_IV_BOLUS};
}

describe('computeNca — clean profiles', () => {
  it('clean PO profile produces all parameters with status "ok"', () => {
    // Build from the synthetic 1-comp model: ka=1.5, ke=0.3, V=0.375, F=0.6, dose=2.5
    // C(t) = F·D/V · ka/(ka−ke) · (exp(−ke·t) − exp(−ka·t))
    const ka = 1.5; const ke = 0.3; const V = 0.375; const F = 0.6; const dose = 2.5;
    const times = [0, 0.25, 0.5, 1, 2, 4, 6, 8, 12];
    const conc = times.map((t) =>
      F * dose / V * ka / (ka - ke) * (Math.exp(-ke * t) - Math.exp(-ka * t)));
    const r = computeNca(poInputs(times, conc, dose), DEFAULT_RULES);
    expect(r.status).toBe('ok');
    expect(r.values.cmax).toBeGreaterThan(0);
    expect(r.values.tmax).toBeGreaterThan(0);
    expect(r.values.aucLast).toBeGreaterThan(0);
    expect(r.values.aucInf).toBeGreaterThan(r.values.aucLast);
    // Best-fit on a profile with an absorption phase doesn't recover ke
    // exactly — the late-tail samples still carry a small absorption
    // remnant. Within ~1% is the realistic expectation here.
    expect(Math.abs(r.values.lambdaZ - ke) / ke).toBeLessThan(0.01);
    expect(Math.abs(r.values.halfLife - Math.LN2 / ke) / (Math.LN2 / ke))
      .toBeLessThan(0.01);
    expect(r.values.cl).toBeCloseTo(dose / r.values.aucInf, 6);
    expect(Math.abs(r.values.vz - dose / (ke * r.values.aucInf)) /
      (dose / (ke * r.values.aucInf))).toBeLessThan(0.02);
  });

  it('clean IV bolus without t=0 → c0 inserted, observed Cmax preserved', () => {
    // C(t) = c0·exp(−ke·t), c0 = 5, ke = 0.4
    const c0True = 5;
    const ke = 0.4;
    const times = [0.25, 0.5, 1, 2, 4, 8];
    const conc = times.map((t) => c0True * Math.exp(-ke * t));
    const r = computeNca(ivInputs(times, conc, 25), DEFAULT_RULES);
    expect(r.status).toBe('ok');
    // Observed Cmax = first observation, NOT the inserted c0.
    expect(r.values.cmax).toBeCloseTo(c0True * Math.exp(-ke * 0.25), 6);
    expect(r.values.tmax).toBe(0.25);
    expect(r.values.lambdaZ).toBeCloseTo(ke, 6);
  });
});

describe('computeNca — BLQ handling', () => {
  it('one BLQ at the start does not break the pipeline (preFirstMeasurable rule)', () => {
    const inputs: ProfileInputs = {
      ...poInputs([0, 0.25, 0.5, 1, 2, 4, 8, 12],
        [0, 0.5, 1.0, 0.8, 0.5, 0.3, 0.1, 0.05], 2.5),
      blqMask: new Uint8Array([1, 0, 0, 0, 0, 0, 0, 0]),
    };
    const r = computeNca(inputs, DEFAULT_RULES);
    expect(r.status).toBe('ok');
    // BLQ at index 0 → preFirstMeasurable set-zero; no NaN propagation.
    expect(Number.isFinite(r.values.aucLast)).toBe(true);
    expect(Number.isFinite(r.values.lambdaZ)).toBe(true);
  });

  it('every point BLQ → status "failed"', () => {
    const inputs: ProfileInputs = {
      ...poInputs([0, 1, 2], [0.005, 0.003, 0.002], 2.5),
      blqMask: new Uint8Array([1, 1, 1]),
    };
    const r = computeNca(inputs, DEFAULT_RULES);
    expect(r.status).toBe('failed');
    expect(Number.isNaN(r.values.cmax)).toBe(true);
  });
});

describe('computeNca — partial fits', () => {
  it('too few post-Cmax points → status "partial", AUClast still reported', () => {
    // 3 points: peak at index 1, only 1 point after → can\'t fit lambda_z
    // (minPoints=3 requires 3 post-Cmax measurable points).
    const r = computeNca(poInputs([0, 1, 2], [0, 5, 3], 2.5), DEFAULT_RULES);
    expect(r.status).toBe('partial');
    expect(Number.isFinite(r.values.cmax)).toBe(true);
    expect(Number.isFinite(r.values.aucLast)).toBe(true);
    expect(Number.isNaN(r.values.lambdaZ)).toBe(true);
    expect(Number.isNaN(r.values.aucInf)).toBe(true);
    expect(Number.isNaN(r.values.cl)).toBe(true);
  });
});

describe('computeNca — warnings', () => {
  it('emits AUC_EXTRAP_HIGH when %extrap > rules.extrapWarnPct', () => {
    // Short tail → high extrapolation. Profile decays slowly relative to
    // the observation window.
    const ke = 0.05; // long half-life
    const times = [0, 1, 2, 3, 4, 5];
    const conc = times.map((t) => Math.exp(-ke * t));
    const r = computeNca(poInputs(times, conc, 2.5), DEFAULT_RULES);
    expect(r.status).toBe('ok');
    expect(r.values.pctExtrap).toBeGreaterThan(20);
    const warning = r.provenance.warnings.find((w) => w.code === 'AUC_EXTRAP_HIGH');
    expect(warning).toBeDefined();
  });

  it('escalates to severity "error" when %extrap > extrapErrorPct', () => {
    const ke = 0.01; // very long half-life → very high extrapolation
    const times = [0, 1, 2, 3, 4, 5];
    const conc = times.map((t) => Math.exp(-ke * t));
    const r = computeNca(poInputs(times, conc, 2.5), DEFAULT_RULES);
    const warning = r.provenance.warnings.find((w) => w.code === 'AUC_EXTRAP_HIGH');
    expect(warning?.severity).toBe('error');
  });
});

describe('computeNca — AUC method selection', () => {
  it('linear method differs from log-linear on a decaying profile', () => {
    const ke = 0.3;
    const times = [0, 0.5, 1, 2, 4, 8];
    const conc = times.map((t) => Math.exp(-ke * t));
    const linRules: NcaRules = {...DEFAULT_RULES, aucMethod: 'linear'};
    const logRules: NcaRules = {...DEFAULT_RULES, aucMethod: 'log-linear'};
    const aucLin = computeNca(poInputs(times, conc, 2.5), linRules).values.aucLast;
    const aucLog = computeNca(poInputs(times, conc, 2.5), logRules).values.aucLast;
    expect(aucLin).not.toBe(aucLog);
    expect(aucLin).toBeGreaterThan(aucLog); // trapezoid over-estimates convex decay
  });

  it('compensatedSummation flag picks the compensated functions', () => {
    const ke = 0.3;
    const times = [0, 0.5, 1, 2, 4, 8];
    const conc = times.map((t) => Math.exp(-ke * t));
    const rNaive = computeNca(poInputs(times, conc, 2.5), DEFAULT_RULES);
    const rComp = computeNca(poInputs(times, conc, 2.5),
      {...DEFAULT_RULES, compensatedSummation: true});
    // On well-conditioned data they agree to round-off.
    expect(Math.abs(rComp.values.aucLast - rNaive.values.aucLast))
      .toBeLessThan(1e-12);
    expect(rComp.provenance.compensated).toBe(true);
  });
});

describe('computeNca — provenance', () => {
  it('reports aucMethod and BLQ trace in provenance', () => {
    const r = computeNca(poInputs([0, 1, 2, 4, 8], [0, 1, 0.7, 0.3, 0.1], 2.5),
      DEFAULT_RULES);
    expect(r.provenance.aucMethod).toBe('linear-up-log-down');
    expect(r.provenance.blqApplied.conc.length).toBe(5);
    expect(r.provenance.lambdaZ).not.toBeNull();
  });
});

import {
  halfLifeFromLambdaZ,
  clearance,
  volumeTerminal,
  pctExtrapolated,
  meanResidenceTime,
  volumeSteadyState,
  pctExtrapolatedAumc,
  tlag,
} from '../derived';

describe('halfLifeFromLambdaZ', () => {
  it('ln(2) / λz', () => {
    expect(halfLifeFromLambdaZ(0.1)).toBeCloseTo(Math.LN2 / 0.1, 12);
    expect(halfLifeFromLambdaZ(Math.LN2)).toBe(1); // λz = ln 2 → t½ = 1
  });

  it('matches the worked example from the PKNCA Indometh subject 1', () => {
    // Reference fixture: lambdaZ = 0.1583204824 → half_life = 4.3781270121.
    expect(halfLifeFromLambdaZ(0.1583204824)).toBeCloseTo(4.3781270121, 8);
  });

  it('returns +Infinity when λz = 0', () => {
    expect(halfLifeFromLambdaZ(0)).toBe(Infinity);
  });

  it('returns NaN for negative λz (caller should not call us with that)', () => {
    expect(Number.isNaN(halfLifeFromLambdaZ(-0.1))).toBe(true);
  });
});

describe('clearance', () => {
  it('dose / AUCinf', () => {
    expect(clearance(25, 2.0351803964)).toBeCloseTo(12.2839233534, 8); // Indometh subj 1
  });

  it('returns +Infinity when AUCinf = 0 and dose > 0', () => {
    expect(clearance(100, 0)).toBe(Infinity);
  });

  it('returns NaN for the 0/0 case', () => {
    expect(Number.isNaN(clearance(0, 0))).toBe(true);
  });

  it('returns 0 when dose = 0', () => {
    expect(clearance(0, 5)).toBe(0);
  });
});

describe('volumeTerminal', () => {
  it('dose / (λz · AUCinf)', () => {
    // Indometh subj 1: dose=25, λz=0.1583204824, AUCinf=2.0351803964
    // Vz = 25 / (0.1583204824 · 2.0351803964) = 77.588971226
    expect(volumeTerminal(25, 0.1583204824, 2.0351803964))
      .toBeCloseTo(77.588971226, 6);
  });

  it('returns +Infinity when λz = 0', () => {
    expect(volumeTerminal(100, 0, 5)).toBe(Infinity);
  });

  it('returns +Infinity when AUCinf = 0', () => {
    expect(volumeTerminal(100, 0.1, 0)).toBe(Infinity);
  });
});

describe('pctExtrapolated', () => {
  it('(AUCinf − AUClast) / AUCinf · 100', () => {
    // Indometh subj 1: AUClast=1.71936529, AUCinf=2.0351803964
    // pct = (2.035 - 1.719) / 2.035 · 100 ≈ 15.5178%
    expect(pctExtrapolated(1.71936529, 2.0351803964))
      .toBeCloseTo(15.5177942452, 6);
  });

  it('returns 0 when AUClast == AUCinf (no extrapolation)', () => {
    expect(pctExtrapolated(5, 5)).toBe(0);
  });

  it('returns 100 when AUClast == 0 and AUCinf > 0', () => {
    expect(pctExtrapolated(0, 5)).toBe(100);
  });

  it('returns -Infinity when AUCinf = 0 with positive AUClast (degenerate)', () => {
    // (0 − 1)/0 = −Infinity. Caller should not reach this with AUCinf=0.
    expect(pctExtrapolated(1, 0)).toBe(-Infinity);
  });

  it('returns NaN for the 0/0 case (AUClast = AUCinf = 0)', () => {
    expect(Number.isNaN(pctExtrapolated(0, 0))).toBe(true);
  });
});

describe('meanResidenceTime', () => {
  it('AUMCinf / AUCinf for bolus/EV (T_inf = 0)', () => {
    expect(meanResidenceTime(10, 2)).toBe(5);
  });

  it('subtracts T_inf/2 for an IV infusion', () => {
    // 10/2 − 2/2 = 5 − 1 = 4
    expect(meanResidenceTime(10, 2, 2)).toBe(4);
  });

  it('returns +Infinity when AUCinf = 0 with positive AUMCinf', () => {
    expect(meanResidenceTime(5, 0)).toBe(Infinity);
  });

  it('returns NaN for the 0/0 case', () => {
    expect(Number.isNaN(meanResidenceTime(0, 0))).toBe(true);
  });
});

describe('volumeSteadyState', () => {
  it('Vss = dose·MRT/AUCinf = dose·AUMCinf/AUCinf² for IV bolus', () => {
    // dose=10, AUMCinf=10, AUCinf=2 → MRT=5 → Vss = 10·5/2 = 25
    expect(volumeSteadyState(10, 10, 2)).toBe(25);
    // identity with the closed form dose·AUMC/AUC²
    expect(volumeSteadyState(10, 10, 2)).toBe(10 * 10 / (2 * 2));
  });

  it('applies the infusion correction via meanResidenceTime', () => {
    // T_inf=2 → MRT = 5 − 1 = 4 → Vss = 10·4/2 = 20
    expect(volumeSteadyState(10, 10, 2, 2)).toBe(20);
    expect(volumeSteadyState(10, 10, 2, 2))
      .toBe(10 * meanResidenceTime(10, 2, 2) / 2);
  });

  it('Indometh subj 1 worked example (bolus, dose = 25)', () => {
    // AUCinf=2.3257135428; pick AUMCinf so MRT≈8h → Vss = 25·8/2.3257 ≈ 86.0
    const aucInf = 2.3257135428;
    const aumcInf = 8 * aucInf; // MRT = 8h by construction
    expect(volumeSteadyState(25, aumcInf, aucInf))
      .toBeCloseTo(25 * 8 / aucInf, 8);
  });
});

describe('pctExtrapolatedAumc', () => {
  it('(AUMCinf − AUMClast) / AUMCinf · 100', () => {
    expect(pctExtrapolatedAumc(8, 10)).toBe(20);
  });

  it('returns 0 when AUMClast == AUMCinf', () => {
    expect(pctExtrapolatedAumc(5, 5)).toBe(0);
  });

  it('returns 100 when AUMClast == 0 and AUMCinf > 0', () => {
    expect(pctExtrapolatedAumc(0, 5)).toBe(100);
  });

  it('returns NaN for the 0/0 case', () => {
    expect(Number.isNaN(pctExtrapolatedAumc(0, 0))).toBe(true);
  });
});

describe('tlag', () => {
  it('time of the last zero before the first positive', () => {
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([0, 0, 5, 3]);
    expect(tlag(time, conc)).toBe(1);
  });

  it('returns the first time when the profile is positive from the start', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([5, 3, 1]);
    expect(tlag(time, conc)).toBe(0);
  });

  it('single zero before the rise', () => {
    const time = new Float64Array([0, 0.5, 1]);
    const conc = new Float64Array([0, 2, 1]);
    expect(tlag(time, conc)).toBe(0);
  });

  it('returns NaN when concentration never rises above zero', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([0, 0, 0]);
    expect(Number.isNaN(tlag(time, conc))).toBe(true);
  });
});

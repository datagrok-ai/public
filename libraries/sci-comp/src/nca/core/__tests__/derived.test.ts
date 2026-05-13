import {
  halfLifeFromLambdaZ,
  clearance,
  volumeTerminal,
  pctExtrapolated,
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

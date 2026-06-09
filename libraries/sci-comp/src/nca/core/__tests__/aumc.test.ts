import {
  aumcLinearNaive,
  aumcLogLinearNaive,
  aumcLinearUpLogDownNaive,
  aumcLinearCompensated,
  aumcLogLinearCompensated,
  aumcLinearUpLogDownCompensated,
  aumcExtrapolateToInfinity,
} from '../aumc';

/**
 * Independent oracle: numerically integrate `t·C(t)` over one log-linear
 * interval, where `C(t) = c1·exp(−k·(t−t1))`, `k = ln(c1/c2)/(t2−t1)`.
 * Composite Simpson with a fine grid — derived from the integrand, NOT from
 * the kernel's closed form, so it genuinely validates the formula.
 */
function simpsonMoment(
  t1: number, c1: number, t2: number, c2: number, n = 100000,
): number {
  const k = Math.log(c1 / c2) / (t2 - t1);
  const f = (t: number) => t * c1 * Math.exp(-k * (t - t1));
  const h = (t2 - t1) / n;
  let s = f(t1) + f(t2);
  for (let i = 1; i < n; i++) s += (i % 2 === 0 ? 2 : 4) * f(t1 + i * h);
  return s * h / 3;
}

describe('aumcLinearNaive', () => {
  it('matches hand-computed value on a 3-point profile', () => {
    // time=[0,1,2], conc=[0,2,4]
    //   (0→1): 1·(0·0 + 1·2)/2 = 1
    //   (1→2): 1·(1·2 + 2·4)/2 = 5
    //   total = 6
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([0, 2, 4]);
    expect(aumcLinearNaive(time, conc, 0, 2)).toBe(6);
  });

  it('handles a single interval', () => {
    // time=[0,2], conc=[3,5]: 2·(0·3 + 2·5)/2 = 10
    const time = new Float64Array([0, 2]);
    const conc = new Float64Array([3, 5]);
    expect(aumcLinearNaive(time, conc, 0, 1)).toBe(10);
  });

  it('returns 0 for an empty range (startIdx == endIdx)', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([1, 2, 3]);
    expect(aumcLinearNaive(time, conc, 1, 1)).toBe(0);
  });

  it('respects startIdx — sub-range integration', () => {
    // Whole: 6 (above). From index 1: 1·(1·2 + 2·4)/2 = 5.
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([0, 2, 4]);
    expect(aumcLinearNaive(time, conc, 1, 2)).toBe(5);
  });
});

describe('aumcLogLinearNaive', () => {
  it('matches the numeric moment integral on a 2-point decay', () => {
    // time=[0,1], conc=[2,1]: log-linear moment of t·C.
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([2, 1]);
    const got = aumcLogLinearNaive(time, conc, 0, 1);
    const expected = simpsonMoment(0, 2, 1, 1);
    expect(Math.abs(got - expected) / expected).toBeLessThan(1e-8);
  });

  it('matches the numeric moment integral on a non-zero start time', () => {
    // The 1/k² term carries a non-trivial t-offset here.
    const time = new Float64Array([3, 7]);
    const conc = new Float64Array([8, 1]);
    const got = aumcLogLinearNaive(time, conc, 0, 1);
    const expected = simpsonMoment(3, 8, 7, 1);
    expect(Math.abs(got - expected) / expected).toBeLessThan(1e-8);
  });

  it('falls back to linear when c1 == c2 (avoids div by k=0)', () => {
    // time=[0,1], conc=[2,2]: linear moment = 1·(0·2 + 1·2)/2 = 1.
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([2, 2]);
    const aumc = aumcLogLinearNaive(time, conc, 0, 1);
    expect(Number.isFinite(aumc)).toBe(true);
    expect(aumc).toBe(1);
  });

  it('falls back to linear when c1 == 0 or c2 == 0', () => {
    // time=[0,1,2], conc=[0,1,0]:
    //   (0→1) c1=0 → linear: 1·(0·0 + 1·1)/2 = 0.5
    //   (1→0) c2=0 → linear: 1·(1·1 + 2·0)/2 = 0.5
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([0, 1, 0]);
    expect(aumcLogLinearNaive(time, conc, 0, 2)).toBe(1);
  });

  it('falls back to linear on negative concentration (defensive)', () => {
    // time=[0,1], conc=[1,-0.5]: linear moment = 1·(0·1 + 1·(−0.5))/2 = −0.25
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([1, -0.5]);
    expect(aumcLogLinearNaive(time, conc, 0, 1)).toBe(-0.25);
  });
});

describe('aumcLinearUpLogDownNaive', () => {
  it('uses linear on ascending intervals', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([0, 2, 4]);
    expect(aumcLinearUpLogDownNaive(time, conc, 0, 2))
      .toBe(aumcLinearNaive(time, conc, 0, 2));
  });

  it('uses log-linear on descending intervals', () => {
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([2, 1]);
    expect(aumcLinearUpLogDownNaive(time, conc, 0, 1))
      .toBe(aumcLogLinearNaive(time, conc, 0, 1));
  });

  it('switches between methods inside one call (up then down)', () => {
    // time=[0,1,2,3], conc=[0,2,1,0.5]
    //   (0→2) ascending → linear moment
    //   (2→1), (1→0.5) descending positive → log moment
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([0, 2, 1, 0.5]);
    const expected =
      aumcLinearNaive(new Float64Array([0, 1]), new Float64Array([0, 2]), 0, 1) +
      simpsonMoment(1, 2, 2, 1) + simpsonMoment(2, 1, 3, 0.5);
    expect(Math.abs(aumcLinearUpLogDownNaive(time, conc, 0, 3) - expected))
      .toBeLessThan(1e-6);
  });

  it('uses linear on flat intervals (c2 == c1)', () => {
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([2, 2]);
    expect(aumcLinearUpLogDownNaive(time, conc, 0, 1)).toBe(1);
  });

  it('uses linear when a descending interval contains a zero', () => {
    // (1→0): linear moment = 1·(0·1 + 1·0)/2 = 0
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([1, 0]);
    expect(aumcLinearUpLogDownNaive(time, conc, 0, 1)).toBe(0);
  });
});

describe('AUMC on a clean mono-exponential — exact moment recovery', () => {
  // C(t) = exp(−k·t). The log-linear interpolant on each interval IS the true
  // exponential (k' = k exactly), so the log-moment kernel integrates t·C
  // exactly. Analytical ∫₀^T t·e^{−kt} dt = 1/k² − e^{−kT}·(T/k + 1/k²).
  const N = 60;
  const k = 0.15;
  const dt = 1;
  const time = new Float64Array(N);
  const conc = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    time[i] = i * dt;
    conc[i] = Math.exp(-k * time[i]);
  }
  const T = (N - 1) * dt;
  const aumcLastAnalytical =
    1 / (k * k) - Math.exp(-k * T) * (T / k + 1 / (k * k));

  it('log-linear AUMClast matches the analytical first moment', () => {
    const got = aumcLogLinearNaive(time, conc, 0, N - 1);
    expect(Math.abs(got - aumcLastAnalytical) / aumcLastAnalytical)
      .toBeLessThan(1e-9);
  });

  it('two-term tail completes AUMClast to the exact infinite moment 1/k²', () => {
    // ∫₀^∞ t·e^{−kt} dt = 1/k². The two-term tail is the EXACT complement of
    // AUMClast — a one-term tail (cLast/k² only) would under-report.
    const aumcLast = aumcLogLinearNaive(time, conc, 0, N - 1);
    const cLast = conc[N - 1];
    const tail = aumcExtrapolateToInfinity(T, cLast, k);
    const aumcInf = aumcLast + tail;
    expect(Math.abs(aumcInf - 1 / (k * k)) / (1 / (k * k))).toBeLessThan(1e-9);
  });

  it('one-term tail (the classic bug) under-reports — guard assertion', () => {
    const aumcLast = aumcLogLinearNaive(time, conc, 0, N - 1);
    const cLast = conc[N - 1];
    const oneTerm = aumcLast + cLast / (k * k); // missing (T·cLast)/k
    expect(oneTerm).toBeLessThan(1 / (k * k));
  });
});

describe('compensated AUMC == naive on well-conditioned data', () => {
  const time = new Float64Array([0, 0.5, 1, 2, 4, 6, 8, 12, 24]);
  const concUp = new Float64Array([0, 1.5, 2.5, 3.0, 2.0, 1.0, 0.5, 0.2, 0.05]);
  const concDown = new Float64Array([10, 8, 6, 4, 2, 1, 0.5, 0.2, 0.05]);

  it('linear', () => {
    const naive = aumcLinearNaive(time, concUp, 0, time.length - 1);
    const comp = aumcLinearCompensated(time, concUp, 0, time.length - 1);
    expect(Math.abs(comp - naive)).toBeLessThan(1e-12 * Math.abs(naive));
  });

  it('log-linear', () => {
    const naive = aumcLogLinearNaive(time, concDown, 0, time.length - 1);
    const comp = aumcLogLinearCompensated(time, concDown, 0, time.length - 1);
    expect(Math.abs(comp - naive)).toBeLessThan(1e-12 * Math.abs(naive));
  });

  it('linear-up/log-down', () => {
    const naive = aumcLinearUpLogDownNaive(time, concUp, 0, time.length - 1);
    const comp = aumcLinearUpLogDownCompensated(time, concUp, 0, time.length - 1);
    expect(Math.abs(comp - naive)).toBeLessThan(1e-12 * Math.abs(naive));
  });

  it('compensated variants honour the same fallback rules as naive', () => {
    // (0→1) ascending, (1→1) flat, (1→0) descending-with-zero — all linear.
    const t = new Float64Array([0, 1, 2, 3]);
    const c = new Float64Array([0, 1, 1, 0]);
    expect(aumcLogLinearCompensated(t, c, 0, 3))
      .toBe(aumcLogLinearNaive(t, c, 0, 3));
    expect(aumcLinearUpLogDownCompensated(t, c, 0, 3))
      .toBe(aumcLinearUpLogDownNaive(t, c, 0, 3));
  });
});

describe('aumcExtrapolateToInfinity (two-term tail)', () => {
  it('(tLast·cLast)/λz + cLast/λz²', () => {
    // tLast=2, cLast=1, λz=0.5 → 2·1/0.5 + 1/0.25 = 4 + 4 = 8
    expect(aumcExtrapolateToInfinity(2, 1, 0.5)).toBe(8);
  });

  it('returns 0 when cLast == 0', () => {
    expect(aumcExtrapolateToInfinity(5, 0, 0.5)).toBe(0);
  });

  it('returns +Infinity when λz == 0 and cLast > 0', () => {
    expect(aumcExtrapolateToInfinity(2, 1, 0)).toBe(Infinity);
  });

  it('returns NaN for the 0/0 case (cLast = 0, λz = 0)', () => {
    expect(Number.isNaN(aumcExtrapolateToInfinity(2, 0, 0))).toBe(true);
  });
});

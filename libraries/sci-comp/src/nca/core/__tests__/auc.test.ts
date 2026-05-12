import {
  aucLinearNaive,
  aucLogLinearNaive,
  aucLinearUpLogDownNaive,
  aucLinearCompensated,
  aucLogLinearCompensated,
  aucLinearUpLogDownCompensated,
  aucExtrapolateToInfinity,
  neumaierSum,
} from '../auc';

/** Build a clean exponential decay profile c(t) = exp(-k·t). */
function exponentialProfile(N: number, k: number, dt: number) {
  const time = new Float64Array(N);
  const conc = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    time[i] = i * dt;
    conc[i] = Math.exp(-k * time[i]);
  }
  return {time, conc};
}

describe('aucLinearNaive', () => {
  it('matches hand-computed value on a 3-point profile', () => {
    // time=[0,1,2], conc=[0,2,4]
    //   AUC = 1·(0+2)/2 + 1·(2+4)/2 = 1 + 3 = 4
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([0, 2, 4]);
    expect(aucLinearNaive(time, conc, 0, 2)).toBe(4);
  });

  it('handles a single interval', () => {
    const time = new Float64Array([0, 0.5]);
    const conc = new Float64Array([2, 4]);
    expect(aucLinearNaive(time, conc, 0, 1)).toBe(1.5); // 0.5·(2+4)/2
  });

  it('returns 0 for an empty range (startIdx == endIdx)', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([1, 2, 3]);
    expect(aucLinearNaive(time, conc, 1, 1)).toBe(0);
  });

  it('respects startIdx — sub-range integration', () => {
    // Whole: 1·(0+2)/2 + 1·(2+4)/2 = 4. From index 1: 1·(2+4)/2 = 3.
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([0, 2, 4]);
    expect(aucLinearNaive(time, conc, 1, 2)).toBe(3);
  });
});

describe('aucLogLinearNaive', () => {
  it('matches hand-computed value on a 2-point decay', () => {
    // time=[0,1], conc=[2,1]: AUC = 1·(2-1)/ln(2) = 1/ln(2)
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([2, 1]);
    expect(aucLogLinearNaive(time, conc, 0, 1)).toBeCloseTo(1 / Math.log(2), 12);
  });

  it('falls back to linear when c1 == c2 (avoids div by ln(1)=0)', () => {
    // time=[0,1], conc=[2,2]: linear gives 1·(2+2)/2 = 2. Log-linear is
    // undefined (ln(1)=0), so the function must fall back to 2, not NaN/Inf.
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([2, 2]);
    const auc = aucLogLinearNaive(time, conc, 0, 1);
    expect(Number.isFinite(auc)).toBe(true);
    expect(auc).toBe(2);
  });

  it('falls back to linear when c1 == 0 or c2 == 0', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc1 = new Float64Array([0, 1, 0]);
    // intervals: (0,1)→linear=0.5; (1,0)→linear=0.5; total=1
    expect(aucLogLinearNaive(time, conc1, 0, 2)).toBe(1);
  });

  it('falls back to linear on negative concentration (defensive)', () => {
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([1, -0.5]);
    // log-linear undefined for negatives → linear: 1·(1-0.5)/2 = 0.25
    expect(aucLogLinearNaive(time, conc, 0, 1)).toBe(0.25);
  });
});

describe('aucLinearUpLogDownNaive', () => {
  it('uses linear on ascending intervals', () => {
    // time=[0,1,2], conc=[0,2,4] — strictly ascending
    // AUC matches aucLinearNaive: 4
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([0, 2, 4]);
    expect(aucLinearUpLogDownNaive(time, conc, 0, 2))
      .toBe(aucLinearNaive(time, conc, 0, 2));
  });

  it('uses log-linear on descending intervals', () => {
    // time=[0,1], conc=[2,1] — descending → log-linear
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([2, 1]);
    expect(aucLinearUpLogDownNaive(time, conc, 0, 1))
      .toBeCloseTo(1 / Math.log(2), 12);
  });

  it('switches between methods inside one call (up then down)', () => {
    // time=[0,1,2,3], conc=[0,2,1,0.5]
    //  (0→2)  ascending → linear: 1·(0+2)/2 = 1
    //  (2→1)  descending, both >0 → log: 1·(2-1)/ln(2) = 1/ln(2)
    //  (1→0.5) descending, both >0 → log: 1·(1-0.5)/ln(2) = 0.5/ln(2)
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([0, 2, 1, 0.5]);
    const expected = 1 + 1 / Math.log(2) + 0.5 / Math.log(2);
    expect(aucLinearUpLogDownNaive(time, conc, 0, 3)).toBeCloseTo(expected, 12);
  });

  it('uses linear on flat intervals (c2 == c1)', () => {
    // [(0,2), (1,2)]: NOT descending (c2 == c1) → linear: 2
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([2, 2]);
    expect(aucLinearUpLogDownNaive(time, conc, 0, 1)).toBe(2);
  });

  it('uses linear when descending interval contains a zero', () => {
    // [(0,1), (1,0)] — descending but c2=0 → linear: 1·(1+0)/2 = 0.5
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([1, 0]);
    expect(aucLinearUpLogDownNaive(time, conc, 0, 1)).toBe(0.5);
  });
});

describe('three methods on a clean exponential decay', () => {
  // c(t) = exp(-k·t), 50 points, dt=1, k=0.1.
  // Analytical AUC from t=0 to t=49 = ∫₀⁴⁹ e^{-0.1·t} dt = 10·(1 - e^{-4.9}).
  const N = 50;
  const k = 0.1;
  const dt = 1;
  const {time, conc} = exponentialProfile(N, k, dt);
  const analytical = (1 - Math.exp(-k * (N - 1) * dt)) / k;

  it('log-linear matches the analytical AUC within 0.01%', () => {
    const auc = aucLogLinearNaive(time, conc, 0, N - 1);
    const relErr = Math.abs(auc - analytical) / analytical;
    expect(relErr).toBeLessThan(1e-4);
  });

  it('linear-up/log-down is identical to log-linear here (all descending)', () => {
    const aucLog = aucLogLinearNaive(time, conc, 0, N - 1);
    const aucMix = aucLinearUpLogDownNaive(time, conc, 0, N - 1);
    expect(aucMix).toBeCloseTo(aucLog, 12);
  });

  it('linear over-estimates the analytical AUC on a convex decay', () => {
    // Trapezoid above an exponential decay → over-estimation.
    const aucLin = aucLinearNaive(time, conc, 0, N - 1);
    expect(aucLin).toBeGreaterThan(analytical);
  });
});

describe('neumaierSum', () => {
  it('handles a basic case', () => {
    expect(neumaierSum([1, 2, 3, 4])).toBe(10);
  });

  it('returns 0 on empty input', () => {
    expect(neumaierSum([])).toBe(0);
  });

  it('cancels round-off on the classic pathological case', () => {
    // Naive: 1 + 1e100 = 1e100 (loses the 1) → +1 = 1e100 → -1e100 = 0.
    // Neumaier: tracks the lost 1's via compensation → 2.
    const naive = [1, 1e100, 1, -1e100].reduce((a, b) => a + b, 0);
    expect(naive).toBe(0); // confirms the reference behaviour we beat
    expect(neumaierSum([1, 1e100, 1, -1e100])).toBe(2);
  });

  it('cancels round-off on N small + 1 large', () => {
    // 1000 contributions of 1, plus ±1e16.
    const N = 1000;
    const values: number[] = [];
    values.push(1e16);
    for (let i = 0; i < N; i++) values.push(1);
    values.push(-1e16);
    const naive = values.reduce((a, b) => a + b, 0);
    const compensated = neumaierSum(values);
    expect(compensated).toBe(N);
    expect(Math.abs(naive - N)).toBeGreaterThan(0); // naive drifts
  });
});

describe('compensated AUC == naive on well-conditioned data', () => {
  it('linear', () => {
    const time = new Float64Array([0, 0.5, 1, 2, 4, 6, 8, 12, 24]);
    const conc = new Float64Array([0, 1.5, 2.5, 3.0, 2.0, 1.0, 0.5, 0.2, 0.05]);
    const naive = aucLinearNaive(time, conc, 0, time.length - 1);
    const comp = aucLinearCompensated(time, conc, 0, time.length - 1);
    expect(Math.abs(comp - naive)).toBeLessThan(1e-15 * Math.abs(naive));
  });

  it('log-linear', () => {
    const time = new Float64Array([0, 0.5, 1, 2, 4, 6, 8, 12, 24]);
    const conc = new Float64Array([10, 8, 6, 4, 2, 1, 0.5, 0.2, 0.05]);
    const naive = aucLogLinearNaive(time, conc, 0, time.length - 1);
    const comp = aucLogLinearCompensated(time, conc, 0, time.length - 1);
    expect(Math.abs(comp - naive)).toBeLessThan(1e-15 * Math.abs(naive));
  });

  it('linear-up/log-down', () => {
    const time = new Float64Array([0, 0.5, 1, 2, 4, 6, 8, 12, 24]);
    const conc = new Float64Array([0, 1.5, 2.5, 3.0, 2.0, 1.0, 0.5, 0.2, 0.05]);
    const naive = aucLinearUpLogDownNaive(time, conc, 0, time.length - 1);
    const comp = aucLinearUpLogDownCompensated(time, conc, 0, time.length - 1);
    expect(Math.abs(comp - naive)).toBeLessThan(1e-15 * Math.abs(naive));
  });

  it('all three compensated honour the same fallback rules as naive', () => {
    // Profile that hits every fallback branch:
    //  (0→1) ascending, both >0  → up/down picks linear
    //  (1→1) flat                → log fallback to linear
    //  (1→0) descending with 0   → up/down fallback to linear
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([0, 1, 1, 0]);
    expect(aucLogLinearCompensated(time, conc, 0, 3))
      .toBe(aucLogLinearNaive(time, conc, 0, 3));
    expect(aucLinearUpLogDownCompensated(time, conc, 0, 3))
      .toBe(aucLinearUpLogDownNaive(time, conc, 0, 3));
  });
});

describe('compensated AUC beats naive on a pathological profile', () => {
  it('1000 small terminal contributions plus a Cmax-like spike', () => {
    // Profile: a single huge contribution sandwiched by many tiny terminal
    // contributions. The huge interval (1e8) followed by the long tail
    // (each ~1e-8) is what trips naive Float64 summation.
    const N = 1002;
    const time = new Float64Array(N);
    const conc = new Float64Array(N);
    // Spike interval: time 0 -> 0.001, conc 0 -> 2e11 (gives dt·c2/2 ≈ 1e8)
    time[0] = 0;
    time[1] = 1e-3;
    conc[0] = 0;
    conc[1] = 2e11;
    // Long tail: each tiny step at small dt with both conc ~ 1e-8.
    // Each interval contributes ~1e-8.
    for (let i = 2; i < N; i++) {
      time[i] = time[i - 1] + 1;
      conc[i] = 1e-8;
    }
    const naive = aucLinearNaive(time, conc, 0, N - 1);
    const comp = aucLinearCompensated(time, conc, 0, N - 1);
    // Both compute the same dominant term (~1e8), but the compensated
    // result preserves the small-tail mass with 0 round-off, while naive
    // drops part of it. So compensated >= naive and the gap is non-zero.
    expect(comp).toBeGreaterThan(naive);
    // The dropped mass is at least ~1e-8 (one tail interval). On a 1e8
    // dominant term that is at the limit of Float64 precision (~1e-16),
    // so even one preserved tail contribution proves the point.
    expect(comp - naive).toBeGreaterThanOrEqual(1e-8);
  });
});

describe('aucExtrapolateToInfinity', () => {
  it('cLast / lambdaZ', () => {
    expect(aucExtrapolateToInfinity(2, 0.5)).toBe(4);
    expect(aucExtrapolateToInfinity(1, 0.1)).toBe(10);
  });

  it('returns 0 when cLast == 0', () => {
    expect(aucExtrapolateToInfinity(0, 0.5)).toBe(0);
  });

  it('returns +Infinity when lambdaZ == 0 and cLast > 0', () => {
    expect(aucExtrapolateToInfinity(2, 0)).toBe(Infinity);
  });

  it('returns NaN for the 0/0 case', () => {
    expect(Number.isNaN(aucExtrapolateToInfinity(0, 0))).toBe(true);
  });
});

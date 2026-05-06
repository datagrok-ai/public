import {mulberry32, deriveWorkerSeeds} from '../prng';

describe('mulberry32', () => {
  it('produces a deterministic regression baseline for seed=0', () => {
    // If this test ever fails after a code change, the PRNG output has
    // shifted — that breaks reproducibility of every bootstrap run that
    // pinned a seed. Update the baseline only on intentional algorithm
    // changes, not as a passing-side patch.
    const expected = [
      0.26642920868471265,
      0.0003297457005828619,
      0.2232720274478197,
      0.1462021479383111,
      0.46732782293111086,
    ];
    const rand = mulberry32(0);
    for (const v of expected)
      expect(rand()).toBe(v);
  });

  it('is deterministic — same seed reproduces the same sequence', () => {
    const a = mulberry32(12345);
    const b = mulberry32(12345);
    for (let i = 0; i < 100; i++)
      expect(a()).toBe(b());
  });

  it('different seeds produce different sequences', () => {
    const a = mulberry32(1);
    const b = mulberry32(2);
    // First few values should diverge — vanishingly unlikely they all match.
    let identical = 0;
    for (let i = 0; i < 10; i++)
      if (a() === b()) identical++;
    expect(identical).toBeLessThan(10);
  });

  it('output stays strictly in [0, 1)', () => {
    const rand = mulberry32(42);
    for (let i = 0; i < 1000; i++) {
      const v = rand();
      expect(v).toBeGreaterThanOrEqual(0);
      expect(v).toBeLessThan(1);
    }
  });

  it('mean of 10 000 samples is approximately 0.5 (±0.01)', () => {
    const rand = mulberry32(123);
    let sum = 0;
    const N = 10_000;
    for (let i = 0; i < N; i++)
      sum += rand();
    const mean = sum / N;
    expect(Math.abs(mean - 0.5)).toBeLessThan(0.01);
  });
});

describe('deriveWorkerSeeds', () => {
  it('produces a deterministic regression baseline for masterSeed=42', () => {
    const seeds = deriveWorkerSeeds(42, 4);
    expect(Array.from(seeds)).toEqual([
      2581720956, 1925393290, 3661312704, 2876485805,
    ]);
  });

  it('returns a Uint32Array of the requested length', () => {
    const seeds = deriveWorkerSeeds(7, 16);
    expect(seeds).toBeInstanceOf(Uint32Array);
    expect(seeds.length).toBe(16);
  });

  it('returns an empty array for workerCount = 0', () => {
    const seeds = deriveWorkerSeeds(7, 0);
    expect(seeds).toBeInstanceOf(Uint32Array);
    expect(seeds.length).toBe(0);
  });

  it('is fully deterministic at fixed (masterSeed, workerCount)', () => {
    const a = deriveWorkerSeeds(99, 8);
    const b = deriveWorkerSeeds(99, 8);
    expect(Array.from(a)).toEqual(Array.from(b));
  });

  it('different master seeds produce different worker seeds', () => {
    const a = deriveWorkerSeeds(1, 4);
    const b = deriveWorkerSeeds(2, 4);
    // All four must differ — vanishing probability of accidental collision.
    let same = 0;
    for (let i = 0; i < 4; i++)
      if (a[i] === b[i]) same++;
    expect(same).toBe(0);
  });
});

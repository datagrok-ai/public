import {estimateC0, insertC0, C0_DEFAULT_METHODS} from '../c0';

describe('estimateC0', () => {
  describe('method "c0" — existing observation at dose time', () => {
    it('uses the conc at t = dose when present and positive', () => {
      const time = new Float64Array([0, 0.25, 0.5]);
      const conc = new Float64Array([2.5, 1.85, 1.39]);
      const blqMask = new Uint8Array( [0, 0, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['c0']})).toBe(2.5);
    });

    it('returns null when the t=0 conc is BLQ', () => {
      const time = new Float64Array([0, 0.25]);
      const conc = new Float64Array([0.04, 1.85]);
      const blqMask = new Uint8Array( [1, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['c0']})).toBeNull();
    });

    it('returns null when no t=0 sample exists', () => {
      const time = new Float64Array([0.25, 0.5]);
      const conc = new Float64Array([1.85, 1.39]);
      const blqMask = new Uint8Array( [0, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['c0']})).toBeNull();
    });
  });

  describe('method "logslope" — back-extrapolation of first two', () => {
    it('matches the formula c0 = c1 · (c1/c2)^(t1/(t2−t1))', () => {
      // Indometh subject 4 reference: c1=1.85, t1=0.25, c2=1.39, t2=0.5
      // c0 expected ≈ 2.46223 (matches PKNCA fixture).
      const time = new Float64Array([0.25, 0.5]);
      const conc = new Float64Array([1.85, 1.39]);
      const blqMask = new Uint8Array( [0, 0]);
      const c0 = estimateC0(time, conc, blqMask, {methods: ['logslope']})!;
      expect(c0).not.toBeNull();
      expect(c0).toBeCloseTo(2.46223, 4);
    });

    it('returns null on ascending profile (c2 >= c1)', () => {
      const time = new Float64Array([1, 2]);
      const conc = new Float64Array([1, 2]);
      const blqMask = new Uint8Array( [0, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['logslope']})).toBeNull();
    });

    it('returns null when c2 == 0 (cannot take log)', () => {
      const time = new Float64Array([1, 2]);
      const conc = new Float64Array([1, 0]);
      const blqMask = new Uint8Array( [0, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['logslope']})).toBeNull();
    });

    it('returns null when fewer than two post-dose measurable points exist', () => {
      const time = new Float64Array([0.25]);
      const conc = new Float64Array([1.85]);
      const blqMask = new Uint8Array( [0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['logslope']})).toBeNull();
    });

    it('skips BLQ points when picking the first two post-dose', () => {
      // First measurable is at index 1 (idx 0 is BLQ).
      const time = new Float64Array([0.25, 0.5, 1]);
      const conc = new Float64Array([0.04, 1.85, 1.39]);
      const blqMask = new Uint8Array( [1, 0, 0]);
      const c0 = estimateC0(time, conc, blqMask, {methods: ['logslope']})!;
      expect(c0).not.toBeNull();
      // Logslope on t=[0.5, 1], c=[1.85, 1.39] → ~slope = ln(1.39/1.85)/0.5
      const slope = (Math.log(1.39) - Math.log(1.85)) / (1 - 0.5);
      const expected = Math.exp(Math.log(1.85) - slope * 0.5);
      expect(c0).toBeCloseTo(expected, 12);
    });

    it('honours timeDose != 0', () => {
      const time = new Float64Array([5.25, 5.5]);
      const conc = new Float64Array([1.85, 1.39]);
      const blqMask = new Uint8Array( [0, 0]);
      const c0 = estimateC0(time, conc, blqMask,
        {methods: ['logslope'], timeDose: 5})!;
      // (t1 − timeDose) = 0.25, same logslope → same numeric result as t1=0.25.
      expect(c0).toBeCloseTo(2.46223, 4);
    });
  });

  describe('method "c1" — first post-dose measurable', () => {
    it('returns the first valid post-dose conc', () => {
      const time = new Float64Array([0.25, 0.5]);
      const conc = new Float64Array([1.85, 1.39]);
      const blqMask = new Uint8Array( [0, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['c1']})).toBe(1.85);
    });

    it('skips BLQ before picking first', () => {
      const time = new Float64Array([0.25, 0.5, 1]);
      const conc = new Float64Array([0.04, 1.85, 1.39]);
      const blqMask = new Uint8Array( [1, 0, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['c1']})).toBe(1.85);
    });

    it('returns null when no post-dose measurable points exist', () => {
      const time = new Float64Array([0.25]);
      const conc = new Float64Array([0.04]);
      const blqMask = new Uint8Array( [1]);
      expect(estimateC0(time, conc, blqMask, {methods: ['c1']})).toBeNull();
    });
  });

  describe('method "cmin" — smallest positive measurable', () => {
    it('returns the minimum positive conc', () => {
      const time = new Float64Array([0.25, 0.5, 1]);
      const conc = new Float64Array([1.85, 0.05, 1.39]);
      const blqMask = new Uint8Array( [0, 0, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['cmin']})).toBe(0.05);
    });

    it('ignores BLQ and zero-valued points', () => {
      const time = new Float64Array([0.25, 0.5, 1]);
      const conc = new Float64Array([1.85, 0.04, 0]);
      const blqMask = new Uint8Array( [0, 1, 0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['cmin']})).toBe(1.85);
    });

    it('returns null when no positive measurable conc exists', () => {
      const time = new Float64Array([0.25, 0.5]);
      const conc = new Float64Array([0.04, 0.03]);
      const blqMask = new Uint8Array( [1, 1]);
      expect(estimateC0(time, conc, blqMask, {methods: ['cmin']})).toBeNull();
    });
  });

  describe('method "set0"', () => {
    it('always returns 0', () => {
      const time = new Float64Array([0.25]);
      const conc = new Float64Array([1.85]);
      const blqMask = new Uint8Array( [0]);
      expect(estimateC0(time, conc, blqMask, {methods: ['set0']})).toBe(0);
    });
  });

  describe('default chain', () => {
    it('PKNCA-equivalent chain: c0 → logslope → c1 → cmin → set0', () => {
      expect(C0_DEFAULT_METHODS).toEqual(['c0', 'logslope', 'c1', 'cmin', 'set0']);
    });

    it('falls through "c0" → "logslope" when no t=0 sample', () => {
      const time = new Float64Array([0.25, 0.5]);
      const conc = new Float64Array([1.85, 1.39]);
      const blqMask = new Uint8Array( [0, 0]);
      const c0 = estimateC0(time, conc, blqMask)!;
      expect(c0).toBeCloseTo(2.46223, 4);
    });

    it('uses "c0" first when t=0 sample is present', () => {
      const time = new Float64Array([0, 0.25, 0.5]);
      const conc = new Float64Array([2.5, 1.85, 1.39]);
      const blqMask = new Uint8Array( [0, 0, 0]);
      expect(estimateC0(time, conc, blqMask)).toBe(2.5);
    });

    it('falls all the way through to "set0" on a degenerate profile', () => {
      // Single ascending pair after dose → c0/logslope/c1/cmin all fail or
      // produce something — but with only one BLQ, all real ones fail and
      // we get 0 from set0.
      const time = new Float64Array([0.25]);
      const conc = new Float64Array([0.04]);
      const blqMask = new Uint8Array( [1]);
      expect(estimateC0(time, conc, blqMask)).toBe(0);
    });
  });
});

describe('insertC0', () => {
  it('prepends an inserted (timeDose, c0) and re-finds Cmax', () => {
    const time = new Float64Array([0.25, 0.5, 1, 2, 4]);
    const conc = new Float64Array([1.85, 1.39, 0.89, 0.4, 0.11]);
    const blqMask = new Uint8Array( [0, 0, 0, 0, 0]);
    const r = insertC0(time, conc, blqMask)!;
    expect(r).not.toBeNull();
    expect(r.time.length).toBe(6);
    expect(r.time[0]).toBe(0);
    expect(r.conc[0]).toBeCloseTo(2.46223, 4);
    expect(r.c0).toBeCloseTo(2.46223, 4);
    // Cmax now at the inserted point (highest concentration).
    expect(r.cmaxIdx).toBe(0);
  });

  it('does not mutate input arrays', () => {
    const time = new Float64Array([0.25, 0.5]);
    const conc = new Float64Array([1.85, 1.39]);
    const blqMask = new Uint8Array( [0, 0]);
    const tCopy = time.slice();
    const cCopy = conc.slice();
    insertC0(time, conc, blqMask);
    expect(Array.from(time)).toEqual(Array.from(tCopy));
    expect(Array.from(conc)).toEqual(Array.from(cCopy));
  });

  it('returns null when c0 estimation fails entirely', () => {
    // No measurable points at all.
    const time = new Float64Array([0.25, 0.5]);
    const conc = new Float64Array([0.04, 0.03]);
    const blqMask = new Uint8Array( [1, 1]);
    // Default chain ends with set0=0 → not null. To force null we restrict
    // methods to logslope only on data that doesn't allow it.
    const r = insertC0(time, conc, blqMask, {methods: ['logslope']});
    expect(r).toBeNull();
  });
});

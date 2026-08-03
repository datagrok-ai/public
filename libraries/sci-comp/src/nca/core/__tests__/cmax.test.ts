import {findCmax} from '../cmax';

describe('findCmax', () => {
  it('finds the unique maximum on a typical absorption profile', () => {
    const time = new Float64Array([0, 0.5, 1, 2, 4, 8]);
    const conc = new Float64Array([0, 1.5, 3, 2, 1, 0.5]);
    const blqMask = new Uint8Array( [0, 0, 0, 0, 0, 0]);
    const r = findCmax(time, conc, blqMask)!;
    expect(r.cmax).toBe(3);
    expect(r.tmax).toBe(1);
    expect(r.cmaxIdx).toBe(2);
  });

  it('PKNCA convention: returns the FIRST occurrence at tied maxima', () => {
    // Two points share Cmax = 3. PKNCA reports the earlier one.
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([2, 3, 3, 1]);
    const blqMask = new Uint8Array( [0, 0, 0, 0]);
    const r = findCmax(time, conc, blqMask)!;
    expect(r.cmax).toBe(3);
    expect(r.tmax).toBe(1);
    expect(r.cmaxIdx).toBe(1);
  });

  it('first occurrence holds even when the tie is at the very start', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([3, 3, 2]);
    const blqMask = new Uint8Array( [0, 0, 0]);
    const r = findCmax(time, conc, blqMask)!;
    expect(r.cmaxIdx).toBe(0);
  });

  it('skips BLQ points', () => {
    // Index 0 has the largest raw value but is BLQ; the actual peak among
    // measurable points is at index 2.
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([10, 3, 5, 2]);
    const blqMask = new Uint8Array( [1, 0, 0, 0]);
    const r = findCmax(time, conc, blqMask)!;
    expect(r.cmax).toBe(5);
    expect(r.tmax).toBe(2);
    expect(r.cmaxIdx).toBe(2);
  });

  it('skips NaN concentrations (e.g. from the `missing` BLQ rule)', () => {
    // Same as above but with NaN at index 0 — also has to be ignored even
    // though blqMask[0] = 0 (the BLQ rule already replaced the value).
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([NaN, 3, 5, 2]);
    const blqMask = new Uint8Array( [0, 0, 0, 0]);
    const r = findCmax(time, conc, blqMask)!;
    expect(r.cmax).toBe(5);
    expect(r.cmaxIdx).toBe(2);
  });

  describe('edge cases', () => {
    it('returns null when every point is BLQ', () => {
      const time = new Float64Array([0, 1, 2]);
      const conc = new Float64Array([0.04, 0.03, 0.02]);
      const blqMask = new Uint8Array( [1, 1, 1]);
      expect(findCmax(time, conc, blqMask)).toBeNull();
    });

    it('returns null when the input is empty', () => {
      expect(findCmax(new Float64Array(0), new Float64Array(0), new Uint8Array(0))).toBeNull();
    });

    it('returns null when every measurable conc is NaN', () => {
      const time = new Float64Array([0, 1]);
      const conc = new Float64Array([NaN, NaN]);
      const blqMask = new Uint8Array( [0, 0]);
      expect(findCmax(time, conc, blqMask)).toBeNull();
    });

    it('handles a single measurable point', () => {
      const time = new Float64Array([2.5]);
      const conc = new Float64Array([7.0]);
      const blqMask = new Uint8Array( [0]);
      const r = findCmax(time, conc, blqMask)!;
      expect(r.cmax).toBe(7.0);
      expect(r.tmax).toBe(2.5);
      expect(r.cmaxIdx).toBe(0);
    });

    it('returns null on a single BLQ point', () => {
      const time = new Float64Array([2.5]);
      const conc = new Float64Array([0.04]);
      const blqMask = new Uint8Array( [1]);
      expect(findCmax(time, conc, blqMask)).toBeNull();
    });

    it('returns the first index when all measurable values are equal', () => {
      const time = new Float64Array([0, 1, 2]);
      const conc = new Float64Array([2, 2, 2]);
      const blqMask = new Uint8Array( [0, 0, 0]);
      const r = findCmax(time, conc, blqMask)!;
      expect(r.cmaxIdx).toBe(0);
      expect(r.tmax).toBe(0);
    });

    it('handles negative concentrations (defensive — not expected in NCA inputs)', () => {
      const time = new Float64Array([0, 1, 2]);
      const conc = new Float64Array([-3, -1, -2]);
      const blqMask = new Uint8Array( [0, 0, 0]);
      const r = findCmax(time, conc, blqMask)!;
      expect(r.cmax).toBe(-1);
      expect(r.cmaxIdx).toBe(1);
    });
  });
});

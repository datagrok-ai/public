import {augmentProfile} from '../augment';
import {estimateC0Detailed, C0_DEFAULT_METHODS} from '../c0';
import {ROUTE_IV_BOLUS, ROUTE_IV_INFUSION, ROUTE_PO} from '../types';
import type {ProfileInputs, BlqStrategy, RouteCode} from '../types';

const SET_ZERO: BlqStrategy = {
  preFirstMeasurable: 'set-zero', embedded: 'set-zero',
  afterLast: 'set-zero', consecutiveAfterLast: 'set-zero',
};
const HALF_LLOQ: BlqStrategy = {
  preFirstMeasurable: 'set-half-lloq', embedded: 'set-half-lloq',
  afterLast: 'set-half-lloq', consecutiveAfterLast: 'set-half-lloq',
};
const EXCLUDE: BlqStrategy = {
  preFirstMeasurable: 'exclude', embedded: 'exclude',
  afterLast: 'exclude', consecutiveAfterLast: 'exclude',
};

function inputs(
  time: number[], conc: number[], route: RouteCode, blq?: number[],
): ProfileInputs {
  return {
    time: Float64Array.from(time),
    conc: Float64Array.from(conc),
    blqMask: blq === undefined ? new Uint8Array(time.length) : Uint8Array.from(blq),
    lloq: 0.05,
    dose: 25, doseUnits: 'mg', concentrationUnits: 'mg/L', timeUnits: 'h',
    route, infusionDuration: route === ROUTE_IV_INFUSION ? 1 : null, bodyWeight: null,
  };
}

// indometh subject 1 (02_indometh.csv) — the modal IV-bolus shape with no t=0 row.
const IND_T = [0.25, 0.5, 0.75, 1, 1.25, 2, 3, 4, 5, 6, 8];
const IND_C = [1.5, 0.94, 0.78, 0.48, 0.37, 0.19, 0.12, 0.11, 0.08, 0.07, 0.05];

describe('augmentProfile — sourceIndex contract (AC-U1.1.2)', () => {
  it('no augmentation → identity map', () => {
    // PO with a native t=0 row (theoph-1 shape: 0.74 at t=0).
    const a = augmentProfile(inputs([0, 0.25, 0.5, 1, 2], [0.74, 2.8, 6.6, 10.5, 9.7], ROUTE_PO), SET_ZERO)!;
    expect(a).not.toBeNull();
    expect(Array.from(a.sourceIndex)).toEqual([0, 1, 2, 3, 4]);
    expect(a.time.length).toBe(5); // same length as input, no synthetic slot
    expect(a.c0).toBeNull();
    expect(a.cmaxIdx).toBe(3);
  });

  it('prepend (extravascular, no t=0 row) → [-1, 0..n-1]', () => {
    const a = augmentProfile(inputs([0.5, 1, 2, 4], [1, 2, 1.5, 0.5], ROUTE_PO), SET_ZERO)!;
    expect(Array.from(a.sourceIndex)).toEqual([-1, 0, 1, 2, 3]);
    expect(a.time[0]).toBe(0);
    expect(a.conc[0]).toBe(0);
    expect(a.cmaxIdx).toBe(2); // observed peak at input index 1, shifted by the prepend
    expect(a.observedCmax.cmaxIdx).toBe(1);
    expect(a.c0).toBeNull();
  });

  it('prepend (IV bolus, no t=0 row) → [-1, 0..n-1] with (0, c0) in slot 0', () => {
    const a = augmentProfile(inputs(IND_T, IND_C, ROUTE_IV_BOLUS), SET_ZERO)!;
    expect(Array.from(a.sourceIndex)).toEqual([-1, ...IND_T.map((_, i) => i)]);
    expect(a.time[0]).toBe(0);
    expect(a.c0).not.toBeNull();
    expect(a.c0!.method).toBe('logslope');
    expect(a.c0!.replacedDoseTimeRow).toBe(false);
    expect(a.conc[0]).toBe(a.c0!.value);
    expect(a.c0!.value).toBeGreaterThan(IND_C[0]);
    expect(a.cmaxIdx).toBe(0); // the inserted c0 is the augmented peak
    // The REPORTED peak is the observed one, never the inserted c0.
    expect(a.observedCmax.cmax).toBe(1.5);
    expect(a.observedCmax.tmax).toBe(0.25);
  });

  it('REPLACE (IV bolus, non-positive t=0 row) → [-1, 1..n-1], input row 0 absent', () => {
    const a = augmentProfile(inputs([0, ...IND_T], [0, ...IND_C], ROUTE_IV_BOLUS), SET_ZERO)!;
    expect(a.time.length).toBe(IND_T.length + 1); // same length as the input
    expect(Array.from(a.sourceIndex)).toEqual([-1, ...IND_T.map((_, i) => i + 1)]);
    expect(a.sourceIndex).not.toContain(0);
    expect(a.c0!.replacedDoseTimeRow).toBe(true);
    expect(a.c0!.method).toBe('logslope');
    expect(a.conc[0]).toBe(a.c0!.value);
    // Byte-identical to the no-row profile in augmented space.
    const noRow = augmentProfile(inputs(IND_T, IND_C, ROUTE_IV_BOLUS), SET_ZERO)!;
    expect(Array.from(a.time)).toEqual(Array.from(noRow.time));
    expect(Array.from(a.conc)).toEqual(Array.from(noRow.conc));
    expect(a.c0!.value).toBe(noRow.c0!.value);
  });

  it('REPLACE (IV bolus, BLQ-FLAGGED t=0 row) — under every rule the row is absent (AD-3)', () => {
    for (const blq of [SET_ZERO, HALF_LLOQ, EXCLUDE]) {
      const a = augmentProfile(
        inputs([0, ...IND_T], [0, ...IND_C], ROUTE_IV_BOLUS, [1, ...IND_T.map(() => 0)]), blq)!;
      expect(a.c0!.replacedDoseTimeRow).toBe(true);
      expect(a.c0!.method).toBe('logslope');
      // Never the substitute: under set-half-lloq the row would read LLOQ/2 = 0.025.
      expect(a.conc[0]).not.toBe(0.025);
      expect(Array.from(a.sourceIndex)).toEqual([-1, ...IND_T.map((_, i) => i + 1)]);
    }
  });
});

describe('augmentProfile — boundary conditions (AC-U1.1.3)', () => {
  it('empty profile → null', () => {
    expect(augmentProfile(inputs([], [], ROUTE_PO), SET_ZERO)).toBeNull();
  });

  it('all-BLQ profile → null', () => {
    expect(augmentProfile(inputs([0, 1, 2], [0.01, 0.01, 0.01], ROUTE_IV_BOLUS, [1, 1, 1]), SET_ZERO))
      .toBeNull();
  });

  it('native positive t=0 on IV bolus → no augmentation, c0.method = "observed"', () => {
    const a = augmentProfile(inputs([0, 0.5, 1, 2, 4], [5, 4, 3.2, 2, 0.8], ROUTE_IV_BOLUS), SET_ZERO)!;
    expect(Array.from(a.sourceIndex)).toEqual([0, 1, 2, 3, 4]);
    expect(a.c0).toEqual({value: 5, method: 'observed', replacedDoseTimeRow: false});
    expect(a.cmaxIdx).toBe(0);
  });

  it('native 0.74 at t=0 on PO (theoph-1 shape) → no augmentation', () => {
    const a = augmentProfile(inputs([0, 0.25, 0.57, 1.12], [0.74, 2.84, 6.57, 10.5], ROUTE_PO), SET_ZERO)!;
    expect(Array.from(a.sourceIndex)).toEqual([0, 1, 2, 3]);
    expect(a.time.length).toBe(4);
  });

  it('native conc=0 at t=0 on PO / IV infusion → kept (a valid pre-dose observation), no prepend', () => {
    for (const route of [ROUTE_PO, ROUTE_IV_INFUSION]) {
      const a = augmentProfile(inputs([0, 1, 2, 4], [0, 3, 2, 1], route), SET_ZERO)!;
      expect(Array.from(a.sourceIndex)).toEqual([0, 1, 2, 3]);
      expect(a.c0).toBeNull();
    }
  });

  it('IV infusion is NOT route-gated like bolus: a t=0 zero is kept, not replaced', () => {
    const a = augmentProfile(inputs([0, 0.5, 1, 2, 4], [0, 4, 3.2, 2, 0.8], ROUTE_IV_INFUSION), SET_ZERO)!;
    expect(a.sourceIndex[0]).toBe(0);
    expect(a.conc[0]).toBe(0);
  });
});

describe('augmentProfile — the IV-bolus gate reads the EFFECTIVE mask + value, never the drop set (AC-U1.2.2)', () => {
  it('unflagged 0 at t=0 and flagged BLQ at t=0 both replace; a positive value never does', () => {
    const zero = augmentProfile(inputs([0, ...IND_T], [0, ...IND_C], ROUTE_IV_BOLUS), SET_ZERO)!;
    const flagged = augmentProfile(
      inputs([0, ...IND_T], [0.9, ...IND_C], ROUTE_IV_BOLUS, [1, ...IND_T.map(() => 0)]), SET_ZERO)!;
    const positive = augmentProfile(inputs([0, ...IND_T], [0.9, ...IND_C], ROUTE_IV_BOLUS), SET_ZERO)!;
    expect(zero.c0!.replacedDoseTimeRow).toBe(true);
    expect(flagged.c0!.replacedDoseTimeRow).toBe(true);
    expect(positive.c0!.replacedDoseTimeRow).toBe(false);
    expect(positive.c0!.method).toBe('observed');
    expect(positive.c0!.value).toBe(0.9);
  });

  it('a non-finite t=0 value is absent too', () => {
    const a = augmentProfile(inputs([0, ...IND_T], [NaN, ...IND_C], ROUTE_IV_BOLUS), SET_ZERO)!;
    expect(a.c0!.replacedDoseTimeRow).toBe(true);
    expect(Number.isFinite(a.conc[0])).toBe(true);
  });
});

describe('augmentProfile — only index 0 is the dose-time row (AC-U1.2.5)', () => {
  it('a second t=0 row is an ordinary observation', () => {
    // Two t=0 rows on IV bolus: row 0 is the dose-time row (non-positive →
    // replaced); row 1 (a replicate at t=0) stays as an ordinary observation
    // and keeps its source index.
    const a = augmentProfile(inputs([0, 0, ...IND_T], [0, 2.1, ...IND_C], ROUTE_IV_BOLUS), SET_ZERO)!;
    expect(a.c0!.replacedDoseTimeRow).toBe(true);
    expect(Array.from(a.sourceIndex)).toEqual([-1, 1, ...IND_T.map((_, i) => i + 2)]);
    expect(a.time[1]).toBe(0);
    expect(a.conc[1]).toBe(2.1);
  });
});

describe('augmentProfile — two masks (AD-2): effective BLQ vs the drop set', () => {
  const MISSING: BlqStrategy = {
    preFirstMeasurable: 'missing', embedded: 'missing',
    afterLast: 'missing', consecutiveAfterLast: 'missing',
  };
  const T = [0.5, 1, 2, 4, 6];
  const C = [0.01, 2, 1.5, 0.5, 0.01];
  const M = [1, 0, 0, 0, 1];

  it('exclude: both masks carry the input mask ∪ excluded, synthetic slot unmasked', () => {
    const a = augmentProfile(inputs(T, C, ROUTE_PO, M), EXCLUDE)!;
    expect(Array.from(a.blqMask)).toEqual([0, 1, 0, 0, 0, 1]);
    expect(Array.from(a.dropMask)).toEqual([0, 1, 0, 0, 0, 1]);
    expect(Array.from(a.blqApplied.excluded)).toEqual([0, 4]); // INPUT index space
  });

  it('missing: NaN lands in the drop set exactly like exclude', () => {
    const a = augmentProfile(inputs(T, C, ROUTE_PO, M), MISSING)!;
    expect(Array.from(a.dropMask)).toEqual([0, 1, 0, 0, 0, 1]);
    expect(Number.isNaN(a.conc[1])).toBe(true);
    expect(a.blqApplied.excluded.length).toBe(0);
  });

  it('set-zero / set-half-lloq: the BLQ mask still flags the point, the drop set does NOT', () => {
    for (const [blq, value] of [[SET_ZERO, 0], [HALF_LLOQ, 0.025]] as const) {
      const a = augmentProfile(inputs(T, C, ROUTE_PO, M), blq)!;
      expect(Array.from(a.blqMask)).toEqual([0, 1, 0, 0, 0, 1]); // Cmax / Tlag semantics
      expect(Array.from(a.dropMask)).toEqual([0, 0, 0, 0, 0, 0]); // integrates the substitute
      expect(a.conc[1]).toBe(value);
      expect(a.conc[5]).toBe(value);
      // The observed peak still skips the flagged points.
      expect(a.observedCmax.cmax).toBe(2);
    }
  });

  it('EV t=0: a substituted BLQ at the dose time is a KEPT dose-time value — no (0, 0) prepend (AC-U2.1.2)', () => {
    // compute-nca.test.ts "one BLQ at the start": BLQ at index 0 under set-zero.
    const t = [0, 0.25, 0.5, 1, 2, 4, 8, 12];
    const c = [0, 0.5, 1.0, 0.8, 0.5, 0.3, 0.1, 0.05];
    const kept = augmentProfile(inputs(t, c, ROUTE_PO, [1, 0, 0, 0, 0, 0, 0, 0]), SET_ZERO)!;
    expect(Array.from(kept.sourceIndex)).not.toContain(-1);
    expect(kept.time.length).toBe(t.length);
    // exclude / missing at t=0 → dropped → prepend as before.
    for (const blq of [EXCLUDE, MISSING]) {
      const dropped = augmentProfile(inputs(t, c, ROUTE_PO, [1, 0, 0, 0, 0, 0, 0, 0]), blq)!;
      expect(dropped.sourceIndex[0]).toBe(-1);
      expect(dropped.time.length).toBe(t.length + 1);
    }
  });
});

describe('estimateC0Detailed — default chain never null given a post-dose measurable point (AD-6)', () => {
  it.each([
    ['decaying pair → logslope', [0.25, 0.5, 1], [1.5, 0.94, 0.48], [0, 0, 0], 'logslope'],
    ['rising pair → c1', [0.25, 0.5, 1], [1.0, 1.2, 0.8], [0, 0, 0], 'c1'],
    ['single post-dose point → c1', [0.25], [1.0], [0], 'c1'],
    ['second point BLQ → c1', [0.25, 0.5], [1.0, 0.01], [0, 1], 'c1'],
  ])('%s', (_label, t, c, m, method) => {
    const r = estimateC0Detailed(Float64Array.from(t), Float64Array.from(c), Uint8Array.from(m));
    expect(r).not.toBeNull();
    expect(r!.method).toBe(method);
    expect(Number.isFinite(r!.value)).toBe(true);
  });

  it('the chain ends in set0, so even an all-zero post-dose profile answers', () => {
    const r = estimateC0Detailed(Float64Array.from([1, 2]), Float64Array.from([0, 0]), new Uint8Array(2));
    expect(r).not.toBeNull();
    // c1 accepts a finite 0 (PKNCA semantics) — the chain answers before set0.
    expect(['c1', 'set0']).toContain(r!.method);
  });

  it('null is reachable only through a caller-supplied chain', () => {
    const r = estimateC0Detailed(
      Float64Array.from([1, 2]), Float64Array.from([1.0, 1.2]), new Uint8Array(2), {methods: ['logslope']});
    expect(r).toBeNull();
    expect(C0_DEFAULT_METHODS[C0_DEFAULT_METHODS.length - 1]).toBe('set0');
  });
});

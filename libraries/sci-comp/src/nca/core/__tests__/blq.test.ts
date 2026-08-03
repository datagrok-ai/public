import {applyBlqStrategy} from '../blq';
import type {BlqRule, BlqStrategy} from '../types';

/** Build a BlqStrategy that applies the same rule to every phase. */
function uniform(rule: BlqRule): BlqStrategy {
  return {
    preFirstMeasurable: rule,
    embedded: rule,
    afterLast: rule,
    consecutiveAfterLast: rule,
  };
}

/** A 5-point profile with one BLQ in each of the four phases.
 *
 * Index:        0    1    2    3    4    5    6    7
 * blqMask:      1    0    1    0    1    1    0    1    1
 *               ^    ^    ^    ^                   ^    ^
 *               pre  M    emb  M    (none)         M    aft+cons
 *
 * Concretely we lay out:
 *   i=0 BLQ  → preFirstMeasurable
 *   i=1 measurable
 *   i=2 BLQ  → embedded
 *   i=3 measurable (last)
 *   i=4 BLQ  → afterLast
 *   i=5 BLQ  → consecutiveAfterLast
 */
const ALL_PHASES = {
  conc: new Float64Array([0.04, 1.0, 0.04, 0.5, 0.04, 0.04]),
  blqMask: new Uint8Array( [1, 0, 1, 0, 1, 1]),
  lloq: 0.05,
};
const PRE_IDX = 0;
const MEASURABLE_1 = 1;
const EMBEDDED_IDX = 2;
const MEASURABLE_2 = 3;
const AFTER_IDX = 4;
const CONSEC_IDX = 5;

describe('applyBlqStrategy', () => {
  // ────────────────────────────────────────────────────────────────────
  // 4 phases × 4 rules = 16 explicit combinations
  // ────────────────────────────────────────────────────────────────────

  describe('preFirstMeasurable phase', () => {
    it('set-zero replaces conc with 0', () => {
      const s: BlqStrategy = {...uniform('exclude'), preFirstMeasurable: 'set-zero'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(r.conc[PRE_IDX]).toBe(0);
      expect(Array.from(r.excluded)).not.toContain(PRE_IDX);
    });

    it('set-half-lloq replaces conc with lloq/2', () => {
      const s: BlqStrategy = {...uniform('exclude'), preFirstMeasurable: 'set-half-lloq'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, 0.05, 0, s);
      expect(r.conc[PRE_IDX]).toBe(0.025);
    });

    it('exclude lists the index in `excluded`', () => {
      const s: BlqStrategy = {...uniform('set-zero'), preFirstMeasurable: 'exclude'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(Array.from(r.excluded)).toContain(PRE_IDX);
      expect(r.conc[PRE_IDX]).toBe(0.04); // value left as-is
    });

    it('missing replaces conc with NaN', () => {
      const s: BlqStrategy = {...uniform('set-zero'), preFirstMeasurable: 'missing'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(Number.isNaN(r.conc[PRE_IDX])).toBe(true);
    });
  });

  describe('embedded phase', () => {
    it('set-zero', () => {
      const s: BlqStrategy = {...uniform('exclude'), embedded: 'set-zero'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(r.conc[EMBEDDED_IDX]).toBe(0);
    });

    it('set-half-lloq', () => {
      const s: BlqStrategy = {...uniform('exclude'), embedded: 'set-half-lloq'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, 0.05, 0, s);
      expect(r.conc[EMBEDDED_IDX]).toBe(0.025);
    });

    it('exclude', () => {
      const s: BlqStrategy = {...uniform('set-zero'), embedded: 'exclude'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(Array.from(r.excluded)).toContain(EMBEDDED_IDX);
    });

    it('missing', () => {
      const s: BlqStrategy = {...uniform('set-zero'), embedded: 'missing'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(Number.isNaN(r.conc[EMBEDDED_IDX])).toBe(true);
    });
  });

  describe('afterLast phase (first BLQ after last measurable)', () => {
    it('set-zero', () => {
      const s: BlqStrategy = {...uniform('exclude'), afterLast: 'set-zero'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(r.conc[AFTER_IDX]).toBe(0);
    });

    it('set-half-lloq', () => {
      const s: BlqStrategy = {...uniform('exclude'), afterLast: 'set-half-lloq'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, 0.05, 0, s);
      expect(r.conc[AFTER_IDX]).toBe(0.025);
    });

    it('exclude', () => {
      const s: BlqStrategy = {...uniform('set-zero'), afterLast: 'exclude'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(Array.from(r.excluded)).toContain(AFTER_IDX);
    });

    it('missing', () => {
      const s: BlqStrategy = {...uniform('set-zero'), afterLast: 'missing'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(Number.isNaN(r.conc[AFTER_IDX])).toBe(true);
    });
  });

  describe('consecutiveAfterLast phase (second+ BLQ in tail)', () => {
    it('set-zero', () => {
      const s: BlqStrategy = {...uniform('exclude'), consecutiveAfterLast: 'set-zero'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(r.conc[CONSEC_IDX]).toBe(0);
    });

    it('set-half-lloq', () => {
      const s: BlqStrategy = {...uniform('exclude'), consecutiveAfterLast: 'set-half-lloq'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, 0.05, 0, s);
      expect(r.conc[CONSEC_IDX]).toBe(0.025);
    });

    it('exclude', () => {
      const s: BlqStrategy = {...uniform('set-zero'), consecutiveAfterLast: 'exclude'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(Array.from(r.excluded)).toContain(CONSEC_IDX);
    });

    it('missing', () => {
      const s: BlqStrategy = {...uniform('set-zero'), consecutiveAfterLast: 'missing'};
      const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0, s);
      expect(Number.isNaN(r.conc[CONSEC_IDX])).toBe(true);
    });
  });

  // ────────────────────────────────────────────────────────────────────
  // Cross-cutting behaviour
  // ────────────────────────────────────────────────────────────────────

  it('does not mutate the input concentration array', () => {
    const original = ALL_PHASES.conc.slice();
    applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0,
      uniform('set-zero'));
    expect(Array.from(ALL_PHASES.conc)).toEqual(Array.from(original));
  });

  it('leaves measurable points untouched', () => {
    const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, ALL_PHASES.lloq, 0,
      uniform('missing'));
    expect(r.conc[MEASURABLE_1]).toBe(1.0);
    expect(r.conc[MEASURABLE_2]).toBe(0.5);
  });

  it('per-row lloq array is honoured by set-half-lloq', () => {
    const lloqArr = new Float64Array([0.10, 0, 0.20, 0, 0.30, 0.40]);
    const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, lloqArr, 0,
      uniform('set-half-lloq'));
    expect(r.conc[PRE_IDX]).toBe(0.05);
    expect(r.conc[EMBEDDED_IDX]).toBe(0.10);
    expect(r.conc[AFTER_IDX]).toBe(0.15);
    expect(r.conc[CONSEC_IDX]).toBe(0.20);
  });

  it('different phases can use different rules independently', () => {
    const s: BlqStrategy = {
      preFirstMeasurable: 'set-zero',
      embedded: 'set-half-lloq',
      afterLast: 'exclude',
      consecutiveAfterLast: 'missing',
    };
    const r = applyBlqStrategy(ALL_PHASES.conc, ALL_PHASES.blqMask, 0.05, 0, s);
    expect(r.conc[PRE_IDX]).toBe(0);
    expect(r.conc[EMBEDDED_IDX]).toBe(0.025);
    expect(Array.from(r.excluded)).toEqual([AFTER_IDX]);
    expect(Number.isNaN(r.conc[CONSEC_IDX])).toBe(true);
  });

  // ────────────────────────────────────────────────────────────────────
  // Edge cases
  // ────────────────────────────────────────────────────────────────────

  describe('edge cases', () => {
    it('all points BLQ → every point treated as preFirstMeasurable', () => {
      const conc = new Float64Array([0.04, 0.03, 0.02]);
      const blqMask = new Uint8Array( [1, 1, 1]);
      const s: BlqStrategy = {
        preFirstMeasurable: 'set-zero',
        embedded: 'set-half-lloq',
        afterLast: 'set-half-lloq',
        consecutiveAfterLast: 'set-half-lloq',
      };
      const r = applyBlqStrategy(conc, blqMask, 0.05, 0, s);
      expect(Array.from(r.conc)).toEqual([0, 0, 0]);
    });

    it('no BLQ → output identical to input, nothing excluded', () => {
      const conc = new Float64Array([1, 2, 3, 4]);
      const blqMask = new Uint8Array( [0, 0, 0, 0]);
      const r = applyBlqStrategy(conc, blqMask, 0.05, 0, uniform('set-zero'));
      expect(Array.from(r.conc)).toEqual([1, 2, 3, 4]);
      expect(r.excluded.length).toBe(0);
    });

    it('single measurable point → no BLQ phase activated', () => {
      const conc = new Float64Array([2.5]);
      const blqMask = new Uint8Array( [0]);
      const r = applyBlqStrategy(conc, blqMask, 0.05, 0, uniform('set-zero'));
      expect(Array.from(r.conc)).toEqual([2.5]);
      expect(r.excluded.length).toBe(0);
    });

    it('single BLQ point → preFirstMeasurable phase', () => {
      const conc = new Float64Array([0.04]);
      const blqMask = new Uint8Array( [1]);
      const r = applyBlqStrategy(conc, blqMask, 0.05, 0,
        {...uniform('set-zero'), preFirstMeasurable: 'set-half-lloq'});
      expect(r.conc[0]).toBe(0.025);
    });

    it('only one measurable, BLQ before and after → preFirst + afterLast', () => {
      const conc = new Float64Array([0.04, 1.0, 0.04, 0.04]);
      const blqMask = new Uint8Array( [1, 0, 1, 1]);
      const s: BlqStrategy = {
        preFirstMeasurable: 'set-zero',
        embedded: 'missing', // not used here
        afterLast: 'set-half-lloq',
        consecutiveAfterLast: 'exclude',
      };
      const r = applyBlqStrategy(conc, blqMask, 0.05, 0, s);
      expect(r.conc[0]).toBe(0); // preFirst
      expect(r.conc[1]).toBe(1.0); // measurable, untouched
      expect(r.conc[2]).toBe(0.025); // afterLast
      expect(Array.from(r.excluded)).toEqual([3]); // consecutiveAfterLast
    });

    it('empty input → empty output, empty excluded', () => {
      const r = applyBlqStrategy(new Float64Array(0), new Uint8Array(0), 0.05, 0,
        uniform('set-zero'));
      expect(r.conc.length).toBe(0);
      expect(r.excluded.length).toBe(0);
    });
  });
});

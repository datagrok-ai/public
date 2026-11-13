/**
 * Unit tests for the DG missing-value boundary helpers
 * (`adapter/adapter-utils.ts`: `f64Column`, `numAt`, `finiteMinMax`).
 *
 * **Rule 16 — these assert against a REAL `DG.Column`, not a mirror of the
 * helper's own ternary.** A mirror test would prove nothing about the platform,
 * and the platform's behaviour is the entire subject: `isNone(rawNaN) === false`
 * while `isNone(FLOAT_NULL) === true`, and `Column.get()` returns the raw
 * sentinel even on a null cell. Every assertion below goes through `isNone` /
 * `get` / `getRawData` on a constructed column.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import {f64Column, numAt, finiteMinMax} from '../adapter/adapter-utils';

category('DG null sentinel: f64Column (F1 write boundary)', () => {
  test('AC-1: NaN becomes a platform null; a real value does not', async () => {
    const col = f64Column('X', [1, NaN, 3]);
    expect(col.isNone(1), true, 'NaN element reads isNone');
    expect(col.isNone(0), false, 'real 1 does not read isNone');
    expect(col.isNone(2), false, 'real 3 does not read isNone');
    // The raw storage really holds the platform sentinel, not NaN. This is the
    // measured half: a raw NaN would have left isNone false (probe PR-1).
    const raw = col.getRawData() as Float64Array;
    expect(raw[1] === DG.FLOAT_NULL, true, 'raw storage holds DG.FLOAT_NULL');
    expect(Number.isNaN(raw[1]), false, 'raw storage is NOT IEEE NaN');
  });

  test('AC-1b: a raw NaN written WITHOUT the helper is NOT null (the defect)', async () => {
    // The control that gives AC-1 its meaning. Same array, platform constructor.
    const col = DG.Column.fromFloat64Array('X', Float64Array.from([1, NaN, 3]));
    expect(col.isNone(1), false,
      'raw NaN is a PRESENT value to DG — this is what f64Column exists to fix');
  });

  test('AC-2: an already-FLOAT_NULL value passes through and still reads isNone', async () => {
    const col = f64Column('X', [DG.FLOAT_NULL, 2, DG.FLOAT_NULL]);
    expect(col.isNone(0), true, 'pre-existing FLOAT_NULL still null');
    expect(col.isNone(2), true, 'pre-existing FLOAT_NULL still null');
    expect(col.isNone(1), false, 'real value untouched');
    const raw = col.getRawData() as Float64Array;
    expect(raw[0] === DG.FLOAT_NULL, true, 'bits unchanged (no double-translation)');
  });

  test('AC-3: ±Infinity passes through UNCHANGED — it is not translated to null', async () => {
    // Regression pin for the reverted first draft of AC-3. The predicate is
    // Number.isNaN, NOT !Number.isFinite: bioavailability can emit +Infinity for
    // a zero-exposure IV reference, that value is aggregated as PRESENT today,
    // and nulling it here would change N_F_PCT and manufacture a false
    // "no IV reference" disclosure. Filed as its own gap; not absorbed here.
    const col = f64Column('X', [Infinity, -Infinity, 5]);
    expect(col.isNone(0), false, '+Infinity is a value, not a null');
    expect(col.isNone(1), false, '-Infinity is a value, not a null');
    const raw = col.getRawData() as Float64Array;
    expect(raw[0] === Infinity, true, '+Infinity survives bit-exact');
    expect(raw[1] === -Infinity, true, '-Infinity survives bit-exact');
  });

  test('AC-4: the helper COPIES — the caller-owned array is never mutated', async () => {
    // Load-bearing: fromFloat64Array ADOPTS its argument as zero-copy storage,
    // and results-writer passes the controller's own `c0s` array.
    const owned = Float64Array.from([1, NaN, 3]);
    const col = f64Column('X', owned);
    expect(Number.isNaN(owned[1]), true, 'caller array still holds NaN');
    expect(owned[1] === DG.FLOAT_NULL, false, 'caller array NOT rewritten in place');
    expect(col.isNone(1), true, 'the column, however, is null there');
    // And the column does not alias the caller's buffer.
    owned[0] = 99;
    expect(col.get(0) as number, 1, 'column storage is independent of the caller array');
  });

  test('AC-4b: a null / undefined / non-numeric element THROWS, never coerces', async () => {
    // Float64Array.from([null]) yields 0 — on a semi-log concentration axis that
    // is a real point at zero rather than a gap. The helper must not decide it.
    let threw = false;
    try {
      f64Column('X', [1, null as unknown as number, 3]);
    } catch {
      threw = true;
    }
    expect(threw, true, 'null element throws');

    threw = false;
    try {
      f64Column('X', [1, undefined as unknown as number]);
    } catch {
      threw = true;
    }
    expect(threw, true, 'undefined element throws');

    threw = false;
    try {
      f64Column('X', ['7' as unknown as number]);
    } catch {
      threw = true;
    }
    expect(threw, true, 'string element throws');
  });

  test('AC-6: real values round-trip BIT-EXACT (NFR-04 parity is not perturbed)', async () => {
    // Equality on the raw backing array, not a tolerance — a PKNCA-parity value
    // must survive the boundary unchanged to the last bit.
    // 2.678934406368464e-34 is the next representable double ABOVE DG.FLOAT_NULL
    // (verified by bit increment, 2026-08-03) — the tightest possible check that
    // the helper keys on exact bits and not on a neighbourhood.
    const vals = [0, -0, 1e-300, 2.678934406368464e-34, 123.456789012345, 1.7976931348623157e308];
    const col = f64Column('X', vals);
    const raw = col.getRawData() as Float64Array;
    for (let i = 0; i < vals.length; i++)
      expect(raw[i] === vals[i], true, `element ${i} bit-exact (${vals[i]})`);
    expect(col.isNone(3), false, 'the value one ULP above FLOAT_NULL is still a value');
  });

  test('AC-6b: an empty array builds an empty column (no special-casing)', async () => {
    const col = f64Column('X', []);
    expect(col.length, 0, 'zero-length column');
  });
});

category('DG null sentinel: GRID TEXT of a null cell (RG-3 / DATA-GAP-40)', () => {
  // The migration's user-visible deliverable is a BLANK grid cell where the literal
  // text "NaN" used to be. PR-4 measured the grid-text path for Float64; PR-7
  // measured only the *Column* API for Int32 (`getString(INT_NULL) === ''`) while
  // live QA observed the GRID rendering `LAMBDAZ_NPOINTS` as the literal "NaN".
  // Two different platform paths — one measured, one inferred. This closes the gap
  // by measuring the Int32 grid path directly, which decides whether the F-11
  // contingency (emit LAMBDAZ_NPOINTS as Float64) needs to fire at all.
  test('RG-3: Float64 null, Int32 null and raw NaN — what does the GRID actually render?', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('F64', Float64Array.from([7.5, DG.FLOAT_NULL, NaN])),
      DG.Column.fromInt32Array('I32', Int32Array.from([11, DG.INT_NULL, DG.INT_NULL])),
    ]);
    const grid = DG.Grid.create(df);
    // `GridCell.cell.valueString` is the text the grid renders for the cell — the
    // grid-side path, distinct from `Column.getString` measured above.
    const txt = (col: string, row: number): string => `${grid.cell(col, row).cell.valueString}`;
    // Column-API baseline (already measured; restated so a platform change shows up
    // as a diff on the pair rather than on one of them).
    expect(df.col('F64')!.getString(1), '', 'Column API: Float64 null → empty string');
    expect(df.col('I32')!.getString(1), '', 'Column API: Int32 null → empty string');
    expect(df.col('F64')!.isNone(2), false, 'a raw NaN is still NOT null to DG');

    // THE MEASUREMENT. Recorded verbatim in the cycle state; a change here is a
    // platform behaviour change and should fail loudly rather than drift.
    const f64Null = txt('F64', 1);
    const i32Null = txt('I32', 1);
    const rawNan = txt('F64', 2);
    console.log(`[RG-3 grid text] F64-null="${f64Null}" I32-null="${i32Null}" rawNaN="${rawNan}"`);
    expect(f64Null, '', 'GRID: a Float64 null cell renders BLANK');
    expect(i32Null, '', 'GRID: an Int32 null cell renders BLANK — this is what RG-3 asked');
    expect(rawNan, 'NaN', 'GRID: a raw NaN renders the literal text "NaN" — the trigger symptom');
  });
});

category('DG null sentinel: numAt / finiteMinMax (F1 read boundary)', () => {
  test('numAt maps a platform null back to NaN, and a raw get() does not', async () => {
    const col = f64Column('X', [7.5, NaN]);
    expect(Number.isNaN(numAt(col, 1)), true, 'null cell reads back as NaN');
    expect(numAt(col, 0), 7.5, 'real cell reads back verbatim');
    // The measured trap (probe PR-2): get() returns the RAW sentinel on a null
    // cell, and that value is FINITE and POSITIVE. This is why every result-column
    // read must go through numAt.
    const bare = col.get(1) as number;
    expect(Number.isFinite(bare), true, 'bare get() on a null cell IS finite — the trap');
    expect(bare > 0, true, 'and positive, so `> 0` / `< threshold` guards fire');
  });

  test('finiteMinMax skips nulls; positiveOnly additionally skips <= 0', async () => {
    const col = f64Column('C', [5, NaN, 20, 0]);
    const all = finiteMinMax(col);
    expect(all !== null, true, 'range found');
    expect(all!.min, 0, 'default mode includes the real 0');
    expect(all!.max, 20, 'max unaffected by the null');
    const pos = finiteMinMax(col, true);
    expect(pos !== null, true, 'positive range found');
    expect(pos!.min, 5, 'positiveOnly excludes the real 0 (log axis has no place for it)');
    expect(pos!.max, 20, 'max unchanged');
  });

  test('finiteMinMax excludes FLOAT_NULL where a bare isFinite && v>0 would admit it', async () => {
    // The R2-3 regression, pinned. A column holding only the sentinel must yield
    // NO range — the naive guard would return min = 2.68e-34 and stretch a
    // semi-log axis ~40 decades.
    const col = f64Column('C', [NaN, NaN]);
    expect(finiteMinMax(col, true), null, 'an all-null column yields no range');
    // And the naive predicate, for contrast — this is what used to run.
    const bare = col.get(0) as number;
    expect(Number.isFinite(bare) && bare > 0, true,
      'the guard this replaced would have admitted the sentinel');
  });

  test('finiteMinMax returns null for an absent column and for an all-non-finite one', async () => {
    expect(finiteMinMax(null), null, 'absent column → null');
    const inf = f64Column('C', [Infinity, -Infinity]);
    expect(finiteMinMax(inf), null, '±Infinity is not a range bound');
  });
});

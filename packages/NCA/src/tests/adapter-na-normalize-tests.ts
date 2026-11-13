/**
 * Adapter — NONMEM NA-marker normalization category.
 *
 * @see src/adapter/na-normalize.ts (implementation under test)
 *
 * Coverage:
 * - Round-trip through the real platform CSV parser (`DG.DataFrame.fromCsv`):
 *   a numeric column carrying `.`/empty markers is parsed STRING, then retyped
 *   FLOAT with markers → null. This is the actual symptom Ed reported, not a
 *   mirror-pattern (CLAUDE.md rule 16).
 * - Non-coercible string columns (text categories, mixed) are left alone.
 * - All-marker columns stay STRING (no phantom empty FLOAT column).
 * - Empty-string markers are treated as NA alongside `.`.
 * - Column order is preserved; converted-name list is accurate.
 *
 * The inline CSVs run through the same `DG.DataFrame.fromCsv` path that the
 * production launcher and `loadFixture` use, so these cover the real platform
 * parse → normalize round-trip. Integration assertions against the bundled
 * 03/06 NONMEM fixtures land with VAL-01, which bundles those files into
 * `files/` alongside their PKNCA reference JSONs.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {normalizeNumericNaMarkers} from '../adapter/na-normalize';

category('Adapter: NA-marker normalization', () => {
  // -------------------------------------------------------------------
  // Round-trip through the real platform CSV parser
  // -------------------------------------------------------------------

  test('NONMEM dot column: STRING on parse → FLOAT after normalize', async () => {
    // DV has `.` on the dose row, real values on observation rows — exactly the
    // NONMEM layout. The platform parses the mixed column as STRING.
    const df = DG.DataFrame.fromCsv(
      'ID,DV\n' +
      '1,.\n' +
      '1,235.16\n' +
      '1,286.29\n');
    expect(df.col('DV')!.type, DG.COLUMN_TYPE.STRING, 'precondition: DV parses STRING');

    const converted = normalizeNumericNaMarkers(df);
    expect(converted.includes('DV'), true, 'DV reported as converted');
    expect(df.col('DV')!.type, DG.COLUMN_TYPE.FLOAT, 'DV retyped FLOAT');

    const dv = df.col('DV')!;
    expect(dv.isNone(0), true, 'dose-row `.` → null');
    expect(dv.isNone(1), false, 'observation value present');
    expect(dv.get(1), 235.16, 'value preserved');
    expect(dv.get(2), 286.29, 'value preserved');
  });

  test('empty-string and dot markers both treated as NA', async () => {
    // A `.` forces the column STRING; an empty field (real empty cell, kept off
    // a blank line by the TAG column) must also normalize to null. The platform
    // handles a bare empty cell natively, so the marker that matters is `.` —
    // this asserts empty is still cleared once the column is being converted.
    const df = DG.DataFrame.fromCsv(
      'AMT,TAG\n' +
      '0.2441,a\n' +
      '.,b\n' +
      ',c\n' +
      '0.5,d\n');
    expect(df.col('AMT')!.type, DG.COLUMN_TYPE.STRING, 'precondition: `.` → STRING');
    const converted = normalizeNumericNaMarkers(df);
    expect(converted.includes('AMT'), true, 'AMT converted');
    const amt = df.col('AMT')!;
    expect(amt.type, DG.COLUMN_TYPE.FLOAT, 'AMT FLOAT');
    expect(amt.isNone(1), true, '`.` cell → null');
    expect(amt.isNone(2), true, 'empty cell → null');
    expect(amt.get(0), 0.2441, 'value preserved');
    expect(amt.get(3), 0.5, 'value preserved');
    expect(df.col('TAG')!.type, DG.COLUMN_TYPE.STRING, 'TAG stays STRING');
  });

  test('non-numeric string column left untouched (ROUTE / subject id)', async () => {
    const df = DG.DataFrame.fromCsv(
      'SubjectID,ROUTE\nR001,IV\nR002,PO\nR003,.\n');
    const converted = normalizeNumericNaMarkers(df);
    expect(converted.length, 0, 'nothing converted');
    expect(df.col('SubjectID')!.type, DG.COLUMN_TYPE.STRING, 'SubjectID stays STRING');
    // ROUTE has a `.` but also real text → not coercible, stays STRING.
    expect(df.col('ROUTE')!.type, DG.COLUMN_TYPE.STRING, 'ROUTE stays STRING');
  });

  test('all-marker column stays STRING (no phantom empty FLOAT)', async () => {
    // File 03's LLOQ is `.` for every row — there is no numeric content.
    const df = DG.DataFrame.fromCsv('ID,LLOQ\n1,.\n1,.\n1,.\n');
    const converted = normalizeNumericNaMarkers(df);
    expect(converted.includes('LLOQ'), false, 'all-`.` LLOQ not converted');
    expect(df.col('LLOQ')!.type, DG.COLUMN_TYPE.STRING, 'LLOQ stays STRING');
  });

  test('clean numeric column is already FLOAT — reported as no-op', async () => {
    const df = DG.DataFrame.fromCsv('conc\n10\n5\n2.5\n');
    expect(df.col('conc')!.type, DG.COLUMN_TYPE.FLOAT, 'clean column parses FLOAT');
    const converted = normalizeNumericNaMarkers(df);
    expect(converted.includes('conc'), false, 'already-FLOAT column not re-converted');
  });

  test('column order preserved across conversion', async () => {
    const df = DG.DataFrame.fromCsv('A,DV,C\nx,.,1\ny,2.0,2\n');
    const before = df.columns.names().join(',');
    normalizeNumericNaMarkers(df);
    const after = df.columns.names().join(',');
    expect(after, before, 'column order unchanged');
  });
});

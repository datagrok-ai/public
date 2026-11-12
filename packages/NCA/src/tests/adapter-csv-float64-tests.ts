/**
 * Adapter — Float64 CSV re-parse category (BUG-01).
 *
 * @see src/adapter/csv-float64.ts (implementation under test)
 *
 * Coverage:
 * - The exact BUG-01 symptom: a value DG parses through Float32 (`1.12` →
 *   `1.1200000047683716`) is recovered to the bit-identical Float64 after
 *   re-parse. Asserted against the real `DG.DataFrame.fromCsv` parser, not a
 *   mirror (CLAUDE.md rule 16).
 * - Empty cells in a FLOAT column → null (not 0).
 * - STRING columns are left untouched (only FLOAT columns re-parsed).
 * - Column order preserved; converted-name list accurate.
 * - Quoted fields with embedded commas tokenize correctly.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {reparseFloatColumnsToFloat64} from '../adapter/csv-float64';

category('Adapter: Float64 CSV re-parse', () => {
  test('BUG-01 symptom: Float32-drifted value recovered to exact Float64', async () => {
    const df = DG.DataFrame.fromCsv('Time\n1.12\n2.02\n');
    const before = df.col('Time')!.get(0) as number;
    // Precondition: DG's Float32 parse drifts 1.12.
    expect(before !== 1.12, true, `precondition: expected drift, got exact ${before}`);

    const converted = reparseFloatColumnsToFloat64(df, 'Time\n1.12\n2.02\n');
    expect(converted.includes('Time'), true, 'Time re-parsed');
    expect(df.col('Time')!.get(0), 1.12, 'value bit-identical to Number("1.12")');
    expect(df.col('Time')!.get(1), 2.02, 'second value exact');
  });

  test('empty cell in FLOAT column → null, not 0', async () => {
    // TAG keeps the empty cell off a blank line so the row survives parsing.
    const df = DG.DataFrame.fromCsv('conc,TAG\n1.12,a\n,b\n2.02,c\n');
    expect(df.col('conc')!.type, DG.COLUMN_TYPE.FLOAT, 'conc parses FLOAT');
    reparseFloatColumnsToFloat64(df, 'conc,TAG\n1.12,a\n,b\n2.02,c\n');
    const conc = df.col('conc')!;
    expect(conc.get(0), 1.12, 'value exact');
    expect(conc.isNone(1), true, 'empty cell → null');
    expect(conc.get(2), 2.02, 'value exact');
  });

  test('STRING column untouched (only FLOAT re-parsed)', async () => {
    const csv = 'SubjectID,conc\nR001,1.12\nR002,2.02\n';
    const df = DG.DataFrame.fromCsv(csv);
    const converted = reparseFloatColumnsToFloat64(df, csv);
    expect(converted.includes('SubjectID'), false, 'string column not re-parsed');
    expect(df.col('SubjectID')!.type, DG.COLUMN_TYPE.STRING, 'SubjectID stays STRING');
    expect(converted.includes('conc'), true, 'conc re-parsed');
  });

  test('column order preserved', async () => {
    const csv = 'A,conc,C\nx,1.12,1\ny,2.02,2\n';
    const df = DG.DataFrame.fromCsv(csv);
    const before = df.columns.names().join(',');
    reparseFloatColumnsToFloat64(df, csv);
    expect(df.columns.names().join(','), before, 'order unchanged');
  });

  test('quoted field with embedded comma tokenizes correctly', async () => {
    // The quoted label column must not desync the numeric column index.
    const csv = 'label,conc\n"alpha, beta",1.12\n"gamma",2.02\n';
    const df = DG.DataFrame.fromCsv(csv);
    reparseFloatColumnsToFloat64(df, csv);
    expect(df.col('conc')!.get(0), 1.12, 'conc value exact despite quoted comma');
    expect(df.col('conc')!.get(1), 2.02, 'second conc value exact');
  });
});

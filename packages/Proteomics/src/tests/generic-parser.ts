import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {SEMTYPE} from '../utils/proteomics-types';
import {
  detectDelimiter,
  autoSuggestProteinIdColumn,
  autoSuggestIntensityColumns,
  detectLog2Status,
  log2TransformColumns,
  copyAsLog2Columns,
  addPrimaryColumnIfNeeded,
} from '../parsers/shared-utils';

/** Builds a CSV or TSV string for testing. */
function makeGenericCsv(params: {headers: string[]; rows: string[][]; delimiter?: string}): string {
  const d = params.delimiter ?? ',';
  return [params.headers.join(d), ...params.rows.map((r) => r.join(d))].join('\n');
}

category('Generic Parser', () => {
  test('Generic: parses CSV with comma delimiter', async () => {
    const text = makeGenericCsv({
      headers: ['Protein', 'Intensity A', 'Intensity B', 'Intensity C'],
      rows: [
        ['P12345', '1000', '2000', '3000'],
        ['P67890', '4000', '5000', '6000'],
      ],
    });
    const delimiter = detectDelimiter(text);
    expect(delimiter, ',');
    const df = DG.DataFrame.fromCsv(text, {delimiter: ','});
    expect(df.columns.length, 4);
    expect(df.rowCount, 2);
  });

  test('Generic: parses TSV with tab delimiter', async () => {
    const text = makeGenericCsv({
      headers: ['Protein', 'Intensity A', 'Intensity B', 'Intensity C'],
      rows: [
        ['P12345', '1000', '2000', '3000'],
        ['P67890', '4000', '5000', '6000'],
      ],
      delimiter: '\t',
    });
    const delimiter = detectDelimiter(text);
    expect(delimiter, '\t');
    const df = DG.DataFrame.fromCsv(text, {delimiter: '\t'});
    expect(df.columns.length, 4);
    expect(df.rowCount, 2);
  });

  test('Generic: auto-suggests protein ID column', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein IDs', ['P12345']),
      DG.Column.fromStrings('Gene names', ['BRCA1']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Intensity A', [1000]),
    ]);
    const suggested = autoSuggestProteinIdColumn(df);
    expect(suggested !== null, true);
    expect(suggested!.name, 'Protein IDs');
  });

  test('Generic: auto-suggests intensity columns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P12345']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Intensity Sample1', [1000]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'LFQ Sample2', [2000]),
      DG.Column.fromStrings('notes', ['some note']),
    ]);
    const suggested = autoSuggestIntensityColumns(df);
    expect(suggested.length, 2);
    expect(suggested.includes('Intensity Sample1'), true);
    expect(suggested.includes('LFQ Sample2'), true);
  });

  test('Generic: detects raw intensities for log2 toggle', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'vals', [50000, 100000, 200000]),
    ]);
    const result = detectLog2Status(df, ['vals']);
    expect(result.isLog2, false);
  });

  test('Generic: detects already-log2 data', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'vals', [15.2, 18.7, 12.1]),
    ]);
    const result = detectLog2Status(df, ['vals']);
    expect(result.isLog2, true);
  });

  test('Generic: log2TransformColumns produces correct values', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Intensity A', [1024, 0, 2048]),
    ]);
    log2TransformColumns(df, ['Intensity A']);
    const log2Col = df.col('log2(Intensity A)');
    expect(log2Col !== null, true);
    if (log2Col) {
      expect(Math.abs(log2Col.get(0)! - 10.0) < 0.001, true);
      expect(log2Col.isNone(1), true); // zero input -> FLOAT_NULL
      expect(Math.abs(log2Col.get(2)! - 11.0) < 0.001, true);
      expect(log2Col.semType, SEMTYPE.INTENSITY);
    }
    const origCol = df.col('Intensity A');
    expect(origCol!.semType, SEMTYPE.INTENSITY);
  });

  test('Generic: copyAsLog2Columns copies values without transform', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Intensity A', [15.2, 18.7]),
    ]);
    copyAsLog2Columns(df, ['Intensity A']);
    const log2Col = df.col('log2(Intensity A)');
    expect(log2Col !== null, true);
    if (log2Col) {
      expect(Math.abs(log2Col.get(0)! - 15.2) < 0.001, true);
      expect(Math.abs(log2Col.get(1)! - 18.7) < 0.001, true);
    }
  });

  test('Generic: addPrimaryColumnIfNeeded creates column for semicolons', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein', ['P12345;Q67890', 'P11111']),
    ]);
    addPrimaryColumnIfNeeded(df, 'Protein', 'Primary Protein', SEMTYPE.PROTEIN_ID);
    const primaryCol = df.col('Primary Protein');
    expect(primaryCol !== null, true);
    if (primaryCol) {
      expect(primaryCol.get(0), 'P12345');
      expect(primaryCol.get(1), 'P11111');
    }
  });

  test('Generic: addPrimaryColumnIfNeeded skips when no semicolons', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein', ['P12345', 'P11111']),
    ]);
    addPrimaryColumnIfNeeded(df, 'Protein', 'Primary Protein', SEMTYPE.PROTEIN_ID);
    const primaryCol = df.col('Primary Protein');
    expect(primaryCol, null);
  });

  test('Generic: assigns semantic types to selected columns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein', ['P12345']),
    ]);
    const col = df.col('Protein')!;
    col.semType = SEMTYPE.PROTEIN_ID;
    expect(col.semType, SEMTYPE.PROTEIN_ID);
  });
});

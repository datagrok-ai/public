import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {df, look} from '../helpers';

// Source: core/client/d4/lib/src/viewers/line_chart/line_chart_look.dart
// `LineChartLook.auto(df)` — runs on viewer attach when xColumnName is not
// pinned via factory options. Four-stage X-axis pick:
//   1. existing df[xColumnName] resolves to a real column → keep it.
//   2. df.columns.numerical.length <= 10000:
//        firstOrDefault(df.columns.numerical, (c) =>
//          c is DateTimeColumn ? c.minDate != c.maxDate
//                              : xKeywords.any((k) =>
//                                  c.name.toLowerCase().startsWith(k) ||
//                                  c.name.toLowerCase().endsWith(k)))
//      where xKeywords = ['time', 'date', 'year', 'day', 'week',
//                         'timestamp', 'offset'].
//   3. firstWhere over ALL df.columns: DateTimeColumn with minDate != maxDate.
//   4. firstOrNull(df.columns.numerical).
//
// DateTimeColumn.isNumerical returns true (column.dart), so DateTime columns
// already participate in stage 2. Stage 3 only fires when stage 2 was skipped
// (numerical.length > 10000). We don't pin that path here.
//
// Iteration is in declaration order: the FIRST numerical column that
// satisfies the predicate wins. Tests deliberately put a non-matching numeric
// column first so a positive result proves the keyword/type predicate fired
// rather than the stage-4 first-numerical fallback.

// Each row: [columnName, dartType, values]
const alphaRow: [string, string, number[]] = ['alpha', 'double', [1, 2, 3, 4, 5]];

function expectAutoX(rows: Array<[string, string, any[]]>, expected: string): void {
  expect(df(rows).plot.line().props['xColumnName'], expected);
}

function addEventCol(d: DG.DataFrame, name: string = 'event_at'): void {
  d.columns.addNew(name, DG.COLUMN_TYPE.DATE_TIME);
  d.set(name, 0, '2020-01-01T00:00:00.000Z');
  d.set(name, 1, '2020-02-01T00:00:00.000Z');
  d.set(name, 2, '2020-03-01T00:00:00.000Z');
  d.set(name, 3, '2020-04-01T00:00:00.000Z');
  d.set(name, 4, '2020-05-01T00:00:00.000Z');
}

function addFlatDateCol(d: DG.DataFrame, name: string = 'flat_date'): void {
  d.columns.addNew(name, DG.COLUMN_TYPE.DATE_TIME);
  for (let i = 0; i < 5; i++) d.set(name, i, '2020-01-01T00:00:00.000Z');
}

category('AI: Viewers: LineChart auto-detect', () => {

  // === Per-keyword startsWith coverage (lowercase, prefix at start) ===
  // Each test prepends a non-matching numeric column so the stage-4 fallback
  // would yield 'alpha' if stage 2 did not match.

  // [keyword, columnName, dartType, values]
  const starts: Array<[string, string, string, any[]]> = [
    ['time', 'time_seconds', 'double', [10, 20, 30, 40, 50]],
    ['date', 'date_index', 'double', [10, 20, 30, 40, 50]],
    ['year', 'year_value', 'int', [2020, 2021, 2022, 2023, 2024]],
    ['day', 'day_index', 'int', [1, 2, 3, 4, 5]],
    ['week', 'week_number', 'int', [1, 2, 3, 4, 5]],
    ['timestamp', 'timestamp_utc', 'double', [1e9, 1.1e9, 1.2e9, 1.3e9, 1.4e9]],
    ['offset', 'offset_ms', 'double', [10, 20, 30, 40, 50]],
  ];
  for (const s of starts) {
    const [kw, name, type, values] = s;
    test(`keyword "${kw}" at start picks numeric column as X`, async () =>
      expectAutoX([alphaRow, [name, type, values]], name));
  }

  // === Per-keyword endsWith coverage (suffix at end) ===

  const ends: Array<[string, string, string, any[]]> = [
    ['time', 'elapsed_time', 'double', [10, 20, 30, 40, 50]],
    ['date', 'birthdate', 'double', [10, 20, 30, 40, 50]],
    ['year', 'fiscal_year', 'int', [2020, 2021, 2022, 2023, 2024]],
    ['day', 'birth_day', 'int', [1, 2, 3, 4, 5]],
    ['week', 'this_week', 'int', [1, 2, 3, 4, 5]],
    ['timestamp', 'creation_timestamp', 'double', [1e9, 1.1e9, 1.2e9, 1.3e9, 1.4e9]],
    ['offset', 'utc_offset', 'double', [10, 20, 30, 40, 50]],
  ];
  for (const e of ends) {
    const [kw, name, type, values] = e;
    test(`keyword "${kw}" at end picks numeric column as X`, async () =>
      expectAutoX([alphaRow, [name, type, values]], name));
  }

  // === Exact-match (full name equals keyword) ===

  test('column named exactly "time" picked as X', async () =>
    expectAutoX([alphaRow, ['time', 'double', [10, 20, 30, 40, 50]]], 'time'));

  test('column named exactly "year" picked as X', async () =>
    expectAutoX([alphaRow, ['year', 'int', [2020, 2021, 2022, 2023, 2024]]], 'year'));

  // === Case insensitivity ===
  // The predicate lower-cases the column name before comparing.

  test('case-insensitive: capitalized "Time" matches keyword', async () =>
    expectAutoX([alphaRow, ['Time', 'double', [10, 20, 30, 40, 50]]], 'Time'));

  test('case-insensitive: uppercase "TIMESTAMP_UTC" matches keyword', async () =>
    expectAutoX([alphaRow, ['TIMESTAMP_UTC', 'double', [1e9, 1.1e9, 1.2e9, 1.3e9, 1.4e9]]], 'TIMESTAMP_UTC'));

  test('case-insensitive: mixed-case "BirthDate" matches keyword', async () =>
    expectAutoX([alphaRow, ['BirthDate', 'double', [10, 20, 30, 40, 50]]], 'BirthDate'));

  // === Negative: no keyword anywhere → stage-4 fallback to first numerical ===

  test('no keyword/DateTime → first numerical column wins (stage-4 fallback)', async () =>
    expectAutoX([
      ['name', 'string', ['a', 'b', 'c', 'd', 'e']],
      ['salary', 'double', [100, 200, 300, 400, 500]],
      ['bonus', 'double', [10, 20, 30, 40, 50]],
    ], 'salary'));

  test('keyword as middle substring does NOT match (must be prefix or suffix)', async () =>
    expectAutoX([alphaRow, ['cdate_id', 'double', [10, 20, 30, 40, 50]]], 'alpha'));

  test('non-keyword name with numeric-only chars does not falsely match', async () =>
    expectAutoX([
      ['value', 'double', [1, 2, 3, 4, 5]],
      ['amount', 'double', [10, 20, 30, 40, 50]],
    ], 'value'));

  // === Column-type autodetection: DateTimeColumn ===

  test('DateTimeColumn with min != max wins over plain numeric', async () => {
    const d = df([alphaRow]);
    addEventCol(d);
    expect(d.plot.line().props['xColumnName'], 'event_at');
  });

  test('DateTimeColumn with min == max is rejected; falls through to keyword', async () => {
    const d = df([alphaRow]);
    addFlatDateCol(d);
    d.columns.add(DG.Column.fromList('double', 'time_ms', [10, 20, 30, 40, 50]));
    expect(d.plot.line().props['xColumnName'], 'time_ms');
  });

  test('DateTimeColumn with min == max and no keyword → stage-4 first numerical', async () => {
    const d = df([alphaRow, ['beta', 'double', [10, 20, 30, 40, 50]]]);
    addFlatDateCol(d);
    expect(d.plot.line().props['xColumnName'], 'alpha');
  });

  // === Iteration order between DateTime and keyword match ===

  test('DateTime column placed first in DF wins over later keyword column', async () => {
    const d = df([alphaRow]);
    addEventCol(d);
    d.columns.add(DG.Column.fromList('double', 'time_ms', [10, 20, 30, 40, 50]));
    expect(d.plot.line().props['xColumnName'], 'event_at');
  });

  test('keyword column placed first in DF wins over later DateTime column', async () => {
    const d = df([alphaRow, ['time_ms', 'double', [10, 20, 30, 40, 50]]]);
    addEventCol(d);
    expect(d.plot.line().props['xColumnName'], 'time_ms');
  });

  // === Existing-name pin overrides keyword/type detection ===

  test('explicit xColumnName option overrides auto-pick', async () => {
    const d = df([alphaRow, ['time_ms', 'double', [10, 20, 30, 40, 50]]]);
    expect(d.plot.line({xColumnName: 'alpha'}).props['xColumnName'], 'alpha');
  });

  // === String column with keyword name is NOT picked (predicate runs only over numerical) ===

  test('string column with keyword-matching name is ignored (predicate is over numerical)', async () =>
    expectAutoX([
      ['time', 'string', ['t1', 't2', 't3', 't4', 't5']],
      ['alpha', 'double', [1, 2, 3, 4, 5]],
      ['beta', 'double', [10, 20, 30, 40, 50]],
    ], 'alpha'));

}, {owner: 'agolovko@datagrok.ai'});

// Trellis path: trellis_plot_core.dart:729 wires the inner-viewer look as
//   `look.innerViewerLook..auto(dataFrame, trellisLook: look)..apply()`.
// `LineChartLook.auto(DataFrame df, {TrellisPlotLook trellisLook})` declares
// `trellisLook` but never reads it inside the body — the same four-stage
// X-axis pick fires whether the line chart is standalone or hosted inside a
// trellis. We mirror a representative subset against the resolved inner look
// surfaced via `getOptions(true).look['innerViewerLook']` (same readback
// pattern as ApiTests/src/ai/reported-issues/grok-19466.ts).
//
// Each dataframe carries a `split` string column passed as the trellis
// outer-axis splitter (`xColumnNames: ['split']`). String columns are
// invisible to LineChartLook.auto stage 2 (`df.columns.numerical`) and to
// stage 3 (DateTimeColumn predicate) and stage 4 (`columns.numerical`), so
// adding `split` cannot perturb the X auto-pick — only the trellis cell grid.

const trellisAlpha: [string, string, number[]] = ['alpha', 'double', [1, 2, 3, 4, 5, 6]];
const trellisSplit: [string, string, string[]] = ['split', 'string', ['a', 'b', 'a', 'b', 'a', 'b']];

function trellisInnerX(rows: Array<[string, string, any[]]>): any {
  const tv = DG.Viewer.trellisPlot(df(rows), {viewerType: DG.VIEWER.LINE_CHART, xColumnNames: ['split']});
  return look(tv)['innerViewerLook']['xColumnName'];
}

category('AI: Viewers: LineChart auto-detect (as trellis inner viewer)', () => {

  test('keyword "time" at start picks numeric column as X (trellis inner)', async () =>
    expect(trellisInnerX(
      [trellisAlpha, ['time_seconds', 'double', [10, 20, 30, 40, 50, 60]], trellisSplit]), 'time_seconds'));

  test('keyword "date" at end picks numeric column as X (trellis inner)', async () =>
    expect(trellisInnerX(
      [trellisAlpha, ['birthdate', 'double', [10, 20, 30, 40, 50, 60]], trellisSplit]), 'birthdate'));

  test('case-insensitive uppercase keyword match (trellis inner)', async () =>
    expect(trellisInnerX(
      [trellisAlpha, ['TIMESTAMP_UTC', 'double', [1e9, 1.1e9, 1.2e9, 1.3e9, 1.4e9, 1.5e9]], trellisSplit]),
    'TIMESTAMP_UTC'));

  test('no keyword/DateTime → first numerical wins (stage-4 fallback, trellis inner)', async () =>
    expect(trellisInnerX([
      ['name', 'string', ['a', 'b', 'c', 'a', 'b', 'c']],
      ['salary', 'double', [100, 200, 300, 400, 500, 600]],
      ['bonus', 'double', [10, 20, 30, 40, 50, 60]],
      trellisSplit,
    ]), 'salary'));

  test('DateTimeColumn with min != max wins over plain numeric (trellis inner)', async () => {
    const d = df([['alpha', 'double', [1, 2, 3, 4, 5]]]);
    addEventCol(d);
    d.columns.add(DG.Column.fromList('string', 'split', ['a', 'b', 'a', 'b', 'a']));
    const tv = DG.Viewer.trellisPlot(d, {viewerType: DG.VIEWER.LINE_CHART, xColumnNames: ['split']});
    expect(look(tv)['innerViewerLook']['xColumnName'], 'event_at');
  });

}, {owner: 'agolovko@datagrok.ai'});

import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

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

category('AI: Viewers: LineChart auto-detect', () => {

  // === Per-keyword startsWith coverage (lowercase, prefix at start) ===
  // Each test prepends a non-matching numeric column so the stage-4 fallback
  // would yield 'alpha' if stage 2 did not match. A pass on the keyword
  // column proves the keyword predicate fired.

  test('keyword "time" at start picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'time_seconds', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'time_seconds');
  });

  test('keyword "date" at start picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'date_index', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'date_index');
  });

  test('keyword "year" at start picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('int', 'year_value', [2020, 2021, 2022, 2023, 2024]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'year_value');
  });

  test('keyword "day" at start picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('int', 'day_index', [1, 2, 3, 4, 5]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'day_index');
  });

  test('keyword "week" at start picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('int', 'week_number', [1, 2, 3, 4, 5]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'week_number');
  });

  test('keyword "timestamp" at start picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'timestamp_utc', [1e9, 1.1e9, 1.2e9, 1.3e9, 1.4e9]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'timestamp_utc');
  });

  test('keyword "offset" at start picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'offset_ms', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'offset_ms');
  });

  // === Per-keyword endsWith coverage (suffix at end) ===

  test('keyword "time" at end picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'elapsed_time', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'elapsed_time');
  });

  test('keyword "date" at end picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'birthdate', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'birthdate');
  });

  test('keyword "year" at end picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('int', 'fiscal_year', [2020, 2021, 2022, 2023, 2024]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'fiscal_year');
  });

  test('keyword "day" at end picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('int', 'birth_day', [1, 2, 3, 4, 5]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'birth_day');
  });

  test('keyword "week" at end picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('int', 'this_week', [1, 2, 3, 4, 5]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'this_week');
  });

  test('keyword "timestamp" at end picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'creation_timestamp', [1e9, 1.1e9, 1.2e9, 1.3e9, 1.4e9]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'creation_timestamp');
  });

  test('keyword "offset" at end picks numeric column as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'utc_offset', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'utc_offset');
  });

  // === Exact-match (full name equals keyword — startsWith and endsWith both true) ===

  test('column named exactly "time" picked as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'time', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'time');
  });

  test('column named exactly "year" picked as X', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('int', 'year', [2020, 2021, 2022, 2023, 2024]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'year');
  });

  // === Case insensitivity ===
  // The predicate lower-cases the column name before comparing. Mixed-case
  // and all-uppercase names must match.

  test('case-insensitive: capitalized "Time" matches keyword', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'Time', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'Time');
  });

  test('case-insensitive: uppercase "TIMESTAMP_UTC" matches keyword', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'TIMESTAMP_UTC', [1e9, 1.1e9, 1.2e9, 1.3e9, 1.4e9]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'TIMESTAMP_UTC');
  });

  test('case-insensitive: mixed-case "BirthDate" matches keyword', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'BirthDate', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'BirthDate');
  });

  // === Negative: no keyword anywhere → stage-4 fallback to first numerical ===

  test('no keyword/DateTime → first numerical column wins (stage-4 fallback)', async () => {
    // None of these names start or end with any xKeyword, so stage 2 yields
    // null. Stage 3 finds no DateTimeColumn. Stage 4 returns
    // firstOrNull(df.columns.numerical) → 'salary' (declaration order).
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'name', ['a', 'b', 'c', 'd', 'e']),
      DG.Column.fromList('double', 'salary', [100, 200, 300, 400, 500]),
      DG.Column.fromList('double', 'bonus', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'salary');
  });

  test('keyword as middle substring does NOT match (must be prefix or suffix)', async () => {
    // 'cdate_id' contains 'date' but neither starts nor ends with it. The
    // predicate is startsWith||endsWith — middle substrings are rejected.
    // With no other keyword/DateTime match, stage-4 fallback picks the first
    // numerical column, which is 'cdate_id' itself (only numeric column).
    // To prove the negative we need TWO numerical columns where the
    // would-be-match comes second — then stage 4 returns the FIRST
    // ('alpha'), confirming the keyword predicate did not fire.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'cdate_id', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'alpha');
  });

  test('non-keyword name with numeric-only chars does not falsely match', async () => {
    // 'value' / 'amount' contain none of the xKeywords as prefix/suffix.
    // Stage 4 returns the first numerical column.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'value', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'amount', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'value');
  });

  // === Column-type autodetection: DateTimeColumn ===

  test('DateTimeColumn with min != max wins over plain numeric', async () => {
    // alpha is iterated first (no keyword, not DateTime) → skip. Then the
    // DateTime column with distinct dates satisfies the predicate.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
    ]);
    df.columns.addNew('event_at', DG.COLUMN_TYPE.DATE_TIME);
    df.set('event_at', 0, '2020-01-01T00:00:00.000Z');
    df.set('event_at', 1, '2020-02-01T00:00:00.000Z');
    df.set('event_at', 2, '2020-03-01T00:00:00.000Z');
    df.set('event_at', 3, '2020-04-01T00:00:00.000Z');
    df.set('event_at', 4, '2020-05-01T00:00:00.000Z');
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'event_at');
  });

  test('DateTimeColumn with min == max is rejected; falls through to keyword', async () => {
    // The DateTime predicate requires minDate != maxDate. All-equal dates
    // make the column fail the predicate, so iteration continues. The next
    // numeric matches 'time' → it wins.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
    ]);
    df.columns.addNew('flat_date', DG.COLUMN_TYPE.DATE_TIME);
    for (let i = 0; i < 5; i++)
      df.set('flat_date', i, '2020-01-01T00:00:00.000Z');
    df.columns.add(DG.Column.fromList('double', 'time_ms', [10, 20, 30, 40, 50]));
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'time_ms');
  });

  test('DateTimeColumn with min == max and no keyword → stage-4 first numerical', async () => {
    // Flat DateTime fails stage 2's DateTime predicate AND stage 3
    // (firstWhere over all columns also requires minDate != maxDate).
    // Falls all the way to stage 4: first numerical column = 'alpha'.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50]),
    ]);
    df.columns.addNew('flat_date', DG.COLUMN_TYPE.DATE_TIME);
    for (let i = 0; i < 5; i++)
      df.set('flat_date', i, '2020-01-01T00:00:00.000Z');
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'alpha');
  });

  // === Iteration order between DateTime and keyword match ===

  test('DateTime column placed first in DF wins over later keyword column', async () => {
    // firstOrDefault returns the first match. DateTime is iterated before
    // 'time_ms' so DateTime wins.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
    ]);
    df.columns.addNew('event_at', DG.COLUMN_TYPE.DATE_TIME);
    df.set('event_at', 0, '2020-01-01T00:00:00.000Z');
    df.set('event_at', 1, '2020-02-01T00:00:00.000Z');
    df.set('event_at', 2, '2020-03-01T00:00:00.000Z');
    df.set('event_at', 3, '2020-04-01T00:00:00.000Z');
    df.set('event_at', 4, '2020-05-01T00:00:00.000Z');
    df.columns.add(DG.Column.fromList('double', 'time_ms', [10, 20, 30, 40, 50]));
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'event_at');
  });

  test('keyword column placed first in DF wins over later DateTime column', async () => {
    // Reverse order: 'time_ms' iterated first, matches keyword → wins.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'time_ms', [10, 20, 30, 40, 50]),
    ]);
    df.columns.addNew('event_at', DG.COLUMN_TYPE.DATE_TIME);
    df.set('event_at', 0, '2020-01-01T00:00:00.000Z');
    df.set('event_at', 1, '2020-02-01T00:00:00.000Z');
    df.set('event_at', 2, '2020-03-01T00:00:00.000Z');
    df.set('event_at', 3, '2020-04-01T00:00:00.000Z');
    df.set('event_at', 4, '2020-05-01T00:00:00.000Z');
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'time_ms');
  });

  // === Existing-name pin overrides keyword/type detection ===

  test('explicit xColumnName option overrides auto-pick', async () => {
    // df[xColumnName] resolves first. Even though 'time_ms' would otherwise
    // be picked, pinning 'alpha' wins.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'time_ms', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line({xColumnName: 'alpha'});
    expect(v.props['xColumnName'], 'alpha');
  });

  // === String column with keyword name is NOT picked (predicate runs only over numerical) ===

  test('string column with keyword-matching name is ignored (predicate is over numerical)', async () => {
    // Stage 2 iterates df.columns.numerical, so a string 'time' column is
    // invisible to the predicate. Stage-4 fallback picks the first numerical
    // column = 'alpha'.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'time', ['t1', 't2', 't3', 't4', 't5']),
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50]),
    ]);
    const v = df.plot.line();
    expect(v.props['xColumnName'], 'alpha');
  });

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

function trellisLineInnerLook(df: DG.DataFrame): {[key: string]: any} {
  const v = DG.Viewer.trellisPlot(df, {
    viewerType: DG.VIEWER.LINE_CHART,
    xColumnNames: ['split'],
  });
  return v.getOptions(true).look['innerViewerLook'];
}

category('AI: Viewers: LineChart auto-detect (as trellis inner viewer)', () => {

  test('keyword "time" at start picks numeric column as X (trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6]),
      DG.Column.fromList('double', 'time_seconds', [10, 20, 30, 40, 50, 60]),
      DG.Column.fromList('string', 'split', ['a', 'b', 'a', 'b', 'a', 'b']),
    ]);
    expect(trellisLineInnerLook(df)['xColumnName'], 'time_seconds');
  });

  test('keyword "date" at end picks numeric column as X (trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6]),
      DG.Column.fromList('double', 'birthdate', [10, 20, 30, 40, 50, 60]),
      DG.Column.fromList('string', 'split', ['a', 'b', 'a', 'b', 'a', 'b']),
    ]);
    expect(trellisLineInnerLook(df)['xColumnName'], 'birthdate');
  });

  test('case-insensitive uppercase keyword match (trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6]),
      DG.Column.fromList('double', 'TIMESTAMP_UTC', [1e9, 1.1e9, 1.2e9, 1.3e9, 1.4e9, 1.5e9]),
      DG.Column.fromList('string', 'split', ['a', 'b', 'a', 'b', 'a', 'b']),
    ]);
    expect(trellisLineInnerLook(df)['xColumnName'], 'TIMESTAMP_UTC');
  });

  test('no keyword/DateTime → first numerical wins (stage-4 fallback, trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'name', ['a', 'b', 'c', 'a', 'b', 'c']),
      DG.Column.fromList('double', 'salary', [100, 200, 300, 400, 500, 600]),
      DG.Column.fromList('double', 'bonus', [10, 20, 30, 40, 50, 60]),
      DG.Column.fromList('string', 'split', ['a', 'b', 'a', 'b', 'a', 'b']),
    ]);
    expect(trellisLineInnerLook(df)['xColumnName'], 'salary');
  });

  test('DateTimeColumn with min != max wins over plain numeric (trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5]),
    ]);
    df.columns.addNew('event_at', DG.COLUMN_TYPE.DATE_TIME);
    df.set('event_at', 0, '2020-01-01T00:00:00.000Z');
    df.set('event_at', 1, '2020-02-01T00:00:00.000Z');
    df.set('event_at', 2, '2020-03-01T00:00:00.000Z');
    df.set('event_at', 3, '2020-04-01T00:00:00.000Z');
    df.set('event_at', 4, '2020-05-01T00:00:00.000Z');
    df.columns.add(DG.Column.fromList('string', 'split', ['a', 'b', 'a', 'b', 'a']));
    expect(trellisLineInnerLook(df)['xColumnName'], 'event_at');
  });

}, {owner: 'agolovko@datagrok.ai'});

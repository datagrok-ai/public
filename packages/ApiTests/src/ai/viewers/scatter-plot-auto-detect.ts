import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Source: core/client/d4/lib/src/viewers/scatterplot/scatterplot_look.dart
// `ScatterPlotLook.auto(df)` — runs on viewer attach when X/Y are not pinned by
// the caller. Two-stage pick:
//   1. If both xColumnName and yColumnName are null AND df.columns.length <
//      10000 → `findMatchingErrorBarColumns(df)` (utils.dart:1068) is consulted
//      first. When it returns a tagged or suffix-matched whisker triple
//      (<base>, <range|null>, <min|null>, <max|null>), the base column wins
//      that axis and the whisker* slots are filled.
//   2. If after step 1 either axis is still null, the auto-pick falls back to
//      the diversity-sorted numeric+categorical concat (tilde-prefixed columns
//      are filtered out via `filterOutTildeColumns`).
// Everything below pins exactly that contract via the public JS surface
// (props / getOptions(true).look). No canvas-geometry getters are touched —
// auto() is settled by the time df.plot.scatter() returns.
//
// Suffix dictionary (lowercase, see utils.dart `completeSuffixes`):
//   range  : ' error range', ' whisker range', ' range', ' std range', ' stddev range'
//   min    : ' error min',   ' whisker min',   ' min',   ' std min',   ' stddev min'
//   max    : ' error max',   ' whisker max',   ' max',   ' std max',   ' stddev max'
//   range* : ' std', ' stddev'  (bare std / stddev are folded into range)
// firstWhere() iterates in declaration order, so the bare ' min' / ' max'
// suffix shadows ' std min' / ' stddev min' / ' std max' / ' stddev max' —
// "a std min" matches ' min' first and yields base "a std" (not "a"), which
// the Phase-3 byName() then fails to resolve unless a literal "a std" column
// also exists. We pin only the suffix families that actually fire (' whisker'
// and ' error' come BEFORE the bare suffix in the variant list).

category('AI: Viewers: ScatterPlot auto-detect', () => {

  // === X/Y auto-pick (no whiskers) ===

  test('two numerical columns pick distinct X and Y', async () => {
    // Auto-pick sorts numerical columns by ColumnHistogram bin diversity desc.
    // With two numerical columns of clearly different diversities, X and Y
    // must end up on the two distinct columns.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50, 60, 70, 80]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const x = v.props['xColumnName'];
    const y = v.props['yColumnName'];
    expect(x === 'alpha' || x === 'beta', true);
    expect(y === 'alpha' || y === 'beta', true);
    expect(x !== y, true);
  });

  test('single numerical column places same column on both axes', async () => {
    // auto() does `yColumnName ??= (columns.length > 1 ? columns[1].name :
    // xColumnName)` — exactly one numeric candidate must collapse Y onto X.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'only', [1, 2, 3, 4, 5, 6, 7, 8]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    expect(v.props['xColumnName'], 'only');
    expect(v.props['yColumnName'], 'only');
  });

  test('numerical preferred over categorical when both present', async () => {
    // The candidate iterable concats numerical (sorted) FIRST, then
    // categorical. With one numerical + one categorical, X must be the
    // numerical and Y falls onto the categorical (length > 1 path).
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'cat', ['a', 'b', 'c', 'a', 'b', 'c', 'a', 'b']),
      DG.Column.fromList('double', 'num', [1, 2, 3, 4, 5, 6, 7, 8]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    expect(v.props['xColumnName'], 'num');
    expect(v.props['yColumnName'], 'cat');
  });

  test('categorical-only dataframe falls back to categoricals', async () => {
    // No numerical columns at all — auto() must still produce a valid X/Y
    // from the categorical fallback.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'cat1', ['a', 'b', 'c', 'a', 'b', 'c']),
      DG.Column.fromList('string', 'cat2', ['x', 'y', 'z', 'x', 'y', 'z']),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const x = v.props['xColumnName'];
    const y = v.props['yColumnName'];
    expect(x === 'cat1' || x === 'cat2', true);
    expect(y === 'cat1' || y === 'cat2', true);
    expect(x !== y, true);
  });

  test('tilde-prefixed columns excluded from auto-pick when others exist', async () => {
    // filterOutTildeColumns: if any non-tilde exists, all '~' columns are
    // dropped from the candidate list. With two tildes + two regular numerics,
    // X and Y must both land on the regulars.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', '~hidden1', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', '~hidden2', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'visible1', [100, 200, 300, 400, 500, 600, 700, 800]),
      DG.Column.fromList('double', 'visible2', [9, 18, 27, 36, 45, 54, 63, 72]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const x = v.props['xColumnName'];
    const y = v.props['yColumnName'];
    expect(x === 'visible1' || x === 'visible2', true);
    expect(y === 'visible1' || y === 'visible2', true);
    expect(x !== y, true);
    expect((x as string).startsWith('~'), false);
    expect((y as string).startsWith('~'), false);
  });

  // === Whisker auto-detect: suffix family coverage ===

  test('suffix " whisker min" / " whisker max" auto-detected', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a whisker min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      DG.Column.fromList('double', 'a whisker max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b whisker min', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'b whisker max', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const xBase = v.props['xColumnName'];
    expect(xBase === 'a' || xBase === 'b', true);
    expect(v.props['xWhiskerMinColumnName'], xBase + ' whisker min');
    expect(v.props['xWhiskerMaxColumnName'], xBase + ' whisker max');
    const yBase = v.props['yColumnName'];
    expect(yBase !== xBase, true);
    expect(v.props['yWhiskerMinColumnName'], yBase + ' whisker min');
    expect(v.props['yWhiskerMaxColumnName'], yBase + ' whisker max');
  });

  test('suffix " error min" / " error max" auto-detected', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a error min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      DG.Column.fromList('double', 'a error max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b error min', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'b error max', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const xBase = v.props['xColumnName'];
    const yBase = v.props['yColumnName'];
    expect(xBase === 'a' || xBase === 'b', true);
    expect(yBase !== xBase, true);
    expect(v.props['xWhiskerMinColumnName'], xBase + ' error min');
    expect(v.props['xWhiskerMaxColumnName'], xBase + ' error max');
    expect(v.props['yWhiskerMinColumnName'], yBase + ' error min');
    expect(v.props['yWhiskerMaxColumnName'], yBase + ' error max');
  });

  test('suffix " whisker range" auto-detected (no min/max)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a whisker range', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b whisker range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const xBase = v.props['xColumnName'];
    const yBase = v.props['yColumnName'];
    expect(xBase === 'a' || xBase === 'b', true);
    expect(yBase !== xBase, true);
    expect(v.props['xWhiskerRangeColumnName'], xBase + ' whisker range');
    expect(v.props['yWhiskerRangeColumnName'], yBase + ' whisker range');
    expect(v.props['xWhiskerMinColumnName'] == null || v.props['xWhiskerMinColumnName'] === '', true);
    expect(v.props['yWhiskerMinColumnName'] == null || v.props['yWhiskerMinColumnName'] === '', true);
  });

  test('suffix " error range" auto-detected (no min/max)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a error range', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b error range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const xBase = v.props['xColumnName'];
    const yBase = v.props['yColumnName'];
    expect(xBase === 'a' || xBase === 'b', true);
    expect(yBase !== xBase, true);
    expect(v.props['xWhiskerRangeColumnName'], xBase + ' error range');
    expect(v.props['yWhiskerRangeColumnName'], yBase + ' error range');
  });

  test('bare " std" suffix folded into range', async () => {
    // Per utils.dart: bare ' std' / ' stddev' suffixes are deposited into
    // rangeColumns (not min/max). So a column "a std" alone resolves the
    // axis 'a' with whisker.range = "a std".
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a std', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b std', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const xBase = v.props['xColumnName'];
    const yBase = v.props['yColumnName'];
    expect(xBase === 'a' || xBase === 'b', true);
    expect(yBase !== xBase, true);
    expect(v.props['xWhiskerRangeColumnName'], xBase + ' std');
    expect(v.props['yWhiskerRangeColumnName'], yBase + ' std');
    expect(v.props['xWhiskerMinColumnName'] == null || v.props['xWhiskerMinColumnName'] === '', true);
    expect(v.props['xWhiskerMaxColumnName'] == null || v.props['xWhiskerMaxColumnName'] === '', true);
  });

  test('bare " stddev" suffix folded into range', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a stddev', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b stddev', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const xBase = v.props['xColumnName'];
    const yBase = v.props['yColumnName'];
    expect(xBase === 'a' || xBase === 'b', true);
    expect(yBase !== xBase, true);
    expect(v.props['xWhiskerRangeColumnName'], xBase + ' stddev');
    expect(v.props['yWhiskerRangeColumnName'], yBase + ' stddev');
  });

  // === Tag-based detection ===

  test('tag "whisker.range" on numeric column drives axis pick + range slot', async () => {
    // findMatchingErrorBarColumns Phase 1 scans column tags; isMatchingTag
    // looks for 'range' (or 'whisker range') in the tag key. 'whisker.range'
    // contains both, so the tagged column wins the axis pick and the tag
    // VALUE becomes the range column name.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'rng1', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      DG.Column.fromList('double', 'rng2', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]);
    df.col('alpha')!.tags['whisker.range'] = 'rng1';
    df.col('beta')!.tags['whisker.range'] = 'rng2';
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    const xCol = v.props['xColumnName'];
    const yCol = v.props['yColumnName'];
    expect(xCol === 'alpha' || xCol === 'beta', true);
    expect(yCol === 'alpha' || yCol === 'beta', true);
    expect(xCol !== yCol, true);
    const expRange = (axis: string) => axis === 'alpha' ? 'rng1' : 'rng2';
    expect(v.props['xWhiskerRangeColumnName'], expRange(xCol as string));
    expect(v.props['yWhiskerRangeColumnName'], expRange(yCol as string));
  });

  // === Per-axis re-detect after options apply ===

  // grok_api.dart:411 wires the JS factory as
  //   `new ScatterPlotCore()..dataFrame = t..look.setOptionsJson(json)`
  // — dataFrame attaches FIRST (bulk auto runs with X/Y null and sets
  // _whiskersAutoDetected = true), THEN options like {x: 'a'} apply and
  // fire onLookChanged($xColumnName), which re-runs the per-axis whisker
  // detector (whiskersFeature.resetWhiskers, scatterplot_core.dart:370–372).
  // So passing {x: <base>} does NOT suppress whiskers — the matching min/max
  // pair on the chosen base must still light up.
  test('per-axis whisker re-detect honors xColumnName passed via options', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'alpha min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      DG.Column.fromList('double', 'alpha max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'beta min', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'beta max', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
    const v = df.plot.scatter({x: 'alpha'}) as DG.ScatterPlotViewer;
    expect(v.props['xColumnName'], 'alpha');
    expect(v.props['xWhiskerMinColumnName'], 'alpha min');
    expect(v.props['xWhiskerMaxColumnName'], 'alpha max');
  });

  test('per-axis whisker re-detect honors yColumnName passed via options', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'alpha min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      DG.Column.fromList('double', 'alpha max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'beta min', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'beta max', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
    const v = df.plot.scatter({y: 'beta'}) as DG.ScatterPlotViewer;
    expect(v.props['yColumnName'], 'beta');
    expect(v.props['yWhiskerMinColumnName'], 'beta min');
    expect(v.props['yWhiskerMaxColumnName'], 'beta max');
  });

}, {owner: 'agolovko@datagrok.ai'});

// Trellis path: trellis_plot_core.dart:729 wires the inner-viewer look as
//   `look.innerViewerLook..auto(dataFrame, trellisLook: look)..apply()`.
// `ScatterPlotLook.auto(DataFrame df, {TrellisPlotLook trellisLook})` declares
// `trellisLook` but never reads it inside the body — Phase 1
// (findMatchingErrorBarColumns: tag scan + suffix scan) and Phase 2
// (numerical-by-bin-diversity sort + categorical fallback) both fire
// identically whether the scatter plot is standalone or hosted inside a
// trellis. We mirror a representative subset against the resolved inner look
// surfaced via `getOptions(true).look['innerViewerLook']` (same readback
// pattern as ApiTests/src/ai/reported-issues/grok-19466.ts).
//
// Each dataframe carries a `split` string column passed as the trellis
// outer-axis splitter (`xColumnNames: ['split']`). It enters Phase 2's
// candidate iterable AFTER all numerical columns (concat is numerical first,
// categorical second) — declaring it last in the dataframe keeps Y on the
// second-most-diverse numerical for the two-numerical case. Phase 1 is
// unaffected: `findMatchingErrorBarColumns` skips non-numeric base columns.

function trellisScatterInnerLook(df: DG.DataFrame): {[key: string]: any} {
  const v = DG.Viewer.trellisPlot(df, {
    viewerType: DG.VIEWER.SCATTER_PLOT,
    xColumnNames: ['split'],
  });
  return v.getOptions(true).look['innerViewerLook'];
}

category('AI: Viewers: ScatterPlot auto-detect (as trellis inner viewer)', () => {

  test('two numerical columns pick distinct X and Y (trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('string', 'split', ['a', 'b', 'a', 'b', 'a', 'b', 'a', 'b']),
    ]);
    const inner = trellisScatterInnerLook(df);
    const x = inner['xColumnName'];
    const y = inner['yColumnName'];
    expect(x === 'alpha' || x === 'beta', true);
    expect(y === 'alpha' || y === 'beta', true);
    expect(x !== y, true);
  });

  test('tilde-prefixed columns excluded from auto-pick (trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', '~hidden1', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', '~hidden2', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'visible1', [100, 200, 300, 400, 500, 600, 700, 800]),
      DG.Column.fromList('double', 'visible2', [9, 18, 27, 36, 45, 54, 63, 72]),
      DG.Column.fromList('string', 'split', ['a', 'b', 'a', 'b', 'a', 'b', 'a', 'b']),
    ]);
    const inner = trellisScatterInnerLook(df);
    const x = inner['xColumnName'];
    const y = inner['yColumnName'];
    expect(x === 'visible1' || x === 'visible2', true);
    expect(y === 'visible1' || y === 'visible2', true);
    expect(x !== y, true);
    expect((x as string).startsWith('~'), false);
    expect((y as string).startsWith('~'), false);
  });

  test('suffix " whisker min" / " whisker max" auto-detected (trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a whisker min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      DG.Column.fromList('double', 'a whisker max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b whisker min', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'b whisker max', [11, 22, 32, 42, 52, 62, 72, 82]),
      DG.Column.fromList('string', 'split', ['s', 't', 's', 't', 's', 't', 's', 't']),
    ]);
    const inner = trellisScatterInnerLook(df);
    const xBase = inner['xColumnName'];
    const yBase = inner['yColumnName'];
    expect(xBase === 'a' || xBase === 'b', true);
    expect(yBase !== xBase, true);
    expect(inner['xWhiskerMinColumnName'], xBase + ' whisker min');
    expect(inner['xWhiskerMaxColumnName'], xBase + ' whisker max');
    expect(inner['yWhiskerMinColumnName'], yBase + ' whisker min');
    expect(inner['yWhiskerMaxColumnName'], yBase + ' whisker max');
  });

  test('suffix " error range" auto-detected (no min/max, trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'a', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'a error range', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      DG.Column.fromList('double', 'b', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'b error range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
      DG.Column.fromList('string', 'split', ['s', 't', 's', 't', 's', 't', 's', 't']),
    ]);
    const inner = trellisScatterInnerLook(df);
    const xBase = inner['xColumnName'];
    const yBase = inner['yColumnName'];
    expect(xBase === 'a' || xBase === 'b', true);
    expect(yBase !== xBase, true);
    expect(inner['xWhiskerRangeColumnName'], xBase + ' error range');
    expect(inner['yWhiskerRangeColumnName'], yBase + ' error range');
  });

  test('tag "whisker.range" drives axis pick + range slot (trellis inner)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'rng1', [0.1, 0.2, 0.15, 0.25, 0.18, 0.22, 0.19, 0.21]),
      DG.Column.fromList('double', 'rng2', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
      DG.Column.fromList('string', 'split', ['s', 't', 's', 't', 's', 't', 's', 't']),
    ]);
    df.col('alpha')!.tags['whisker.range'] = 'rng1';
    df.col('beta')!.tags['whisker.range'] = 'rng2';
    const inner = trellisScatterInnerLook(df);
    const xCol = inner['xColumnName'];
    const yCol = inner['yColumnName'];
    expect(xCol === 'alpha' || xCol === 'beta', true);
    expect(yCol === 'alpha' || yCol === 'beta', true);
    expect(xCol !== yCol, true);
    const expRange = (axis: string) => axis === 'alpha' ? 'rng1' : 'rng2';
    expect(inner['xWhiskerRangeColumnName'], expRange(xCol as string));
    expect(inner['yWhiskerRangeColumnName'], expRange(yCol as string));
  });

}, {owner: 'agolovko@datagrok.ai'});

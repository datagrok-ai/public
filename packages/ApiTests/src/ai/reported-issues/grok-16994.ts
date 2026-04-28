import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-16994: Scatter plot whisker columns
// (StdDev / errorbars) — public IScatterPlotSettings now exposes
// xWhisker{Min,Max,Range}ColumnName / yWhisker{Min,Max,Range}ColumnName,
// and the viewer auto-detects them from name suffixes (` min`, ` max`,
// ` range`, ` whisker`, ` error`) when the user does not set X/Y
// columns explicitly. We pin the property bag and the auto-detection
// behavior. We deliberately do not pass {x, y} when exercising
// auto-detect: ScatterPlotLook.auto runs the bulk
// findMatchingErrorBarColumns only when both xColumnName and
// yColumnName are unset (see scatterplot_look.dart `auto()`),
// otherwise the auto path is skipped and the test would falsely fail.
category('AI: GROK-16994: Scatter plot whiskers columns', () => {
  test('auto-detects min+max suffixes for both X and Y axes', async () => {
    // Build numeric columns with the canonical "<base> min" / "<base> max"
    // suffixes on both axes. Categorical noise columns are left out to
    // keep the auto-axis-pick deterministic.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'x', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'x min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      DG.Column.fromList('double', 'x max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      DG.Column.fromList('double', 'y', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'y min', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'y max', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    expect(v instanceof DG.ScatterPlotViewer, true);
    expect(v.type, DG.VIEWER.SCATTER_PLOT);
    // auto() should have picked one numeric base column per axis and
    // populated the whisker-min/max names. The exact pairing of x↔y is
    // implementation-defined, so accept either ordering.
    const xMin = v.props['xWhiskerMinColumnName'];
    const xMax = v.props['xWhiskerMaxColumnName'];
    const yMin = v.props['yWhiskerMinColumnName'];
    const yMax = v.props['yWhiskerMaxColumnName'];
    const xPair = (xMin === 'x min' && xMax === 'x max') ||
      (xMin === 'y min' && xMax === 'y max');
    const yPair = (yMin === 'y min' && yMax === 'y max') ||
      (yMin === 'x min' && yMax === 'x max');
    expect(xPair, true);
    expect(yPair, true);
    // The two axes should have picked different base columns.
    expect(xMin !== yMin, true);
    expect(xMax !== yMax, true);
    // Range stays empty when min+max are both present.
    expect(v.props['xWhiskerRangeColumnName'] == null ||
      v.props['xWhiskerRangeColumnName'] === '', true);
    expect(v.props['yWhiskerRangeColumnName'] == null ||
      v.props['yWhiskerRangeColumnName'] === '', true);
  });

  test('auto-detects range suffix when no min/max columns present', async () => {
    // Two axis-candidates, only y has a "<base> range" column.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'x', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'y', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'y range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    expect(v instanceof DG.ScatterPlotViewer, true);
    // The 'y' base column should resolve a range whisker.
    // It may land on either axis depending on heuristics, so check
    // whichever side picked 'y' as its base.
    const yIsX = v.props['xColumnName'] === 'y';
    const rangeProp = yIsX ? 'xWhiskerRangeColumnName' : 'yWhiskerRangeColumnName';
    const minProp = yIsX ? 'xWhiskerMinColumnName' : 'yWhiskerMinColumnName';
    const maxProp = yIsX ? 'xWhiskerMaxColumnName' : 'yWhiskerMaxColumnName';
    expect(v.props[rangeProp], 'y range');
    expect(v.props[minProp] == null || v.props[minProp] === '', true);
    expect(v.props[maxProp] == null || v.props[maxProp] === '', true);
  });

  test('manual whisker column names with non-conventional names round-trip', async () => {
    // Names that don't match the suffix heuristics — exercises the
    // manual / programmatic API path (no auto-detection).
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'lo', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]),
      DG.Column.fromList('double', 'hi', [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]),
    ]);
    const v = DG.Viewer.scatterPlot(df, {
      x: 'alpha',
      y: 'beta',
      xWhiskerMinColumnName: 'lo',
      xWhiskerMaxColumnName: 'hi',
    });
    expect(v instanceof DG.ScatterPlotViewer, true);
    expect(v.props['xColumnName'], 'alpha');
    expect(v.props['yColumnName'], 'beta');
    expect(v.props['xWhiskerMinColumnName'], 'lo');
    expect(v.props['xWhiskerMaxColumnName'], 'hi');
    const look = v.getOptions(true).look;
    expect(look['xWhiskerMinColumnName'], 'lo');
    expect(look['xWhiskerMaxColumnName'], 'hi');
    // Round-trip via setOptions: swap the pair.
    v.setOptions({
      xWhiskerMinColumnName: 'hi',
      xWhiskerMaxColumnName: 'lo',
    });
    expect(v.props['xWhiskerMinColumnName'], 'hi');
    expect(v.props['xWhiskerMaxColumnName'], 'lo');
    expect(v.getOptions(true).look['xWhiskerMinColumnName'], 'hi');
    expect(v.getOptions(true).look['xWhiskerMaxColumnName'], 'lo');
  });

  test('min+max win over range when both are available for the same axis', async () => {
    // y axis has BOTH a min/max pair AND a range column. The ticket
    // contract says the pair takes priority and range stays unset.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'x', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'y', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'y min', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'y max', [11, 22, 32, 42, 52, 62, 72, 82]),
      DG.Column.fromList('double', 'y range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]);
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    expect(v instanceof DG.ScatterPlotViewer, true);
    // Find which axis took 'y' as its base.
    const yIsX = v.props['xColumnName'] === 'y';
    const yIsY = v.props['yColumnName'] === 'y';
    expect(yIsX || yIsY, true);
    const minProp = yIsX ? 'xWhiskerMinColumnName' : 'yWhiskerMinColumnName';
    const maxProp = yIsX ? 'xWhiskerMaxColumnName' : 'yWhiskerMaxColumnName';
    const rangeProp = yIsX ? 'xWhiskerRangeColumnName' : 'yWhiskerRangeColumnName';
    expect(v.props[minProp], 'y min');
    expect(v.props[maxProp], 'y max');
    // Range must remain unset when the pair won.
    expect(v.props[rangeProp] == null || v.props[rangeProp] === '', true);
  });

  test('tag-based detection: whisker.min / whisker.max tags on axis column', async () => {
    // findErrorColumnsForAxis (utils.dart) iterates the AXIS column's
    // tags and matches keys whose lowercase form contains 'min', 'max',
    // or 'range'. 'whisker.min' / 'whisker.max' satisfy that. Bulk
    // findMatchingErrorBarColumns also has a Phase-1 tag scan with the
    // same isMatchingTag predicate, so attaching the tags to a numeric
    // column makes that column win the auto-axis-pick.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'alpha', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('double', 'beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('double', 'lo1', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]),
      DG.Column.fromList('double', 'hi1', [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]),
      DG.Column.fromList('double', 'lo2', [9, 18, 28, 38, 48, 58, 68, 78]),
      DG.Column.fromList('double', 'hi2', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
    df.col('alpha')!.tags['whisker.min'] = 'lo1';
    df.col('alpha')!.tags['whisker.max'] = 'hi1';
    df.col('beta')!.tags['whisker.min'] = 'lo2';
    df.col('beta')!.tags['whisker.max'] = 'hi2';
    const v = df.plot.scatter() as DG.ScatterPlotViewer;
    expect(v instanceof DG.ScatterPlotViewer, true);
    // alpha and beta both have whisker tags, so each should land on
    // one axis with its tagged whisker columns.
    const xCol = v.props['xColumnName'];
    const yCol = v.props['yColumnName'];
    expect(xCol === 'alpha' || xCol === 'beta', true);
    expect(yCol === 'alpha' || yCol === 'beta', true);
    expect(xCol !== yCol, true);
    const expected = (axisColName: string) => axisColName === 'alpha'
      ? {min: 'lo1', max: 'hi1'} : {min: 'lo2', max: 'hi2'};
    const expX = expected(xCol);
    const expY = expected(yCol);
    expect(v.props['xWhiskerMinColumnName'], expX.min);
    expect(v.props['xWhiskerMaxColumnName'], expX.max);
    expect(v.props['yWhiskerMinColumnName'], expY.min);
    expect(v.props['yWhiskerMaxColumnName'], expY.max);
  });
});

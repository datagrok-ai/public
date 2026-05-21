import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {df, expectCleared, expectPropAndLook} from '../helpers';

// Regression coverage for GROK-16994: Scatter plot whisker columns
// (StdDev / errorbars) — public IScatterPlotSettings now exposes
// xWhisker{Min,Max,Range}ColumnName / yWhisker{Min,Max,Range}ColumnName,
// and the viewer auto-detects them from name suffixes (` min`, ` max`,
// ` range`, ` whisker`, ` error`) when the user does not set X/Y
// columns explicitly. We pin the property bag and the auto-detection
// behavior. We deliberately do not pass {x, y} when exercising
// auto-detect: ScatterPlotLook.auto runs the bulk
// findMatchingErrorBarColumns only when both xColumnName and
// yColumnName are unset, otherwise the auto path is skipped and the
// test would falsely fail.
category('AI: GROK-16994: Scatter plot whiskers columns', () => {
  function nums(name: string, vals: number[]): [string, string, number[]] {
    return [name, 'double', vals];
  }

  test('auto-detects min+max suffixes for both X and Y axes', async () => {
    const v = df([
      nums('x', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('x min', [0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.1]),
      nums('x max', [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9]),
      nums('y', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('y min', [9, 18, 28, 38, 48, 58, 68, 78]),
      nums('y max', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]).plot.scatter() as DG.ScatterPlotViewer;
    const xMin = v.props['xWhiskerMinColumnName'];
    const xMax = v.props['xWhiskerMaxColumnName'];
    const yMin = v.props['yWhiskerMinColumnName'];
    const yMax = v.props['yWhiskerMaxColumnName'];
    const xPair = (xMin === 'x min' && xMax === 'x max') || (xMin === 'y min' && xMax === 'y max');
    const yPair = (yMin === 'y min' && yMax === 'y max') || (yMin === 'x min' && yMax === 'x max');
    expect(xPair, true);
    expect(yPair, true);
    expect(xMin !== yMin && xMax !== yMax, true);
    expectCleared(v.props['xWhiskerRangeColumnName']);
    expectCleared(v.props['yWhiskerRangeColumnName']);
  });

  test('auto-detects range suffix when no min/max columns present', async () => {
    const v = df([
      nums('x', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('y', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('y range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]).plot.scatter() as DG.ScatterPlotViewer;
    const yIsX = v.props['xColumnName'] === 'y';
    const ax = yIsX ? 'x' : 'y';
    expect(v.props[(ax + 'WhiskerRangeColumnName') as keyof typeof v.props], 'y range');
    expectCleared(v.props[(ax + 'WhiskerMinColumnName') as keyof typeof v.props]);
    expectCleared(v.props[(ax + 'WhiskerMaxColumnName') as keyof typeof v.props]);
  });

  test('manual whisker column names with non-conventional names round-trip', async () => {
    const v = DG.Viewer.scatterPlot(df([
      nums('alpha', [1, 2, 3, 4, 5, 6, 7, 8]), nums('beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('lo', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]), nums('hi', [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]),
    ]), {x: 'alpha', y: 'beta', xWhiskerMinColumnName: 'lo', xWhiskerMaxColumnName: 'hi'});
    expectPropAndLook(v, {xColumnName: 'alpha', yColumnName: 'beta', xWhiskerMinColumnName: 'lo', xWhiskerMaxColumnName: 'hi'});
    v.setOptions({xWhiskerMinColumnName: 'hi', xWhiskerMaxColumnName: 'lo'});
    expectPropAndLook(v, {xWhiskerMinColumnName: 'hi', xWhiskerMaxColumnName: 'lo'});
  });

  test('min+max win over range when both are available for the same axis', async () => {
    const v = df([
      nums('x', [1, 2, 3, 4, 5, 6, 7, 8]),
      nums('y', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('y min', [9, 18, 28, 38, 48, 58, 68, 78]),
      nums('y max', [11, 22, 32, 42, 52, 62, 72, 82]),
      nums('y range', [1, 2, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1]),
    ]).plot.scatter() as DG.ScatterPlotViewer;
    const yIsX = v.props['xColumnName'] === 'y';
    const ax = yIsX ? 'x' : 'y';
    expect(v.props[(ax + 'WhiskerMinColumnName') as keyof typeof v.props], 'y min');
    expect(v.props[(ax + 'WhiskerMaxColumnName') as keyof typeof v.props], 'y max');
    expectCleared(v.props[(ax + 'WhiskerRangeColumnName') as keyof typeof v.props]);
  });

  test('tag-based detection: whisker.min / whisker.max tags on axis column', async () => {
    const d = df([
      nums('alpha', [1, 2, 3, 4, 5, 6, 7, 8]), nums('beta', [10, 20, 30, 40, 50, 60, 70, 80]),
      nums('lo1', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]), nums('hi1', [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]),
      nums('lo2', [9, 18, 28, 38, 48, 58, 68, 78]), nums('hi2', [11, 22, 32, 42, 52, 62, 72, 82]),
    ]);
    d.col('alpha')!.tags['whisker.min'] = 'lo1';
    d.col('alpha')!.tags['whisker.max'] = 'hi1';
    d.col('beta')!.tags['whisker.min'] = 'lo2';
    d.col('beta')!.tags['whisker.max'] = 'hi2';
    const v = d.plot.scatter() as DG.ScatterPlotViewer;
    const xCol = v.props['xColumnName'] as string;
    const yCol = v.props['yColumnName'] as string;
    expect(xCol !== yCol, true);
    const exp = (a: string) => a === 'alpha' ? {min: 'lo1', max: 'hi1'} : {min: 'lo2', max: 'hi2'};
    expect(v.props['xWhiskerMinColumnName'], exp(xCol).min);
    expect(v.props['xWhiskerMaxColumnName'], exp(xCol).max);
    expect(v.props['yWhiskerMinColumnName'], exp(yCol).min);
    expect(v.props['yWhiskerMaxColumnName'], exp(yCol).max);
  });
});

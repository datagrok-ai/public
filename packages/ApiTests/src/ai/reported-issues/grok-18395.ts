import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-18395: switching the X axis from a numeric
// column to a string/categorical column on a Line chart caused
// `Unsupported operation: NaN.round()` in Axis.vert during paint. The bug
// only triggers on a real paint, so we attach the viewer to a TableView
// and wait 300ms after the prop swap so the chart can attempt a repaint.
// Dart errors propagate via console only — we do not subscribe to global
// error events; the assertion is "prop swap did not throw + viewer is
// still a LineChartViewer afterwards".
category('AI: GROK-18395: Line chart string X-axis error', () => {
  test('numeric to string X column swap does not throw', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'numCol', [1, 2, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromList('int', 'numX', [10, 20, 30, 40, 50, 60, 70, 80]),
      DG.Column.fromList('string', 'strCol',
        ['a', 'b', 'c', 'a', 'b', 'c', 'a', 'b']),
    ]);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.LINE_CHART,
        {y: 'numCol', xColumnName: 'numX'}) as DG.LineChartViewer;
      expect(v.type, DG.VIEWER.LINE_CHART);
      expect(v instanceof DG.LineChartViewer, true);
      var threw = false;
      try {
        v.props['xColumnName'] = 'strCol';
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      await DG.delay(300);
      expect(v instanceof DG.LineChartViewer, true);
      expect(v.props['xColumnName'], 'strCol');
    }
    finally {
      tv.close();
    }
  });

  test('swap to empty string column does not throw', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'numCol', []),
      DG.Column.fromList('int', 'numX', []),
      DG.Column.fromList('string', 'strCol', []),
    ]);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.LINE_CHART,
        {y: 'numCol', xColumnName: 'numX'}) as DG.LineChartViewer;
      expect(v instanceof DG.LineChartViewer, true);
      var threw = false;
      try {
        v.props['xColumnName'] = 'strCol';
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      await DG.delay(300);
      expect(v instanceof DG.LineChartViewer, true);
      expect(v.props['xColumnName'], 'strCol');
    }
    finally {
      tv.close();
    }
  });

  test('swap to single-value string column does not throw', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'numCol', [1, 2, 3, 4, 5]),
      DG.Column.fromList('int', 'numX', [10, 20, 30, 40, 50]),
      DG.Column.fromList('string', 'strCol',
        ['only', 'only', 'only', 'only', 'only']),
    ]);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.LINE_CHART,
        {y: 'numCol', xColumnName: 'numX'}) as DG.LineChartViewer;
      expect(v instanceof DG.LineChartViewer, true);
      var threw = false;
      try {
        v.props['xColumnName'] = 'strCol';
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      await DG.delay(300);
      expect(v instanceof DG.LineChartViewer, true);
      expect(v.props['xColumnName'], 'strCol');
    }
    finally {
      tv.close();
    }
  });

  test('swap to 1000-row mixed string column does not throw', async () => {
    const n = 1000;
    const nums = new Array<number>(n);
    const xs = new Array<number>(n);
    const cats = new Array<string>(n);
    const labels = ['alpha', 'beta', 'gamma', 'delta', 'epsilon'];
    for (var i = 0; i < n; i++) {
      nums[i] = i;
      xs[i] = i * 2;
      cats[i] = labels[i % labels.length];
    }
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'numCol', nums),
      DG.Column.fromList('int', 'numX', xs),
      DG.Column.fromList('string', 'strCol', cats),
    ]);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.LINE_CHART,
        {y: 'numCol', xColumnName: 'numX'}) as DG.LineChartViewer;
      expect(v instanceof DG.LineChartViewer, true);
      var threw = false;
      try {
        v.props['xColumnName'] = 'strCol';
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      await DG.delay(300);
      expect(v instanceof DG.LineChartViewer, true);
      expect(v.props['xColumnName'], 'strCol');
    }
    finally {
      tv.close();
    }
  });
});

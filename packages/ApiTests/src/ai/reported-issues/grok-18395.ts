import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {df, expectNoThrow, wait, withAttachedViewer} from '../helpers';

// Regression coverage for GROK-18395: switching the X axis from a numeric
// column to a string/categorical column on a Line chart caused
// `Unsupported operation: NaN.round()` in Axis.vert during paint. The bug
// only triggers on a real paint, so we attach the viewer to a TableView
// and wait 300ms after the prop swap so the chart can attempt a repaint.
// Dart errors propagate via console only — we do not subscribe to global
// error events; the assertion is "prop swap did not throw + viewer is
// still a LineChartViewer afterwards".
category('AI: GROK-18395: Line chart string X-axis error', () => {
  async function runSwap(rows: Array<[string, string, any[]]>): Promise<void> {
    await withAttachedViewer<DG.LineChartViewer>(df(rows), DG.VIEWER.LINE_CHART,
      {y: 'numCol', xColumnName: 'numX'}, async (v) => {
        expectNoThrow(() => {v.props['xColumnName'] = 'strCol';});
        await wait();
        expect(v instanceof DG.LineChartViewer, true);
        expect(v.props['xColumnName'], 'strCol');
      });
  }

  test('numeric to string X column swap does not throw', async () =>
    runSwap([
      ['numCol', 'int', [1, 2, 3, 4, 5, 6, 7, 8]],
      ['numX', 'int', [10, 20, 30, 40, 50, 60, 70, 80]],
      ['strCol', 'string', ['a', 'b', 'c', 'a', 'b', 'c', 'a', 'b']],
    ]));

  test('swap to empty string column does not throw', async () =>
    runSwap([['numCol', 'int', []], ['numX', 'int', []], ['strCol', 'string', []]]));

  test('swap to single-value string column does not throw', async () =>
    runSwap([
      ['numCol', 'int', [1, 2, 3, 4, 5]],
      ['numX', 'int', [10, 20, 30, 40, 50]],
      ['strCol', 'string', ['only', 'only', 'only', 'only', 'only']],
    ]));

  test('swap to 1000-row mixed string column does not throw', async () => {
    const n = 1000;
    const nums: number[] = new Array(n);
    const xs: number[] = new Array(n);
    const cats: string[] = new Array(n);
    const labels = ['alpha', 'beta', 'gamma', 'delta', 'epsilon'];
    for (let i = 0; i < n; i++) {
      nums[i] = i;
      xs[i] = i * 2;
      cats[i] = labels[i % labels.length];
    }
    await runSwap([['numCol', 'int', nums], ['numX', 'int', xs], ['strCol', 'string', cats]]);
  });
});

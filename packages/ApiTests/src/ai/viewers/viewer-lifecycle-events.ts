import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {demog, df, expectFiresWithin, expectNoThrow, subscribeAll, withAttachedViewer} from '../helpers';

// Base Viewer lifecycle events
category('AI: Viewers: Lifecycle Events', () => {
  test('base onAfterDrawScene fires on an attached ScatterPlot', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'},
      async (v) => {
        await expectFiresWithin(v.onAfterDrawScene, () => v.invalidateCanvas(), 3000);
      });
  });

  test('base onViewerRendered fires after a re-render on an attached BarChart', async () => {
    await withAttachedViewer<DG.BarChartViewer>(demog(), DG.VIEWER.BAR_CHART, {valueColumnName: 'age'},
      async (v) => {
        await expectFiresWithin(v.onViewerRendered, () => v.setOptions({valueAggrType: 'avg'}), 3000);
      });
  });

  test('base lifecycle getters are Observables inherited by Histogram (subscribe-smoke)', async () => {
    // The remaining base Observables don't deterministically fire on a setOptions toggle headless,
    // so prove base-class inheritance via subscribe-smoke instead.
    await withAttachedViewer<DG.HistogramViewer>(demog(), DG.VIEWER.HISTOGRAM, {value: 'age'}, (v) => {
      expectNoThrow(() => subscribeAll([
        v.onBeforeDrawScene, v.onAfterDrawScene, v.onBeforeDrawOverlay,
        v.onAfterDrawOverlay, v.onAfterLayout, v.onViewerRendered])());
    });
  });

  test('lifecycle getters do not throw on a standalone (un-attached) viewer', async () => {
    const v = df([['a', DG.COLUMN_TYPE.INT, [1, 2, 3]], ['b', DG.COLUMN_TYPE.INT, [3, 2, 1]]])
      .plot.scatter({x: 'a', y: 'b'});
    expectNoThrow(() => subscribeAll([
      v.onBeforeDrawScene, v.onAfterDrawScene, v.onBeforeDrawOverlay,
      v.onAfterDrawOverlay, v.onViewerRendered])());
  });
}, {owner: 'agolovko@datagrok.ai'});

import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, withAttachedViewer} from '../helpers';

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

  test('base onBeforeDrawScene fires before a re-render on an attached ScatterPlot', async () => {
    // onBeforeDrawScene/onAfterDrawScene are defined on the base Viewer; invalidateCanvas() forces a
    // deterministic redraw (Histogram has no base redraw trigger and ignores a no-op fireValuesChanged).
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'},
      async (v) => {
        await expectFiresWithin(v.onBeforeDrawScene, () => v.invalidateCanvas(), 3000);
      });
  });
}, {owner: 'agolovko@datagrok.ai'});

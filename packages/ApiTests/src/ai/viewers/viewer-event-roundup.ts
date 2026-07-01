import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, subscribeAll, withAttachedViewer} from '../helpers';

// Round-up of base-Viewer lifecycle and per-viewer event getters as rxjs Observables.
category('AI: Viewers: Event Roundup', () => {
  test('Histogram missing-values/filter/visibility events are rxjs Observables', async () => {
    await withAttachedViewer<DG.HistogramViewer>(demog(), DG.VIEWER.HISTOGRAM, {value: 'age'}, (v) => {
      expect(v instanceof DG.HistogramViewer, true);
      subscribeAll([v.onMissingValuesFilteredOut, v.onFilterValuesChanged, v.onVisibilityToggled])();
    });
  });

  test('base onDataFrameChanged fires when the DataFrame is replaced', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'},
      async (v) => {
        await expectFiresWithin(v.onDataFrameChanged, () => {
          v.dataFrame = demog(20);
        });
      });
  });

  test('base onResized and onAfterLayout are Observables on an attached ScatterPlot', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'}, (v) => {
      subscribeAll([v.onResized, v.onAfterLayout])();
    });
  });

  test('ScatterPlot onPointsSelected is an rxjs Observable', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'}, (v) => {
      subscribeAll([v.onPointsSelected])();
    });
  });

  test('DensityPlot onZoomed is an rxjs Observable', async () => {
    await withAttachedViewer<DG.DensityPlotViewer>(demog(), DG.VIEWER.DENSITY_PLOT, {x: 'age', y: 'height'}, (v) => {
      expect(v instanceof DG.DensityPlotViewer, true);
      subscribeAll([v.onZoomed])();
    });
  });

  test('TrellisPlot onViewerPropertiesOpen is an rxjs Observable', async () => {
    await withAttachedViewer<DG.TrellisPlotViewer>(demog(), DG.VIEWER.TRELLIS_PLOT, {}, (v) => {
      expect(v instanceof DG.TrellisPlotViewer, true);
      subscribeAll([v.onViewerPropertiesOpen])();
    });
  });
}, {owner: 'agolovko@datagrok.ai'});

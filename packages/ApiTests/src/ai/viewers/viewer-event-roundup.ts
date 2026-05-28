import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, subscribeAll, withAttachedViewer} from '../helpers';

// DG.Viewer / DG.HistogramViewer / DG.ScatterPlotViewer / DG.DensityPlotViewer /
// DG.TrellisPlotViewer — core/client/d4/lib/src/viewer_base/{data_frame_viewer,canvas_viewer_mixin}.dart
// + viewers/{histogram,scatterplot,density_plot,trellis_plot}/*_core.dart
// (scenario: viewer-event-roundup). Round-up of nine event getters reachable via
// onEvent('...'): the three base-Viewer lifecycle streams (onDataFrameChanged,
// onResized, onAfterLayout) plus per-viewer streams that previously had no JS
// binding. onDataFrameChanged fires when the table is replaced; the rest are
// user-interaction streams not headlessly triggerable, so they are proven to be
// well-formed rxjs Observables via subscribeAll.
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

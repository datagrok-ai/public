import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, subscribeAll, until, withAttachedViewer} from '../helpers';

// DG.CorrelationPlot / DG.PivotViewer / DG.TrellisPlotViewer / DG.PcPlot —
// core/client/d4/lib/src/viewers/{correlation_plot,pivot_viewer,trellis_plot,pc_plot}/*_core.dart
// (scenario: viewer-events). Last core-viewer member gap: the stub-viewer event
// observables (onCorrCellClicked, onAggregationChanged, the three TrellisPlot
// events) plus PcPlot.isFiltering. Each event is an onEvent('...') getter — not
// headlessly triggerable — so we prove the typed wrapper is reachable (className
// round-trips through withAttachedViewer) and the streams are well-formed
// Observables (subscribe/unsubscribe), without faking fires. isFiltering is read
// on a fresh and on an attached plot; its true-branch needs a private slider drag
// not reachable from JS, so only the false path is asserted.
category('AI: Viewers: Events', () => {
  test('CorrelationPlot typed wrapper exposes onCorrCellClicked as an rxjs Observable', async () => {
    await withAttachedViewer<DG.CorrelationPlot>(demog(), DG.VIEWER.CORR_PLOT, {}, async (v) => {
      expect(v instanceof DG.CorrelationPlot, true);
      subscribeAll([v.onCorrCellClicked])();
    });
  });

  test('PivotViewer typed wrapper exposes onAggregationChanged as an rxjs Observable', async () => {
    await withAttachedViewer<DG.PivotViewer>(demog(), DG.VIEWER.PIVOT_TABLE, {}, async (v) => {
      expect(v instanceof DG.PivotViewer, true);
      subscribeAll([v.onAggregationChanged])();
    });
  });

  test('TrellisPlotViewer typed wrapper exposes its three event Observables', async () => {
    await withAttachedViewer<DG.TrellisPlotViewer>(demog(), DG.VIEWER.TRELLIS_PLOT, {}, async (v) => {
      expect(v instanceof DG.TrellisPlotViewer, true);
      subscribeAll([v.onInnerViewerClicked, v.onTrellisCurrentCellChanged, v.onViewerTypeChanged])();
    });
  });

  test('event observables are individually well-formed (subscribe/unsubscribe fns)', async () => {
    await withAttachedViewer<DG.CorrelationPlot>(demog(), DG.VIEWER.CORR_PLOT, {}, async (corr) => {
      await withAttachedViewer<DG.PivotViewer>(demog(), DG.VIEWER.PIVOT_TABLE, {}, async (pivot) => {
        await withAttachedViewer<DG.TrellisPlotViewer>(demog(), DG.VIEWER.TRELLIS_PLOT, {}, async (trellis) => {
          subscribeAll([corr.onCorrCellClicked, pivot.onAggregationChanged,
            trellis.onInnerViewerClicked, trellis.onTrellisCurrentCellChanged, trellis.onViewerTypeChanged])();
        });
      });
    });
  });

  test('PcPlot.isFiltering is false on a fresh plot', async () => {
    const v = DG.Viewer.pcPlot(demog(), {columnNames: ['age', 'height', 'weight']});
    expect(v instanceof DG.PcPlot, true);
    expect((v as DG.PcPlot).isFiltering, false);
  });

  test('PcPlot.isFiltering stays false on an attached plot after layout', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT,
      {columnNames: ['age', 'height', 'weight']}, async (v) => {
        await until(() => v.chartBox != null);
        expect(v.isFiltering, false);
      });
  });
}, {owner: 'agolovko@datagrok.ai'});

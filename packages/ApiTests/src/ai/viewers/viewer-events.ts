import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, subscribeAll, until, withAttachedViewer} from '../helpers';

// Stub-viewer event Observables (CorrelationPlot, PivotViewer, TrellisPlot) plus PcPlot.isFiltering.
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

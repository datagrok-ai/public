import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, until, withAttachedViewer} from '../helpers';

// Axis range slider visibility behaviour for ScatterPlot and DensityPlot, asserted on the real slider element.
category('AI: Viewers: Axis Range Sliders', () => {

  // Narrows to the middle 50% so isCustomRange becomes true.
  function narrow(s: DG.RangeSlider): void {
    const lo = s.minRange;
    const hi = s.maxRange;
    const span = hi - lo;
    s.setValues(lo, hi, lo + span * 0.25, hi - span * 0.25);
  }

  test('ScatterPlot: sliders hidden at full range', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        await expectFiresWithin(v.onAfterDrawScene, () => v.invalidateCanvas(), 2000);
        expect(v.xAxisSlider.isCustomRange, false);
        expect(v.yAxisSlider.isCustomRange, false);
        expect(v.xAxisSlider.visible, false);
        expect(v.yAxisSlider.visible, false);
      });
  });

  test('ScatterPlot: Y slider becomes visible when narrowed (default layout)', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        narrow(v.yAxisSlider);
        await expectFiresWithin(v.onAfterDrawScene, () => v.invalidateCanvas(), 2000);
        expect(v.yAxisSlider.isCustomRange, true);
        expect(v.yAxisSlider.visible, true);
      });
  });

  test('ScatterPlot: showYAxis=false hides the Y slider even when narrowed', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        narrow(v.yAxisSlider);
        v.setOptions({showYAxis: false});
        await expectFiresWithin(v.onAfterDrawScene, () => v.invalidateCanvas(), 2000);
        expect(v.yAxisSlider.isCustomRange, true);
        expect(v.yAxisSlider.visible, false);
      });
  });

  test('ScatterPlot: autoLayout=false keeps the Y slider visible regardless of showYAxis', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        v.setOptions({autoLayout: false, showYAxis: false});
        narrow(v.yAxisSlider);
        await expectFiresWithin(v.onAfterDrawScene, () => v.invalidateCanvas(), 2000);
        expect(v.yAxisSlider.visible, true);
      });
  });

  test('ScatterPlot: X and Y slider visibility are independent', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        narrow(v.xAxisSlider);
        await expectFiresWithin(v.onAfterDrawScene, () => v.invalidateCanvas(), 2000);
        expect(v.xAxisSlider.isCustomRange, true);
        expect(v.xAxisSlider.visible, true);
        expect(v.yAxisSlider.isCustomRange, false);
        expect(v.yAxisSlider.visible, false);
      });
  });

  test('DensityPlot: same render-time visibility rule (hidden at full, visible when narrowed)', async () => {
    await withAttachedViewer<DG.DensityPlotViewer>(demog(), DG.VIEWER.DENSITY_PLOT,
      {xColumnName: 'age', yColumnName: 'height'}, async (v) => {
        await until(() => v.xAxisSlider != null);
        expect(v.xAxisSlider.visible, false);
        narrow(v.xAxisSlider);
        // DensityPlot has no invalidateCanvas; setOptions forces a re-render.
        v.setOptions({showXAxis: true});
        await until(() => v.xAxisSlider.visible === true);
        expect(v.xAxisSlider.isCustomRange, true);
        expect(v.xAxisSlider.visible, true);
      });
  });
}, {owner: 'agolovko@datagrok.ai'});

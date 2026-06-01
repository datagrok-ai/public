import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, until, withAttachedViewer} from '../helpers';

// Axis range slider VISIBILITY behaviour, asserted on the real slider element.
//
// ScatterPlot and DensityPlot set slider visibility during RENDER:
//   htmlSetVisible([slider.g], isCustomRange && (!autoLayout || showAxis))
//   (scatterplot_core.dart:1264/1266, density_plot_core.dart:402/405)
// so the outcome is observable headlessly (grok test runs in a real browser):
// we narrow a slider / flip props, force a render, and assert RangeSlider.visible
// against the full rule — including BOTH gating terms (showXAxis/showYAxis and
// autoLayout), per axis, independently.
//
// NOT asserted (documented):
// - BarChart's axis-slider visibility is driven ONLY by mouse enter/leave
//   (bar_chart_core.dart:1296-1305), not by render — so it is not testable
//   headlessly without synthesizing mouse events. Its determinant is still
//   isCustomRange (+ showValueAxis for the value axis).
// - The hover-to-reveal-at-full-range behaviour (slider shows on mouse-in even
//   at full range, gated by D4Settings.showRangeSlidersOnViewers == 'Auto') is
//   mouse-driven and likewise not asserted here.
category('AI: Viewers: Axis Range Sliders', () => {

  // Narrow the selected window to the middle 50% of the slider's full range,
  // keeping the range bounds: makes isCustomRange true. notify -> viewer re-render.
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
        // full range => isCustomRange false => visibility expression false
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
        // isCustomRange true, autoLayout true (default), showYAxis true (default)
        // => true && (false || true) => visible
        expect(v.yAxisSlider.isCustomRange, true);
        expect(v.yAxisSlider.visible, true);
      });
  });

  test('ScatterPlot: showYAxis=false hides the Y slider even when narrowed', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        narrow(v.yAxisSlider);
        v.setOptions({showYAxis: false}); // autoLayout stays true
        await expectFiresWithin(v.onAfterDrawScene, () => v.invalidateCanvas(), 2000);
        // isCustomRange true, but (!autoLayout || showYAxis) = (false || false) = false
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
        // isCustomRange true, (!autoLayout || showYAxis) = (true || false) = true => visible
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
        // Y untouched -> full range -> hidden
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

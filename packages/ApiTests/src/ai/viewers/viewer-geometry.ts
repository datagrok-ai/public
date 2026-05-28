import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, wait, withAttachedViewer} from '../helpers';

// DG.ScatterPlotViewer / DG.LineChartViewer — core/client/d4/lib/src/viewers/scatterplot/*.dart,
// core/client/d4/lib/src/viewers/line_chart/*.dart (scenario: viewer-geometry).
// Smoke-breadth coverage of the canvas-geometry surface that only resolves after
// the first render: ScatterPlot xAxisBox/yAxisBox (Rect), worldToScreen/pointToScreen
// (Point), hitTest (row index or -1), render into an offscreen 2d context, and
// LineChart's chart-indexed worldToScreen. Getters are null pre-first-render, so
// every case attaches the viewer via withAttachedViewer and waits for layout.
category('AI: Viewers: Geometry', () => {
  test('xAxisBox / yAxisBox are sane Rects after layout', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        await wait(200);
        const xb = v.xAxisBox;
        const yb = v.yAxisBox;
        expect(xb != null, true);
        expect(yb != null, true);
        expect(xb.width >= 0, true);
        expect(xb.height >= 0, true);
        expect(yb.width >= 0, true);
        expect(yb.height >= 0, true);
      });
  });

  test('worldToScreen returns a finite Point', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        await wait(200);
        const p = v.worldToScreen(40, 170);
        expect(p != null, true);
        expect(Number.isFinite(p.x), true);
        expect(Number.isFinite(p.y), true);
      });
  });

  test('pointToScreen returns a Point for row 0', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        await wait(200);
        const p = v.pointToScreen(0);
        expect(p != null, true);
        // Clipped rows may map to NaN — assert numeric type, not finiteness.
        expect(typeof p.x, 'number');
        expect(typeof p.y, 'number');
      });
  });

  test('hitTest returns a row index or -1', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        await wait(200);
        const p = v.worldToScreen(40, 170);
        const hit = v.hitTest(p.x, p.y);
        expect(typeof hit, 'number');
        expect(hit === -1 || hit >= 0, true);
        expect(v.hitTest(-10000, -10000), -1);
      });
  });

  test('render into an offscreen 2d context does not throw', async () => {
    await withAttachedViewer<DG.ScatterPlotViewer>(demog(), DG.VIEWER.SCATTER_PLOT,
      {x: 'age', y: 'height'}, async (v) => {
        await wait(200);
        const g = document.createElement('canvas').getContext('2d')!;
        expectNoThrow(() => v.render(g));
      });
  });

  test('LineChart worldToScreen(x, y, chartIdx) maps a point for chart 0', async () => {
    await withAttachedViewer<DG.LineChartViewer>(demog(), DG.VIEWER.LINE_CHART,
      {x: 'age', yColumnNames: ['height']}, async (v) => {
        await wait(200);
        const p = v.worldToScreen(40, 170, 0);
        expect(p != null, true);
        expect(Number.isFinite(p.x), true);
        expect(Number.isFinite(p.y), true);
      // NOTE: an out-of-range chartIdx indexes the chart array and throws — not a
      // tolerated boundary, so it's intentionally not exercised here.
      });
  });
}, {owner: 'agolovko@datagrok.ai'});

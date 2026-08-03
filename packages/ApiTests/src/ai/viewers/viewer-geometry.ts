import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, until, withAttachedViewer} from '../helpers';

// Canvas-geometry surface (axis boxes, worldToScreen, hitTest, render) after first render.
// The ScatterPlot checks are read-only, so they share one attached+rendered plot (clear:false);
// the LineChart case keeps its own viewer.
category('AI: Viewers: Geometry', () => {
  let src: DG.DataFrame;
  let tv: DG.TableView;
  let v: DG.ScatterPlotViewer;

  before(async () => {
    src = demog(50);
    tv = grok.shell.addTableView(src);
    v = tv.addViewer(DG.VIEWER.SCATTER_PLOT, {x: 'age', y: 'height'}) as DG.ScatterPlotViewer;
    await until(() => v.xAxisBox != null);
  });

  after(async () => {
    tv.close();
  });

  test('xAxisBox / yAxisBox are sane Rects after layout', async () => {
    const xb = v.xAxisBox;
    const yb = v.yAxisBox;
    expect(xb != null, true);
    expect(yb != null, true);
    expect(xb.width >= 0, true);
    expect(xb.height >= 0, true);
    expect(yb.width >= 0, true);
    expect(yb.height >= 0, true);
  });

  test('worldToScreen returns a finite Point', async () => {
    const p = v.worldToScreen(40, 170);
    expect(p != null, true);
    expect(Number.isFinite(p.x), true);
    expect(Number.isFinite(p.y), true);
  });

  test('pointToScreen returns a Point for row 0', async () => {
    const p = v.pointToScreen(0);
    expect(p != null, true);
    // Clipped rows may map to NaN — assert numeric type, not finiteness.
    expect(typeof p.x, 'number');
    expect(typeof p.y, 'number');
  });

  test('hitTest returns a row index or -1', async () => {
    const p = v.worldToScreen(40, 170);
    const hit = v.hitTest(p.x, p.y);
    expect(typeof hit, 'number');
    expect(hit === -1 || hit >= 0, true);
    expect(v.hitTest(-10000, -10000), -1);
  });

  test('render into an offscreen 2d context does not throw', async () => {
    const g = document.createElement('canvas').getContext('2d')!;
    expectNoThrow(() => v.render(g));
  });

  test('screenToWorld inverts worldToScreen within tolerance', async () => {
    const screen = v.worldToScreen(40, 170);
    const world = v.screenToWorld(screen.x, screen.y);
    expect(world != null, true);
    expect(Number.isFinite(world.x), true);
    expect(Number.isFinite(world.y), true);
    expect(Math.abs(world.x - 40) < 5, true);
    expect(Math.abs(world.y - 170) < 5, true);
  });

  test('marker accessors return per-row sizes, colors, and types', async () => {
    await until(() => v.getMarkerColors().length === src.rowCount);
    expect(v.getMarkerSizes().length, src.rowCount);
    expect(v.getMarkerColors().length, src.rowCount);
    expect(v.getMarkerTypes().length, src.rowCount);
    expect(Number.isFinite(v.getMarkerSize(0)), true);
    expect(v.getMarkerSize(0) > 0, true);
    expect(typeof v.getMarkerColor(0), 'number');
    expect(typeof v.getMarkerType(0), 'string');
  });

  test('getRowTooltip returns an HTML element for a row', async () => {
    expect(v.getRowTooltip(0) instanceof HTMLElement, true);
  });

  test('LineChart worldToScreen(x, y, chartIdx) maps a point for chart 0', async () => {
    await withAttachedViewer<DG.LineChartViewer>(demog(), DG.VIEWER.LINE_CHART,
      {x: 'age', yColumnNames: ['height']}, async (lc) => {
        await until(() => lc.worldToScreen(40, 170, 0) != null);
        const p = lc.worldToScreen(40, 170, 0);
        expect(p != null, true);
        expect(Number.isFinite(p.x), true);
        expect(Number.isFinite(p.y), true);
        // Out-of-range chartIdx throws, so it's intentionally not exercised here.
      });
  });
}, {owner: 'agolovko@datagrok.ai', clear: false});

/* ---
realizes: [scatterplot.cp.select-and-zoom, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as sp from './scatterplot-shared';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const JITTER_X = 20;
const JITTER_Y = 15;
const JITTER_X_CHANGED = 30;

async function hoverPlot(page: Page): Promise<void> {
  const r = await sp.canvasRect(page);
  const shown = await v.armEvent(page, 'grok.events.onTooltipShown', 200);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  await shown();
}

const viewerProps = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {
    x: s.props.xColumnName, y: s.props.yColumnName,
    zoomAndFilter: s.props.zoomAndFilter,
    jitter: s.props.jitterSize, jitterY: s.props.jitterSizeY,
    lasso: s.props.lassoTool,
    resetOnBackgroundClick: s.props.resetSelectionOnBackgroundClick,
  };
});

function expectSameViewport(actual: sp.Rect, reference: sp.Rect, tolerance = 0.02): void {
  expect(Math.abs(actual.width - reference.width)).toBeLessThan(reference.width * tolerance);
  expect(Math.abs(actual.height - reference.height)).toBeLessThan(reference.height * tolerance);
  expect(Math.abs(actual.x - reference.x)).toBeLessThan(reference.width * tolerance);
  expect(Math.abs(actual.y - reference.y)).toBeLessThan(reference.height * tolerance);
}

const rowCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);

async function lassoCanvas(page: Page, points: sp.Frac[], mods: string[] = ['Shift']): Promise<void> {
  const r = await sp.canvasRect(page);
  const path = points.map((p) => sp.at(r, p));
  const shown1 = await v.armEvent(page, 'grok.events.onTooltipShown', 100);
  await page.mouse.move(path[0].x, path[0].y);
  await shown1();
  for (const m of mods) await page.keyboard.down(m);
  await page.mouse.down();
  // a lasso is a PACED sequence of moves: the pacing is the gesture, not a wait
  for (const p of path.slice(1)) {
    await page.mouse.move(p.x, p.y, {steps: 6});
    await page.waitForTimeout(40);
  }
  const shown2 = await v.armEvent(page, 'grok.events.onTooltipShown', 150);
  await page.mouse.move(path[0].x, path[0].y, {steps: 6});
  await shown2();
  await page.mouse.up();
  for (const m of [...mods].reverse()) await page.keyboard.up(m);
  await v.waitForViewerRendered(page, sp.SP_TYPE, 400);
}

async function wheelZoomIn(page: Page, from: sp.Rect, steps = 1): Promise<sp.Rect> {
  const r = await sp.canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  for (let i = 0; i < steps; i++) {
    await page.mouse.wheel(0, -300);
    await v.waitForViewerRendered(page, sp.SP_TYPE, 300);
  }
  return sp.viewportMoved(page, from);
}

async function resetViewFromMenu(page: Page, from?: sp.Rect): Promise<void> {
  await sp.openPlotContextMenu(page);
  const item = page.locator('.d4-menu-popup [name="div-Reset-View"]').last();
  await item.waitFor({state: 'visible', timeout: 5000});
  await item.click();
  await sp.dismissMenu(page);
  if (from) await sp.viewportMoved(page, from);
  else await v.waitForViewerQuiet(page, sp.SP_TYPE, {gapMs: 200, capMs: 1500});
}

// A point of the plot field with no marker within `clearancePx`, computed from the rows'
// screen positions rather than probed by trial hovers.
const markerFreePoint = (page: Page, clearancePx = 12) => page.evaluate((clearance: number) => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const df = grok.shell.tv.dataFrame;
  const xc = df.col(s.props.xColumnName);
  const yc = df.col(s.props.yColumnName);
  const canvas = s.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
  const b = canvas.getBoundingClientRect();
  const pts: {x: number; y: number}[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    if (!df.filter.get(i) || xc.isNone(i) || yc.isNone(i)) continue;
    const p = s.worldToScreen(xc.get(i), yc.get(i));
    if (p.x >= 0 && p.y >= 0 && p.x <= b.width && p.y <= b.height) pts.push(p);
  }
  const step = 10;
  for (let fy = 0.12; fy <= 0.88; fy += 0.04) {
    for (let fx = 0.12; fx <= 0.88; fx += 0.04) {
      const x = b.width * fx;
      const y = b.height * fy;
      let clear = true;
      for (const p of pts)
        if (Math.abs(p.x - x) < clearance + step && Math.abs(p.y - y) < clearance + step) { clear = false; break; }
      if (clear) return {x: b.x + x, y: b.y + y};
    }
  }
  return null;
}, clearancePx);

async function focusCanvas(page: Page): Promise<void> {
  const r = await sp.canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2);
  await v.waitForViewerRendered(page, sp.SP_TYPE, 200);
}

const sliderHandlePoint = (page: Page, slider: string, handle: string) =>
  page.evaluate(({s, h}: {s: string; h: string}) => {
    const root = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot').root as HTMLElement;
    const el = root.querySelector(`svg[name="${s}"] [name="${h}"]`) as SVGGraphicsElement | null;
    if (!el) return null;
    const b = el.getBoundingClientRect();
    if (b.width === 0 || b.height === 0) return null;
    return {x: b.x + b.width / 2, y: b.y + b.height / 2};
  }, {s: slider, h: handle});

test('Scatter Plot — Point Selection and Viewport Navigation', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const errors = sp.trackErrors(page);
  const errCount = errors.count;

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await sp.addScatterPlot(page);
  await sp.pickOnViewer(page, 'x', 'WEIGHT');
  await sp.pickOnViewer(page, 'y', 'HEIGHT');
  const fullRowCount = await rowCount(page);
  expect(fullRowCount).toBeGreaterThan(0);

  let filteredRows = fullRowCount;

  await softStep('Select points by dragging, add a band, deselect, clear', async () => {
    const errBefore = errCount();
    expect(await sp.selectionHeld(page)).toBe(0);

    await sp.dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const first = await sp.selectionMoved(page, 0);
    expect(first).toBeGreaterThan(0);

    await sp.dragCanvas(page, {fx: 0.05, fy: 0.05}, {fx: 0.27, fy: 0.95}, ['Shift']);
    const second = await sp.selectionMoved(page, first);
    expect(second).toBeGreaterThan(first);

    await sp.dragCanvas(page, {fx: 0.35, fy: 0.35}, {fx: 0.55, fy: 0.55}, ['Control', 'Shift']);
    const deselected = await sp.selectionMoved(page, second);
    expect(deselected).toBeLessThan(second);
    expect(deselected).toBeGreaterThan(0);

    expect((await viewerProps(page)).resetOnBackgroundClick).toBe(true);
    await sp.clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await sp.selectionMoved(page, deselected)).toBe(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Selection survives a jitter change', async () => {
    const errBefore = errCount();
    await sp.setNumericProp(page, 'prop-jitter-size', 'marker', JITTER_X);
    await sp.setNumericProp(page, 'prop-jitter-size-y', 'marker', JITTER_Y);
    const jittered = await viewerProps(page);
    expect(jittered.jitter).toBe(JITTER_X);
    expect(jittered.jitterY).toBe(JITTER_Y);

    await sp.dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const selected = await sp.selectionMoved(page, 0);
    expect(selected).toBeGreaterThan(0);

    await sp.setNumericProp(page, 'prop-jitter-size', 'marker', JITTER_X_CHANGED);
    expect((await viewerProps(page)).jitter).toBe(JITTER_X_CHANGED);
    await v.waitForViewerQuiet(page, sp.SP_TYPE, {gapMs: 200, capMs: 1500});

    expect(await sp.selectionHeld(page)).toBe(selected);

    await sp.dragCanvas(page, {fx: 0.35, fy: 0.35}, {fx: 0.55, fy: 0.55}, ['Control', 'Shift']);
    const narrowed = await sp.selectionMoved(page, selected);
    expect(narrowed).toBeLessThan(selected);
    expect(errCount()).toBe(errBefore);

    await sp.clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await sp.selectionMoved(page, narrowed)).toBe(0);
    await sp.setNumericProp(page, 'prop-jitter-size', 'marker', 0);
    await sp.setNumericProp(page, 'prop-jitter-size-y', 'marker', 0);
    const reverted = await viewerProps(page);
    expect(reverted.jitter).toBe(0);
    expect(reverted.jitterY).toBe(0);
  });

  await softStep('Zoom, pan and reset the viewport', async () => {
    const errBefore = errCount();
    if (await sp.readProp(page, 'zoomAndFilter') !== 'no action')
      await sp.setChoiceProp(page, 'prop-zoom-and-filter', 'data', 'no action');
    expect(await sp.readProp(page, 'zoomAndFilter')).toBe('no action');

    await resetViewFromMenu(page);
    const baseline = await sp.viewport(page);
    expect(baseline.width).toBeGreaterThan(0);
    expect(baseline.height).toBeGreaterThan(0);

    await sp.dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.6, fy: 0.6}, ['Alt']);
    const zoomed = await sp.viewportMoved(page, baseline);
    expect(zoomed.width).toBeLessThan(baseline.width);
    expect(zoomed.height).toBeLessThan(baseline.height);

    expect(await sp.filterCount(page)).toBe(fullRowCount);

    await sp.dragCanvas(page, {fx: 0.5, fy: 0.5}, {fx: 0.3, fy: 0.3});
    const panned = await sp.viewportMoved(page, zoomed);
    expect(Math.abs(panned.x - zoomed.x)).toBeGreaterThan(zoomed.width * 0.02);
    expect(Math.abs(panned.width - zoomed.width)).toBeLessThan(zoomed.width * 0.02);

    await resetViewFromMenu(page, panned);
    expectSameViewport(await sp.viewport(page), baseline);

    const wheeled = await wheelZoomIn(page, baseline);
    expect(wheeled.width).toBeLessThan(baseline.width);
    expect(wheeled.height).toBeLessThan(baseline.height);
    await resetViewFromMenu(page, wheeled);
    expectSameViewport(await sp.viewport(page), baseline);

    await hoverPlot(page);
    const handle = await v.pollValue(() => sliderHandlePoint(page, 'x-slider', 'min-handle'), (p) => p !== null, 3000, 50);
    expect(handle).not.toBeNull();
    const beforeSlider = await sp.viewport(page);

    const shown1 = await v.armEvent(page, 'grok.events.onTooltipShown', 150);
    await page.mouse.move(handle!.x, handle!.y);
    await shown1();
    await page.mouse.down();
    await page.mouse.move(handle!.x + 15, handle!.y, {steps: 6});
    const shown2 = await v.armEvent(page, 'grok.events.onTooltipShown', 150);
    await page.mouse.move(handle!.x + 30, handle!.y, {steps: 6});
    await shown2();
    await page.mouse.up();
    const afterSlider = await sp.viewportMoved(page, beforeSlider);

    expect(afterSlider.x).toBeGreaterThan(beforeSlider.x);
    expect(afterSlider.width).toBeLessThan(beforeSlider.width);

    await resetViewFromMenu(page, afterSlider);
    expectSameViewport(await sp.viewport(page), baseline);

    const wheeledAgain = await wheelZoomIn(page, baseline);
    expect(wheeledAgain.width).toBeLessThan(baseline.width);
    const gap = await markerFreePoint(page);
    expect(gap).not.toBeNull();
    await page.mouse.move(gap!.x, gap!.y);
    expect(await v.pollValue(() => page.evaluate(() => grok.shell.tv.dataFrame.mouseOverRowIdx as number),
      (r) => r === -1, 1000, 50)).toBe(-1);
    await page.mouse.dblclick(gap!.x, gap!.y);
    expectSameViewport(await sp.viewportMoved(page, wheeledAgain), baseline);
    expect(errCount()).toBe(errBefore);

    await sp.setChoiceProp(page, 'prop-zoom-and-filter', 'data', 'filter by zoom');
    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
  });

  await softStep('Keyboard selection and view shortcuts', async () => {
    const errBefore = errCount();
    await sp.openFilterPanel(page);
    await v.applyCategoricalFilter(page, 'SEX', ['F'], 300);
    const filtered = await sp.filterMoved(page, fullRowCount);
    filteredRows = filtered;
    expect(filtered).toBeLessThan(fullRowCount);
    expect(filtered).toBeGreaterThan(0);

    await focusCanvas(page);
    await page.keyboard.press('Control+a');
    expect(await sp.selectionMoved(page, 0)).toBe(filtered);

    await page.keyboard.press('Control+Shift+a');
    expect(await sp.selectionMoved(page, filtered)).toBe(0);

    await sp.dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const dragged = await sp.selectionMoved(page, 0);
    expect(dragged).toBeGreaterThan(0);
    await page.keyboard.press('Escape');
    expect(await sp.selectionMoved(page, dragged)).toBe(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Keyboard selection and view shortcuts — H restores the view', async () => {
    const errBefore = errCount();
    await resetViewFromMenu(page);
    const beforeZoom = await sp.viewport(page);

    const zoomed = await wheelZoomIn(page, beforeZoom);
    expect(zoomed.width).toBeLessThan(beforeZoom.width);

    await focusCanvas(page);
    await page.keyboard.press('h');
    expectSameViewport(await sp.viewportMoved(page, zoomed), beforeZoom);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Keyboard selection and view shortcuts — the Lasso Tool', async () => {
    const errBefore = errCount();
    await focusCanvas(page);
    await page.keyboard.press('l');
    expect(await sp.propIs(page, 'lassoTool', true, 4000)).toBe(true);

    await sp.clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await sp.selectionHeld(page)).toBe(0);
    await lassoCanvas(page, [
      {fx: 0.35, fy: 0.35}, {fx: 0.65, fy: 0.32}, {fx: 0.7, fy: 0.62},
      {fx: 0.5, fy: 0.72}, {fx: 0.33, fy: 0.6},
    ]);
    const lassoed = await sp.selectionMoved(page, 0);
    expect(lassoed).toBeGreaterThan(0);

    await focusCanvas(page);
    await page.keyboard.press('l');
    expect(await sp.propIs(page, 'lassoTool', false, 4000)).toBe(false);
    expect(errCount()).toBe(errBefore);

    await sp.clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await sp.selectionMoved(page, lassoed)).toBe(0);
    await page.locator('[name="viewer-Filters"] [name="icon-arrow-rotate-left"]').click();
    expect(await sp.filterMoved(page, filteredRows)).toBe(fullRowCount);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});

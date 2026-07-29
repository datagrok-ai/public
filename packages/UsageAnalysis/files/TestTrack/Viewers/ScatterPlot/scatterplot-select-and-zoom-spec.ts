/* ---
realizes: [scatterplot.cp.select-and-zoom, viewers.scatter-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const JITTER_X = 20;
const JITTER_Y = 15;
const JITTER_X_CHANGED = 30;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

// This file opens no window in which a Dart NullError or a "Stack trace" line is
// expected, so nothing beyond the shared-server ambient noise is whitelisted —
// such a message during a selection, a jitter change or a viewport operation is
// exactly the regression signal these guards look for.
const isBenignError = (text: string) => isAmbientError(text);

interface Rect {x: number; y: number; width: number; height: number}

/** Rect of the scatter-plot data canvas, resolved on the viewer root that is
 * NOT inside a dialog (a dialog can embed its own preview plot). */
const canvasRect = (page: Page): Promise<Rect> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

/** Move the real pointer over the plot. The range-slider handles carry
 * `visibility: hidden` until the pointer is over the viewer, so the slider drag
 * starts here. */
async function hoverPlot(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  await page.waitForTimeout(200);
}

const pickOnViewer = (page: Page, role: string, column: string) =>
  v.pickColumnViaSelectorTrusted(page, {role, columnName: column});

/** Open the viewer's settings panel and wait for its property grid to build. */
async function openSettings(page: Page): Promise<void> {
  for (let i = 0; i < 4; i++) {
    const built = await page.evaluate(() => !!document.querySelector('[name="prop-category-data"]'));
    if (built) return;
    await v.openViewerGear(page, 'Scatter plot');
    await page.waitForFunction(() => !!document.querySelector('[name="prop-category-data"]'),
      null, {timeout: 2500}).catch(() => {});
  }
  throw new Error('the scatter plot settings panel did not build');
}

/** Make a property editor reachable. A row inside a collapsed category keeps its
 * DOM node with an empty box, so readiness is measured on the editor's own
 * rectangle; the category header is a toggle, so an attempt is simply retried. */
async function revealPropEditor(page: Page, editorSelector: string, category: string): Promise<void> {
  for (let i = 0; i < 8; i++) {
    const ready = await page.waitForFunction((sel: string) => {
      const el = document.querySelector(sel) as HTMLElement | null;
      if (!el || !el.offsetParent) return false;
      const b = el.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    }, editorSelector, {timeout: i === 0 ? 300 : 900}).then(() => true).catch(() => false);
    if (ready) return;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();
  }
  throw new Error(`property editor ${editorSelector} never became reachable`);
}

/** Set a choice property from the settings panel: the value cell turns into a
 * select once clicked. */
async function setChoiceProp(
  page: Page, rowName: string, viewCell: string, category: string, value: string,
): Promise<void> {
  await openSettings(page);
  await revealPropEditor(page, `[name="${viewCell}"]`, category);
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  const editor = page.locator(`[name="${rowName}"] select.property-grid-item-editor-spinner`);
  await editor.waitFor({state: 'visible', timeout: 4000});
  await editor.selectOption(value);
  for (let i = 0; i < 25; i++) {
    if ((await propCellText(page, viewCell)).toLowerCase() === value.toLowerCase()) break;
    await page.waitForTimeout(100);
  }
  await page.waitForTimeout(200);
}

/** Property name a settings-panel row edits: `prop-jitter-size-y` → `jitterSizeY`. */
const propNameOfRow = (rowName: string) =>
  rowName.replace(/^prop-/, '').replace(/-([a-z])/g, (_, c: string) => c.toUpperCase());

/** Set a numeric slider property from its textbox. */
async function setSliderProp(page: Page, rowName: string, category: string, value: number): Promise<void> {
  await openSettings(page);
  const box = `[name="${rowName}"] input.property-grid-slider-textbox`;
  await revealPropEditor(page, box, category);
  const locator = page.locator(box);
  await locator.scrollIntoViewIfNeeded();
  await locator.click();
  await page.keyboard.press('Control+A');
  await page.keyboard.type(String(value));
  await page.keyboard.press('Enter');
  const prop = propNameOfRow(rowName);
  for (let i = 0; i < 25; i++) {
    const cur = await page.evaluate((n: string) => {
      const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
      return sp.props[n];
    }, prop);
    if (cur === value) break;
    await page.waitForTimeout(100);
  }
  await page.waitForTimeout(200);
}

/** Text the settings panel shows for a property's current value. */
const propCellText = (page: Page, viewCell: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`)?.textContent ?? '').trim(), viewCell);

const viewerProps = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {
    x: sp.props.xColumnName, y: sp.props.yColumnName,
    zoomAndFilter: sp.props.zoomAndFilter,
    jitter: sp.props.jitterSize, jitterY: sp.props.jitterSizeY,
    lasso: sp.props.lassoTool,
    resetOnBackgroundClick: sp.props.resetSelectionOnBackgroundClick,
  };
});

/** The viewer's own viewport rectangle — the product-state signal for every
 * zoom, pan and reset step. `props.viewport` is null on this viewer and is not
 * a substitute. */
const viewport = (page: Page): Promise<Rect> => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const vp = sp.viewport;
  return {x: vp.x, y: vp.y, width: vp.width, height: vp.height};
});

/** Two viewport rectangles are the same view when every side agrees to within a
 * small fraction of the reference rectangle — a reset recomputes the rect and
 * is not required to be bit-exact. */
function expectSameViewport(actual: Rect, reference: Rect, tolerance = 0.02): void {
  expect(Math.abs(actual.width - reference.width)).toBeLessThan(reference.width * tolerance);
  expect(Math.abs(actual.height - reference.height)).toBeLessThan(reference.height * tolerance);
  expect(Math.abs(actual.x - reference.x)).toBeLessThan(reference.width * tolerance);
  expect(Math.abs(actual.y - reference.y)).toBeLessThan(reference.height * tolerance);
}

const rowCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);
const filterCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount as number);
const selectionCount = (page: Page) =>
  page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount as number);

const sameRect = (a: Rect, b: Rect) =>
  a.x === b.x && a.y === b.y && a.width === b.width && a.height === b.height;

/** The viewport once two consecutive readings agree — a view operation is
 * applied asynchronously and this is the anchor every viewport assert waits on. */
async function settledViewport(page: Page): Promise<Rect> {
  let prev = await viewport(page);
  for (let i = 0; i < 24; i++) {
    await page.waitForTimeout(150);
    const cur = await viewport(page);
    if (sameRect(cur, prev)) return cur;
    prev = cur;
  }
  return prev;
}

/** The viewport once it has left a known rectangle and stopped moving. Used
 * where the step must move the view: a rectangle that never leaves `from` comes
 * back unchanged and the caller's assert fails on it. */
async function settledViewportAfterChange(page: Page, from: Rect): Promise<Rect> {
  let prev = from;
  for (let i = 0; i < 50; i++) {
    await page.waitForTimeout(120);
    const cur = await viewport(page);
    if (sameRect(cur, from)) continue;
    if (sameRect(cur, prev)) return cur;
    prev = cur;
  }
  return prev;
}

/** Selected-row count, read once it stops moving: the selection is recomputed
 * asynchronously after a gesture. Used where the step is asserted NOT to change
 * the count, so there is no known value to wait away from. */
async function settledSelection(page: Page): Promise<number> {
  let last = -1;
  let stable = 0;
  for (let i = 0; i < 30; i++) {
    const c = await selectionCount(page);
    if (c === last) {
      stable++;
      if (stable >= 2) return c;
    } else {
      stable = 0;
      last = c;
    }
    await page.waitForTimeout(200);
  }
  return last;
}

/** Selected-row count once it has left a known value and stopped moving. A count
 * that never leaves `from` is returned as it stands, so the caller's assert
 * fails on the real reading. */
async function settledSelectionAfterChange(page: Page, from: number): Promise<number> {
  let prev = from;
  for (let i = 0; i < 60; i++) {
    const c = await selectionCount(page);
    if (c !== from) {
      if (c === prev) return c;
      prev = c;
    }
    await page.waitForTimeout(120);
  }
  return await selectionCount(page);
}

/** Filtered-row count once it has left a known value and stopped moving. */
async function settledFilterAfterChange(page: Page, from: number): Promise<number> {
  let prev = from;
  for (let i = 0; i < 60; i++) {
    const c = await filterCount(page);
    if (c !== from) {
      if (c === prev) return c;
      prev = c;
    }
    await page.waitForTimeout(120);
  }
  return await filterCount(page);
}

interface Frac {fx: number; fy: number}

const at = (r: Rect, p: Frac) => ({x: r.x + r.width * p.fx, y: r.y + r.height * p.fy});

/** Drag on the data canvas with the real pointer, optionally holding modifier
 * keys. Intermediate moves are required — the viewer tracks the gesture through
 * its own move handler and a straight down/up pair reads as a click. */
async function dragCanvas(page: Page, from: Frac, to: Frac, mods: string[] = []): Promise<void> {
  const r = await canvasRect(page);
  const p1 = at(r, from);
  const p2 = at(r, to);
  await page.mouse.move(p1.x, p1.y);
  await page.waitForTimeout(100);
  for (const m of mods) await page.keyboard.down(m);
  await page.mouse.down();
  await page.mouse.move((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, {steps: 8});
  await page.mouse.move(p2.x, p2.y, {steps: 8});
  await page.waitForTimeout(150);
  await page.mouse.up();
  for (const m of [...mods].reverse()) await page.keyboard.up(m);
  await page.waitForTimeout(300);
}

/** Drag a closed polygon on the data canvas — the lasso gesture. */
async function lassoCanvas(page: Page, points: Frac[], mods: string[] = ['Shift']): Promise<void> {
  const r = await canvasRect(page);
  const path = points.map((p) => at(r, p));
  await page.mouse.move(path[0].x, path[0].y);
  await page.waitForTimeout(100);
  for (const m of mods) await page.keyboard.down(m);
  await page.mouse.down();
  for (const p of path.slice(1)) {
    await page.mouse.move(p.x, p.y, {steps: 6});
    await page.waitForTimeout(40);
  }
  await page.mouse.move(path[0].x, path[0].y, {steps: 6});
  await page.waitForTimeout(150);
  await page.mouse.up();
  for (const m of [...mods].reverse()) await page.keyboard.up(m);
  await page.waitForTimeout(400);
}

async function clickCanvas(page: Page, p: Frac): Promise<void> {
  const r = await canvasRect(page);
  const pt = at(r, p);
  await page.mouse.click(pt.x, pt.y);
  await page.waitForTimeout(250);
}

/** Wheel-zoom into the middle of the plot and return the viewport it settles on,
 * waited for as a move away from the rectangle the plot showed before. */
async function wheelZoomIn(page: Page, from: Rect, steps = 1): Promise<Rect> {
  const r = await canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  for (let i = 0; i < steps; i++) {
    await page.mouse.wheel(0, -300);
    await page.waitForTimeout(150);
  }
  return await settledViewportAfterChange(page, from);
}

/**
 * Restore the full-data view from the plot's context menu, where Reset View is a
 * top-level entry. This is the only reset route the spec leans on as a utility:
 * the double-click bound to the same command is position sensitive, and it and
 * the keyboard route are exercised as their own steps instead.
 */
async function resetViewFromMenu(page: Page, from?: Rect): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  const item = page.locator('.d4-menu-popup [name="div-Reset-View"]').last();
  await item.waitFor({state: 'visible', timeout: 5000});
  await item.click();
  await page.keyboard.press('Escape');
  for (let i = 0; i < 20; i++) {
    if (await page.locator('.d4-menu-popup:visible').count() === 0) break;
    await page.waitForTimeout(100);
  }
  if (from) await settledViewportAfterChange(page, from);
  else {
    await page.waitForTimeout(250);
    await settledViewport(page);
  }
}

/** Fractions probed for a spot with no marker under the pointer, all kept away
 * from the edges, where a gesture misses the plot area instead of a marker. */
const MARKER_FREE_CANDIDATES: Frac[] = [
  {fx: 0.14, fy: 0.14}, {fx: 0.86, fy: 0.14}, {fx: 0.14, fy: 0.86}, {fx: 0.86, fy: 0.86},
  {fx: 0.22, fy: 0.16}, {fx: 0.78, fy: 0.84}, {fx: 0.16, fy: 0.5}, {fx: 0.84, fy: 0.5},
  {fx: 0.5, fy: 0.14}, {fx: 0.5, fy: 0.86}, {fx: 0.3, fy: 0.85}, {fx: 0.7, fy: 0.15},
];

/** The plot is canvas-drawn, so the table's mouse-over row index is the only way
 * to tell a gap from a marker. */
async function markerFreePoint(page: Page): Promise<{x: number; y: number} | null> {
  const r = await canvasRect(page);
  for (const c of MARKER_FREE_CANDIDATES) {
    const p = at(r, c);
    await page.mouse.move(p.x, p.y);
    if (await settledMouseOverRow(page) === -1) return p;
  }
  return null;
}

/** The row under the pointer, read once two consecutive readings agree — the
 * viewer resolves the hit a frame after the pointer lands. */
async function settledMouseOverRow(page: Page): Promise<number> {
  let prev = -2;
  for (let i = 0; i < 10; i++) {
    await page.waitForTimeout(120);
    const cur = await page.evaluate(() => grok.shell.tv.dataFrame.mouseOverRowIdx as number);
    if (cur === prev) return cur;
    prev = cur;
  }
  return prev;
}

/** Give the canvas keyboard focus with a real click. Every shortcut below acts
 * on the focused plot, and a click on a marker leaves the selection alone. */
async function focusCanvas(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2);
  await page.waitForTimeout(200);
}

/** Wait for one of the viewer's properties to take a value, so a step that
 * toggles it through a shortcut is not timed by a fixed pause. */
async function waitProp(page: Page, name: string, value: any, timeoutMs = 4000): Promise<void> {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    const cur = await page.evaluate((n: string) => {
      const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
      return sp.props[n];
    }, name);
    if (cur === value) return;
    await page.waitForTimeout(120);
  }
}

/** Page coordinates of a range-slider handle. The sliders are hidden until the
 * pointer is over the viewer, so the caller hovers the plot first. */
const sliderHandlePoint = (page: Page, slider: string, handle: string) =>
  page.evaluate(({s, h}: {s: string; h: string}) => {
    const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
      .find((e) => !e.closest('.d4-dialog'))!;
    const el = root.querySelector(`svg[name="${s}"] [name="${h}"]`) as SVGGraphicsElement | null;
    if (!el) return null;
    const b = el.getBoundingClientRect();
    if (b.width === 0 || b.height === 0) return null;
    return {x: b.x + b.width / 2, y: b.y + b.height / 2};
  }, {s: slider, h: handle});

/** The handle's point once the hover has revealed it. */
async function waitSliderHandle(
  page: Page, slider: string, handle: string,
): Promise<{x: number; y: number} | null> {
  for (let i = 0; i < 20; i++) {
    const p = await sliderHandlePoint(page, slider, handle);
    if (p) return p;
    await page.waitForTimeout(150);
  }
  return null;
}

test('Scatter Plot — Point Selection and Viewport Navigation', async ({page}: {page: Page}) => {
  test.setTimeout(1_500_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
  await page.waitForTimeout(500);
  await pickOnViewer(page, 'x', 'WEIGHT');
  await pickOnViewer(page, 'y', 'HEIGHT');
  const fullRowCount = await rowCount(page);
  expect(fullRowCount).toBeGreaterThan(0);
  // Row count the SEX filter leaves, carried to the step that resets that filter.
  let filteredRows = fullRowCount;

  await softStep('Select points by dragging, add a band, deselect, clear', async () => {
    const errBefore = errCount();
    expect(await settledSelection(page)).toBe(0);

    await dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const first = await settledSelectionAfterChange(page, 0);
    expect(first).toBeGreaterThan(0);

    // A band outside the first rectangle, so the rows it covers can only be new
    // ones — a Shift-drag adds to the selection instead of replacing it.
    await dragCanvas(page, {fx: 0.05, fy: 0.05}, {fx: 0.27, fy: 0.95}, ['Shift']);
    const second = await settledSelectionAfterChange(page, first);
    expect(second).toBeGreaterThan(first);

    await dragCanvas(page, {fx: 0.35, fy: 0.35}, {fx: 0.55, fy: 0.55}, ['Control', 'Shift']);
    const deselected = await settledSelectionAfterChange(page, second);
    expect(deselected).toBeLessThan(second);
    expect(deselected).toBeGreaterThan(0);

    // The clearing click works only while the viewer keeps that behaviour on.
    expect((await viewerProps(page)).resetOnBackgroundClick).toBe(true);
    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledSelectionAfterChange(page, deselected)).toBe(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Selection survives a jitter change', async () => {
    const errBefore = errCount();
    await setSliderProp(page, 'prop-jitter-size', 'marker', JITTER_X);
    await setSliderProp(page, 'prop-jitter-size-y', 'marker', JITTER_Y);
    const jittered = await viewerProps(page);
    expect(jittered.jitter).toBe(JITTER_X);
    expect(jittered.jitterY).toBe(JITTER_Y);

    await dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const selected = await settledSelectionAfterChange(page, 0);
    expect(selected).toBeGreaterThan(0);

    await setSliderProp(page, 'prop-jitter-size', 'marker', JITTER_X_CHANGED);
    expect((await viewerProps(page)).jitter).toBe(JITTER_X_CHANGED);
    await page.waitForTimeout(1500);
    // Markers re-drawn at new offsets, but the rows keep their selected state.
    expect(await settledSelection(page)).toBe(selected);

    await dragCanvas(page, {fx: 0.35, fy: 0.35}, {fx: 0.55, fy: 0.55}, ['Control', 'Shift']);
    const narrowed = await settledSelectionAfterChange(page, selected);
    expect(narrowed).toBeLessThan(selected);
    expect(errCount()).toBe(errBefore);

    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledSelectionAfterChange(page, narrowed)).toBe(0);
    await setSliderProp(page, 'prop-jitter-size', 'marker', 0);
    await setSliderProp(page, 'prop-jitter-size-y', 'marker', 0);
    const reverted = await viewerProps(page);
    expect(reverted.jitter).toBe(0);
    expect(reverted.jitterY).toBe(0);
  });

  await softStep('Zoom, pan and reset the viewport', async () => {
    const errBefore = errCount();
    await openSettings(page);
    await revealPropEditor(page, '[name="prop-view-zoom-and-filter"]', 'data');
    if (await propCellText(page, 'prop-view-zoom-and-filter') !== 'no action')
      await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data', 'no action');
    expect(await propCellText(page, 'prop-view-zoom-and-filter')).toBe('no action');

    await resetViewFromMenu(page);
    const baseline = await viewport(page);
    expect(baseline.width).toBeGreaterThan(0);
    expect(baseline.height).toBeGreaterThan(0);

    await dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.6, fy: 0.6}, ['Alt']);
    const zoomed = await settledViewportAfterChange(page, baseline);
    expect(zoomed.width).toBeLessThan(baseline.width);
    expect(zoomed.height).toBeLessThan(baseline.height);
    // With this mode the viewport moves on its own — no row leaves the table.
    expect(await filterCount(page)).toBe(fullRowCount);

    await dragCanvas(page, {fx: 0.5, fy: 0.5}, {fx: 0.3, fy: 0.3});
    const panned = await settledViewportAfterChange(page, zoomed);
    expect(Math.abs(panned.x - zoomed.x)).toBeGreaterThan(zoomed.width * 0.02);
    expect(Math.abs(panned.width - zoomed.width)).toBeLessThan(zoomed.width * 0.02);

    await resetViewFromMenu(page, panned);
    expectSameViewport(await viewport(page), baseline);

    const wheeled = await wheelZoomIn(page, baseline);
    expect(wheeled.width).toBeLessThan(baseline.width);
    expect(wheeled.height).toBeLessThan(baseline.height);
    await resetViewFromMenu(page, wheeled);
    expectSameViewport(await viewport(page), baseline);

    await hoverPlot(page);
    const handle = await waitSliderHandle(page, 'x-slider', 'min-handle');
    expect(handle).not.toBeNull();
    const beforeSlider = await viewport(page);
    await page.mouse.move(handle!.x, handle!.y);
    await page.waitForTimeout(150);
    await page.mouse.down();
    await page.mouse.move(handle!.x + 15, handle!.y, {steps: 6});
    await page.mouse.move(handle!.x + 30, handle!.y, {steps: 6});
    await page.waitForTimeout(150);
    await page.mouse.up();
    const afterSlider = await settledViewportAfterChange(page, beforeSlider);
    // Pushing the X minimum right cuts the low end of the shown range away.
    expect(afterSlider.x).toBeGreaterThan(beforeSlider.x);
    expect(afterSlider.width).toBeLessThan(beforeSlider.width);

    await resetViewFromMenu(page, afterSlider);
    expectSameViewport(await viewport(page), baseline);

    // The double-click resets only clear of a marker — on one it is a point
    // double-click — so the spot is established as marker-free beforehand.
    const wheeledAgain = await wheelZoomIn(page, baseline);
    expect(wheeledAgain.width).toBeLessThan(baseline.width);
    const gap = await markerFreePoint(page);
    expect(gap).not.toBeNull();
    expect(await page.evaluate(() => grok.shell.tv.dataFrame.mouseOverRowIdx as number)).toBe(-1);
    await page.mouse.dblclick(gap!.x, gap!.y);
    expectSameViewport(await settledViewportAfterChange(page, wheeledAgain), baseline);
    expect(errCount()).toBe(errBefore);

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data', 'filter by zoom');
    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
  });

  await softStep('Keyboard selection and view shortcuts', async () => {
    const errBefore = errCount();
    await v.openFilterPanel(page);
    await v.applyCategoricalFilter(page, 'SEX', ['F'], 300);
    const filtered = await settledFilterAfterChange(page, fullRowCount);
    filteredRows = filtered;
    expect(filtered).toBeLessThan(fullRowCount);
    expect(filtered).toBeGreaterThan(0);

    await focusCanvas(page);
    await page.keyboard.press('Control+a');
    // Select-all covers the rows that pass the filter, not the whole table.
    expect(await settledSelectionAfterChange(page, 0)).toBe(filtered);

    await page.keyboard.press('Control+Shift+a');
    expect(await settledSelectionAfterChange(page, filtered)).toBe(0);

    await dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const dragged = await settledSelectionAfterChange(page, 0);
    expect(dragged).toBeGreaterThan(0);
    await page.keyboard.press('Escape');
    expect(await settledSelectionAfterChange(page, dragged)).toBe(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Keyboard selection and view shortcuts — H restores the view', async () => {
    const errBefore = errCount();
    await resetViewFromMenu(page);
    const beforeZoom = await viewport(page);

    const zoomed = await wheelZoomIn(page, beforeZoom);
    expect(zoomed.width).toBeLessThan(beforeZoom.width);

    await focusCanvas(page);
    await page.keyboard.press('h');
    expectSameViewport(await settledViewportAfterChange(page, zoomed), beforeZoom);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Keyboard selection and view shortcuts — the Lasso Tool', async () => {
    const errBefore = errCount();
    await focusCanvas(page);
    await page.keyboard.press('l');
    await waitProp(page, 'lassoTool', true);
    expect((await viewerProps(page)).lasso).toBe(true);

    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledSelection(page)).toBe(0);
    await lassoCanvas(page, [
      {fx: 0.35, fy: 0.35}, {fx: 0.65, fy: 0.32}, {fx: 0.7, fy: 0.62},
      {fx: 0.5, fy: 0.72}, {fx: 0.33, fy: 0.6},
    ]);
    const lassoed = await settledSelectionAfterChange(page, 0);
    expect(lassoed).toBeGreaterThan(0);

    await focusCanvas(page);
    await page.keyboard.press('l');
    await waitProp(page, 'lassoTool', false);
    expect((await viewerProps(page)).lasso).toBe(false);
    expect(errCount()).toBe(errBefore);

    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledSelectionAfterChange(page, lassoed)).toBe(0);
    await page.locator('[name="viewer-Filters"] [name="icon-arrow-rotate-left"]').click();
    expect(await settledFilterAfterChange(page, filteredRows)).toBe(fullRowCount);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

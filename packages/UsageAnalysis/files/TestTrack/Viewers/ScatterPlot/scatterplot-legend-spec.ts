/* ---
realizes: [scatterplot.cp.legend, viewers.scatter-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const RACE_CATEGORIES = ['Asian', 'Black', 'Caucasian', 'Other'];
const SEX_CATEGORIES = ['F', 'M'];
const CONDITIONAL_RANGES = ['634783-634820', '634820-634885'];
const EMPTY_PROBE_COLUMN = 'ZZ_EMPTY_PROBE';
const NO_VALUE_LABEL = '(no value)';

// Ink bounds are margins, never settled readings: each only has to clear
// re-layout noise. The two fractions bound a collapse against the step's own
// baseline; the panel-respect bound is an absolute floor instead, because the
// two gaps it has to cover differ by an order of magnitude.
const SELECTION_INK_MAX_FRACTION = 0.5;
const PANEL_RESPECT_MARGIN_PX = 400;
const EMPTY_SATURATED_MAX_FRACTION = 0.1;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

// NullError / "Stack trace" / cloned-iframe messages are benign ONLY while the
// ribbon project save is running; everywhere else — notably the legend-click
// guards — that same Dart class is the regression signal.
let inProjectSaveWindow = false;
const isBenignError = (text: string) => {
  if (isAmbientError(text)) return true;
  if (inProjectSaveWindow)
    return /Unable to find element in cloned iframe/.test(text) ||
      /Stack trace [A-Za-z]+/.test(text) ||
      /NullError: method not found: '\w+' on null/.test(text);
  return false;
};

/** Rect of the scatter-plot data canvas, resolved on the viewer root that is
 * NOT inside a dialog (the Formula Lines dialog embeds its own preview plot). */
const canvasRect = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

/** The viewer's DOM node attaches before its canvas is laid out, so readiness
 * is a drawn canvas rather than a fixed pause after the add. */
const waitPlotCanvas = (page: Page) => page.waitForFunction(() => {
  return [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .filter((e) => !e.closest('.d4-dialog'))
    .some((e) => {
      const c = e.querySelector('canvas[name="canvas"]');
      return !!c && c.getBoundingClientRect().width > 0 && c.getBoundingClientRect().height > 0;
    });
}, null, {timeout: 20000});

// The two helpers below serve the PROPERTY-PANEL Markers row only. The shared
// pickColumnViaSelectorTrusted covers on-viewer selectors, which it reveals by
// hovering the plot — a panel row has nothing to hover and opens its own popup.
async function waitBackdrop(page: Page, timeout = 5000): Promise<boolean> {
  return await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout}).then(() => true).catch(() => false);
}

/** Poll a condition to a cap. The caller's own assertion, not the return
 * value, is what grades the outcome. */
async function waitFor(page: Page, probe: () => Promise<boolean>, capMs: number): Promise<boolean> {
  const deadline = Date.now() + capMs;
  for (;;) {
    if (await probe()) return true;
    if (Date.now() >= deadline) return false;
    await page.waitForTimeout(100);
  }
}

/** The popup grid is canvas-rendered, so the name goes in with real keys; the
 * first key is separated because the popup focuses its search input a tick
 * later. Enter alone commits the match, and the Marker property carries the
 * commit, so it is what the wait is anchored on. */
async function commitColumn(page: Page, column: string): Promise<void> {
  const text = column.toLowerCase();
  await page.keyboard.press(text[0]);
  await page.waitForTimeout(150);
  if (text.length > 1) await page.keyboard.type(text.slice(1));
  await page.waitForTimeout(200);
  await page.keyboard.press('Enter');
  await waitFor(page, async () => (await viewerColumns(page)).markers === column, 4000);
  await settledLegend(page);
}

/** Pick a column through an on-viewer selector with trusted input only. */
const pickOnViewer = (page: Page, role: string, column: string) =>
  v.pickColumnViaSelectorTrusted(page, {role, columnName: column});

const MARKER_SELECTOR = '[name="prop-markers"] [name="div-column-combobox-markers"]';

const markerSelector = (page: Page) => page.locator(MARKER_SELECTOR);

async function pickMarkerColumn(page: Page, column: string): Promise<void> {
  await v.openViewerGear(page, 'Scatter plot');
  if (!await revealPropEditor(page, MARKER_SELECTOR, 'marker'))
    throw new Error('the Markers column row never became reachable');
  const sel = markerSelector(page);
  await sel.scrollIntoViewIfNeeded();
  await sel.click();
  if (!await waitBackdrop(page)) throw new Error('marker column popup did not open');
  await commitColumn(page, column);
}

/** Clear the Marker column: the popup's first data row is the empty entry, and
 * the grid being canvas-rendered it is reached by coordinate — one row height
 * below the popup's own header. */
async function clearMarkerColumn(page: Page): Promise<void> {
  await v.openViewerGear(page, 'Scatter plot');
  if (!await revealPropEditor(page, MARKER_SELECTOR, 'marker'))
    throw new Error('the Markers column row never became reachable');
  const sel = markerSelector(page);
  await sel.scrollIntoViewIfNeeded();
  await sel.click();
  if (!await waitBackdrop(page)) throw new Error('marker column popup did not open');
  const pt = await page.evaluate(() => {
    const b = document.querySelector('.d4-column-selector-backdrop')!.getBoundingClientRect();
    return {x: b.x + b.width / 2, y: b.y + 30};
  });
  await page.mouse.click(pt.x, pt.y);
  await waitFor(page, async () => (await viewerColumns(page)).markers === '', 4000);
  await settledLegend(page);
}

interface Legend {
  present: boolean;
  all: number;
  coloring: number;
  extra: number;
  labels: string[];
  colorLabels: string[];
  /** Entries carrying a marker glyph: the separate `-extra` entries when Color
   * and Marker sit on different columns, the color entries themselves when the
   * two share a column and the legend is joint. */
  glyphLabels: string[];
  current: string[];
  dimmed: string[];
}

const readLegend = (page: Page): Promise<Legend> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const empty = {present: false, all: 0, coloring: 0, extra: 0, labels: [] as string[],
    colorLabels: [] as string[], glyphLabels: [] as string[], current: [] as string[], dimmed: [] as string[]};
  if (!root) return empty;
  const items = [...root.querySelectorAll('[name="legend"] .d4-legend-item')] as HTMLElement[];
  const txt = (i: Element) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim();
  const colorItems = items.filter((i) => i.classList.contains('d4-legend-item-coloring'));
  return {
    present: !!root.querySelector('[name="legend"]'),
    all: items.length,
    coloring: colorItems.length,
    extra: items.filter((i) => i.classList.contains('d4-legend-item-extra')).length,
    labels: items.map(txt),
    colorLabels: colorItems.map(txt),
    glyphLabels: items.filter((i) => i.querySelector('i[name="legend-icon-color-picker"]')).map(txt),
    current: items.filter((i) => i.classList.contains('d4-legend-item-current')).map(txt),
    // Dimming comes from the legend's own stylesheet, not from an inline style.
    dimmed: items.filter((i) => parseFloat(getComputedStyle(i).opacity) < 0.9).map(txt),
  };
});

/** The legend is rebuilt a tick after the property that drives it, so it is
 * read only once two consecutive reads agree. */
async function settledLegend(page: Page): Promise<Legend> {
  let prev = await readLegend(page);
  for (let i = 0; i < 30; i++) {
    await page.waitForTimeout(100);
    const cur = await readLegend(page);
    if (JSON.stringify(cur) === JSON.stringify(prev)) return cur;
    prev = cur;
  }
  return prev;
}

interface Ink {
  nonWhite: number;
  saturated: number;
  /** Near-neutral pixels painted pale — the shade empty-value markers use. */
  pale: number;
}

/** Ink on the plot's data canvas. Colour buckets are NOT one-to-one with
 * categories (overlap shading smears a category across neighbouring hues), so
 * the split is saturated-versus-pale rather than per-category. */
const readInk = (page: Page): Promise<Ink> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const canvas = root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
  const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
  const d = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
  let nonWhite = 0; let saturated = 0; let pale = 0;
  for (let i = 0; i < d.length; i += 4) {
    if (d[i + 3] === 0) continue;
    const r = d[i]; const g = d[i + 1]; const b = d[i + 2];
    if (r >= 250 && g >= 250 && b >= 250) continue;
    nonWhite++;
    const spread = Math.max(r, g, b) - Math.min(r, g, b);
    if (spread >= 12) saturated++;
    else if (Math.max(r, g, b) >= 200) pale++;
  }
  return {nonWhite, saturated, pale};
});

/** Park the pointer clear of the plot: a hovered point paints its mouse-over
 * indicator into the frame and spoils the measurement. */
async function parkPointer(page: Page): Promise<void> {
  await page.mouse.move(4, 4);
  await page.waitForTimeout(250);
}

const sameInk = (a: Ink, b: Ink) =>
  a.nonWhite === b.nonWhite && a.saturated === b.saturated && a.pale === b.pale;

/** The canvas does not drift at rest, so two consecutive identical reads are
 * the settle gate and no tolerance is applied. */
async function settledInk(page: Page): Promise<Ink> {
  await parkPointer(page);
  let prev = await readInk(page);
  for (let i = 0; i < 40; i++) {
    await page.waitForTimeout(250);
    const cur = await readInk(page);
    if (sameInk(cur, prev)) return cur;
    prev = cur;
  }
  return prev;
}

/** Ink after an action that must repaint the plot. Identical reads alone would
 * settle on the frame the action has not reached yet, so the gate is the canvas
 * first leaving the frame the caller knows and only then holding still. A
 * canvas that never repaints is returned as it is, so the caller's assertion
 * fails on the real frame. */
async function settledInkAfterChange(page: Page, from: Ink): Promise<Ink> {
  await parkPointer(page);
  let prev = await readInk(page);
  for (let i = 0; i < 60; i++) {
    await page.waitForTimeout(250);
    const cur = await readInk(page);
    if (!sameInk(cur, from) && sameInk(cur, prev)) return cur;
    prev = cur;
  }
  return prev;
}

/** Which categories the legend marks as selected — the state a click flips. */
const legendSelection = async (page: Page) => (await readLegend(page)).current.join('|');

/** Real click on a legend entry, which is the only input this surface takes.
 * The wait is anchored on the selection actually flipping. */
async function clickLegendEntry(page: Page, label: string): Promise<void> {
  const before = await legendSelection(page);
  const pt = await legendPoint(page, label);
  expect(pt, `legend entry ${label}`).not.toBeNull();
  await page.mouse.click(pt!.x, pt!.y);
  await waitFor(page, async () => await legendSelection(page) !== before, 4000);
}

/** Release whatever legend category is selected. Tolerant by design: it runs in
 * teardowns, where a missing entry must not mask the original failure. */
async function clearLegendSelection(page: Page): Promise<void> {
  if (!await legendFiltering(page)) return;
  const current = (await readLegend(page)).current;
  if (!current.length) return;
  const pt = await legendPoint(page, current[0]);
  if (!pt) return;
  await page.mouse.click(pt.x, pt.y);
  await waitFor(page, async () => !await legendFiltering(page), 4000);
}

/** Whether the legend is in its filtering state — the class the host gains once
 * a category is selected, and the DOM support for the ink measurement. */
const legendFiltering = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const host = root?.querySelector('[name="legend"]');
  return !!host?.classList.contains('d4-legend-filtering');
});

/** Page coordinates of a legend entry by its text. The legend only responds to
 * trusted input, so callers drive it with `page.mouse.click` on this point. */
const legendPoint = (page: Page, label: string) => page.evaluate((l: string) => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const item = [...(root?.querySelectorAll('[name="legend"] .d4-legend-item') ?? [])]
    .find((i) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim() === l);
  if (!item) return null;
  const el = item.querySelector('.d4-legend-value') ?? item;
  const b = el.getBoundingClientRect();
  return {x: b.x + b.width / 2, y: b.y + b.height / 2};
}, label);

const filterCount = (page: Page) =>
  page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount as number);

/** Filtered row count where a NEW value is expected. The count lags the action
 * that changes it, so the settle gate is it LEAVING the value the caller knows
 * and only then holding still — agreeing samples on the previous value are not
 * a settled reading. A count that never moves is returned as it is, so the
 * caller's assertion fails on the real value. */
async function settledFilterCountAfterChange(page: Page, from: number): Promise<number> {
  let prev = -1;
  let stable = 0;
  for (let i = 0; i < 70; i++) {
    const c = await filterCount(page);
    if (c === prev) {
      stable++;
      if (c !== from && stable >= 4) return c;
    } else {
      stable = 0;
      prev = c;
    }
    await page.waitForTimeout(250);
  }
  return prev;
}

/** Filtered row count where the expectation is that it does NOT move — the
 * legend-click guards, whose whole point is that the table filter stays put.
 * "Did not change" cannot be told apart from "has not changed yet", so
 * sampling starts only after a time floor that outlasts a filter update. */
async function settledFilterCountUnchanged(page: Page, floorMs = 2000): Promise<number> {
  await page.waitForTimeout(floorMs);
  let last = -1;
  let stable = 0;
  for (let i = 0; i < 25; i++) {
    const c = await filterCount(page);
    if (c === last) {
      stable++;
      if (stable >= 3) return c;
    } else {
      stable = 0;
      last = c;
    }
    await page.waitForTimeout(400);
  }
  return last;
}

const viewerColumns = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {color: sp.props.colorColumnName, markers: sp.props.markersColumnName};
});

/** Count of legend hosts inside the viewer root: Legend Visibility Never
 * removes the host outright, so this is the presence signal. */
const legendHostCount = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  return root ? root.querySelectorAll('[name="legend"]').length : 0;
});

const legendVisibility = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp.props.legendVisibility as string;
});

const VISIBILITY_MENU = ['div-Properties...', 'div-Properties...---Legend',
  'div-Properties...---Legend---Legend-Visibility'];

/** Wait for the viewer's properties in the Context Panel, clicking the gear only
 * when they are not there yet. The gear does not toggle, so the guard saves a
 * redundant click rather than preventing a close. */
async function settingsPanelReady(page: Page): Promise<boolean> {
  for (let i = 0; i < 4; i++) {
    if (await page.evaluate(() => !!document.querySelector('[name="prop-category-legend"]'))) return true;
    await v.openViewerGear(page, 'Scatter plot');
    await page.waitForTimeout(1000);
  }
  return false;
}

/** Make a property editor reachable. A row inside a collapsed category keeps its
 * DOM node with an empty box, so readiness is measured on the editor's own
 * rectangle; the category header is a toggle, so an attempt is simply retried. */
async function revealPropEditor(page: Page, editorSelector: string, category: string): Promise<boolean> {
  for (let i = 0; i < 8; i++) {
    const ready = await page.evaluate((sel: string) => {
      const el = document.querySelector(sel) as HTMLElement | null;
      if (!el || !el.offsetParent) return false;
      const b = el.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    }, editorSelector);
    if (ready) return true;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();
    await page.waitForTimeout(800);
  }
  return false;
}

/** The plot's context menu only opens on a real right-click. */
async function openPlotContextMenu(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await page.waitForTimeout(500);
}

/**
 * Set Legend Visibility through the real UI. Both routes — the settings-panel
 * select and the context-menu mirror — enumerate their options rather than
 * assume them, and both end in a read-back of the property, so a route that
 * silently does nothing fails here instead of letting the presence assert pass
 * vacuously.
 */
async function setLegendVisibility(page: Page, value: string): Promise<void> {
  let driven = false;
  const editor = '[name="prop-view-legend-visibility"]';
  if (await settingsPanelReady(page) && await revealPropEditor(page, editor, 'legend')) {
    const cell = page.locator(editor);
    await cell.scrollIntoViewIfNeeded();
    await cell.click();
    const select = page.locator('[name="prop-legend-visibility"] select.property-grid-item-editor-spinner');
    await select.waitFor({timeout: 3000}).catch(() => {});
    if (await select.count() > 0) {
      const options = (await select.locator('option').allTextContents()).map((o) => o.trim());
      expect(options).toContain(value);
      await select.selectOption(value);
      await waitFor(page, async () => await legendVisibility(page) === value, 4000);
      driven = true;
    }
  }
  if (!driven) {
    await openPlotContextMenu(page);
    for (const g of VISIBILITY_MENU) {
      await page.locator(`.d4-menu-popup [name="${g}"]`).last().hover();
      await page.waitForTimeout(700);
    }
    const leaf = page.locator(`.d4-menu-popup [name="${VISIBILITY_MENU[2]}---${value}"]`);
    expect(await leaf.count()).toBeGreaterThan(0);
    await leaf.last().click();
    await waitFor(page, async () => await legendVisibility(page) === value, 4000);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(500);
  }
  expect(await legendVisibility(page)).toBe(value);
}

/** Grid header cells are canvas-drawn, so the right-click is positioned through
 * the grid's own column geometry. */
async function openColumnHeaderMenu(page: Page, column: string): Promise<void> {
  const pt = await page.evaluate((name: string) => {
    const grid = grok.shell.tv.grid;
    const r = grid.root.getBoundingClientRect();
    const gc = grid.columns.byName(name);
    return {x: r.x + gc.left + gc.width / 2, y: r.y + (grid.colHeaderH || 20) / 2};
  }, column);
  await page.mouse.click(pt.x, pt.y, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await page.waitForTimeout(400);
}

/** Submenu leaves stay hidden until their parent group is hovered. */
async function clickHeaderMenuLeaf(page: Page, group: string, leaf: string): Promise<void> {
  await page.locator(`.d4-menu-popup [name="${group}"]`).last().hover();
  await page.waitForTimeout(600);
  await page.locator(`.d4-menu-popup [name="${leaf}"]`).last().click();
  await page.waitForTimeout(1200);
}

const conditionalRanges = (page: Page, column: string) => page.evaluate((c: string) => {
  const col = grok.shell.tv.dataFrame.col(c);
  const raw = col.tags['.color-coding-conditional'];
  return {type: col.tags['.color-coding-type'] ?? null, ranges: raw ? Object.keys(JSON.parse(raw)) : []};
}, column);

test('Scatter Plot — Legend Lifecycle, Filter Interplay, Persistence', async ({page}: {page: Page}) => {
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
  await waitPlotCanvas(page);
  await pickOnViewer(page, 'x', 'WEIGHT');
  await pickOnViewer(page, 'y', 'HEIGHT');
  const fullRowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);
  expect(fullRowCount).toBeGreaterThan(0);

  await softStep('Color legend and the joint Color plus Marker legend', async () => {
    const errBefore = errCount();

    await pickOnViewer(page, 'color', 'RACE');
    const colorOnly = await readLegend(page);
    expect(colorOnly.present).toBe(true);
    expect(colorOnly.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(colorOnly.all).toBe(RACE_CATEGORIES.length);

    await pickMarkerColumn(page, 'SEX');
    const joint = await readLegend(page);
    expect(joint.coloring).toBe(RACE_CATEGORIES.length);
    expect(joint.extra).toBe(SEX_CATEGORIES.length);
    expect(joint.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(joint.all).toBe(RACE_CATEGORIES.length + SEX_CATEGORIES.length);

    // A numeric Color column must not take the marker glyph entries with it.
    await pickOnViewer(page, 'color', 'AGE');
    const numericColor = await readLegend(page);
    expect(numericColor.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(numericColor.extra).toBe(SEX_CATEGORIES.length);

    await pickOnViewer(page, 'color', 'RACE');
    const back = await readLegend(page);
    expect(back.coloring).toBe(RACE_CATEGORIES.length);
    expect(back.extra).toBe(SEX_CATEGORIES.length);
    expect(back.all).toBe(RACE_CATEGORIES.length + SEX_CATEGORIES.length);
    expect(back.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(new Set(back.labels).size).toBe(back.labels.length);

    // Color and Marker on one column form a joint legend: the glyphs sit on the
    // color entries, so both halves cover the same category set.
    await pickMarkerColumn(page, 'RACE');
    const same = await readLegend(page);
    expect(same.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(same.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(same.all).toBe(RACE_CATEGORIES.length);

    await pickMarkerColumn(page, 'SEX');
    expect((await viewerColumns(page)).markers).toBe('SEX');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Clearing the Marker column removes its glyph entries', async () => {
    const errBefore = errCount();
    await pickOnViewer(page, 'color', 'SEX');
    await pickMarkerColumn(page, 'SEX');
    const before = await readLegend(page);
    expect(before.colorLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(before.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());

    await clearMarkerColumn(page);
    const after = await readLegend(page);
    // Clearing a column sets the property to the empty string, not to null.
    expect((await viewerColumns(page)).markers).toBe('');
    expect(after.glyphLabels).toEqual([]);
    expect(after.colorLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(after.present).toBe(true);

    await pickOnViewer(page, 'color', 'RACE');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Filtered-out categories are absent from the marker legend', async () => {
    const errBefore = errCount();
    await v.openFilterPanel(page);

    // Marker column set AFTER the filter.
    await v.applyCategoricalFilter(page, 'RACE', ['Asian', 'Caucasian'], 600);
    const narrowed = await settledFilterCountAfterChange(page, fullRowCount);
    expect(narrowed).toBeLessThan(fullRowCount);
    await pickMarkerColumn(page, 'RACE');
    const afterFilter = await readLegend(page);
    expect(afterFilter.glyphLabels.sort()).toEqual(['Asian', 'Caucasian']);

    // Same subset, opposite ordering: Marker column set BEFORE the filter.
    await v.resetFilters(page);
    await settledFilterCountAfterChange(page, narrowed);
    await clearMarkerColumn(page);
    await pickMarkerColumn(page, 'RACE');
    const unfiltered = await readLegend(page);
    expect(unfiltered.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    await v.applyCategoricalFilter(page, 'RACE', ['Asian', 'Caucasian'], 600);
    await settledFilterCountAfterChange(page, fullRowCount);
    await settledLegend(page);
    const beforeFilter = await readLegend(page);
    expect(beforeFilter.glyphLabels.sort()).toEqual(['Asian', 'Caucasian']);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Filtered-out categories are absent from the marker legend — filter cleared', async () => {
    const before = await filterCount(page);
    await v.resetFilters(page);
    const restored = await settledFilterCountAfterChange(page, before);
    expect(restored).toBe(fullRowCount);
    const legend = await settledLegend(page);
    expect(legend.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
  });

  await softStep('Clicking a legend category hides the other categories on the plot', async () => {
    const errBefore = errCount();
    await pickOnViewer(page, 'color', 'RACE');
    await pickMarkerColumn(page, 'SEX');
    try {
      const baseline = await settledInk(page);
      expect(baseline.nonWhite).toBeGreaterThan(0);
      const baseCount = await settledFilterCountUnchanged(page, 500);
      expect(baseCount).toBe(fullRowCount);

      // The legend filters what the viewer DRAWS, so the table's own filter has
      // to stay where it was.
      await clickLegendEntry(page, 'Asian');
      const legendOnly = await settledInkAfterChange(page, baseline);
      expect(legendOnly.nonWhite)
        .toBeLessThan(Math.round(baseline.nonWhite * SELECTION_INK_MAX_FRACTION));
      expect(await legendFiltering(page)).toBe(true);
      const marked = await readLegend(page);
      expect(marked.current).toContain('Asian');
      expect(marked.dimmed.length).toBeGreaterThan(0);
      expect(await settledFilterCountUnchanged(page)).toBe(baseCount);

      // GROK-14940: a Filter-Panel filter on top must not bring back what the
      // legend removed. The reference is the legend-only plot, NOT the
      // panel-only one: dropping the other categories also drops marker
      // overlap, so the selected category's own ink grows.
      await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
      const narrowed = await settledFilterCountAfterChange(page, fullRowCount);
      expect(narrowed).toBeLessThan(fullRowCount);
      const legendThenPanel = await settledInkAfterChange(page, legendOnly);
      expect(legendThenPanel.nonWhite)
        .toBeLessThan(legendOnly.nonWhite - PANEL_RESPECT_MARGIN_PX);
      expect((await readLegend(page)).labels).not.toContain('M');

      // The same two operations in the opposite order reach the same frame.
      await clearLegendSelection(page);
      await v.resetFilters(page);
      await settledFilterCountAfterChange(page, narrowed);
      await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
      await settledFilterCountAfterChange(page, fullRowCount);
      const panelOnly = await settledInk(page);
      await clickLegendEntry(page, 'Asian');
      const panelThenLegend = await settledInkAfterChange(page, panelOnly);
      expect(panelThenLegend).toEqual(legendThenPanel);
      expect(errCount()).toBe(errBefore);
    } finally {
      await clearLegendSelection(page);
      await v.resetFilters(page);
      await settledFilterCountUnchanged(page, 500);
    }
  });

  await softStep('Clicking a legend category hides the other categories on the plot — second click restores the plot exactly',
    async () => {
      const errBefore = errCount();
      const baseline = await settledInk(page);
      try {
        await clickLegendEntry(page, 'Asian');
        expect(await legendFiltering(page)).toBe(true);
        const selected = await settledInkAfterChange(page, baseline);
        await clickLegendEntry(page, 'Asian');

        // The canvas does not drift at rest, so the restore is exact.
        expect(await settledInkAfterChange(page, selected)).toEqual(baseline);
        const legend = await readLegend(page);
        expect(legend.current).toEqual([]);
        expect(legend.dimmed).toEqual([]);
        expect(await legendFiltering(page)).toBe(false);
        expect(errCount()).toBe(errBefore);
      } finally {
        await clearLegendSelection(page);
      }
    });

  await softStep('Clicking a legend category hides the other categories on the plot — the entry for empty values',
    async () => {
      const errBefore = errCount();
      const baseCount = await settledFilterCountUnchanged(page, 500);
      try {
        // The demog categorical columns carry no empty values, so the column
        // that produces the empty-value legend entry is added here. Its markers
        // are painted pale rather than in a category colour.
        await page.evaluate((name: string) => {
          const df = grok.shell.tv.dataFrame;
          const col = df.columns.addNewString(name);
          for (let i = 0; i < df.rowCount; i++)
            col.set(i, i % 3 === 0 ? null : (i % 3 === 1 ? 'alpha' : 'beta'));
        }, EMPTY_PROBE_COLUMN);
        await page.waitForTimeout(1200);
        await pickOnViewer(page, 'color', EMPTY_PROBE_COLUMN);
        expect((await settledLegend(page)).colorLabels).toContain(NO_VALUE_LABEL);

        const baseline = await settledInk(page);
        expect(baseline.saturated).toBeGreaterThan(0);

        await clickLegendEntry(page, NO_VALUE_LABEL);
        const selected = await settledInkAfterChange(page, baseline);
        expect(selected.saturated)
          .toBeLessThan(Math.round(baseline.saturated * EMPTY_SATURATED_MAX_FRACTION));
        expect(selected.pale).toBeGreaterThan(baseline.pale);
        expect(await legendFiltering(page)).toBe(true);
        expect((await readLegend(page)).current).toContain(NO_VALUE_LABEL);
        expect(await settledFilterCountUnchanged(page)).toBe(baseCount);

        // GROK-20228 respects the panel the same way a valued category does.
        await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
        expect(await settledFilterCountAfterChange(page, baseCount)).toBeLessThan(baseCount);
        const withPanel = await settledInkAfterChange(page, selected);
        expect(withPanel.nonWhite).toBeLessThan(selected.nonWhite - PANEL_RESPECT_MARGIN_PX);

        await clearLegendSelection(page);
        const filtered = await filterCount(page);
        await v.resetFilters(page);
        await settledFilterCountAfterChange(page, filtered);
        expect(await settledInkAfterChange(page, withPanel)).toEqual(baseline);
        expect(errCount()).toBe(errBefore);
      } finally {
        await clearLegendSelection(page);
        await v.resetFilters(page);
        await page.evaluate((name: string) => {
          const df = grok.shell.tv.dataFrame;
          if (df.col(name)) df.columns.remove(name);
        }, EMPTY_PROBE_COLUMN);
        await pickOnViewer(page, 'color', 'RACE');
        await pickMarkerColumn(page, 'SEX');
      }
      const peak = await viewerColumns(page);
      expect(peak.color).toBe('RACE');
      expect(peak.markers).toBe('SEX');
    });

  await softStep('A column with conditional color coding still produces a legend', async () => {
    const errBefore = errCount();
    try {
      await page.evaluate(async (path: string) => {
        const df = await grok.dapi.files.readCsv(path);
        grok.shell.addTableView(df);
        await new Promise((resolve) => {
          const s = df.onSemanticTypeDetected.subscribe(() => { s.unsubscribe(); resolve(null); });
          setTimeout(resolve, 5000);
        });
      }, spgiPath);
      await page.waitForFunction(() =>
        !!grok.shell.tv?.dataFrame?.col('CAST Idea ID'), null, {timeout: 30000});

      // Both ranges must occur in the data — bounds outside it collapse into a
      // single no-value entry.
      await openColumnHeaderMenu(page, 'CAST Idea ID');
      await clickHeaderMenuLeaf(page, 'div-Color-Coding', 'div-Color-Coding---Conditional');
      await openColumnHeaderMenu(page, 'CAST Idea ID');
      await clickHeaderMenuLeaf(page, 'div-Color-Coding', 'div-Color-Coding---Edit...');

      const dialog = page.locator('[name^="dialog-Color-coding"]');
      await dialog.waitFor({timeout: 10000});
      const rangeInputs = dialog.locator('input.ui-input-editor');
      for (let guard = 0; guard < 10 && await rangeInputs.count() > CONDITIONAL_RANGES.length; guard++) {
        await dialog.locator('[name="button-Remove-row"]').last().click();
        await page.waitForTimeout(600);
      }
      expect(await rangeInputs.count()).toBe(CONDITIONAL_RANGES.length);
      for (let i = 0; i < CONDITIONAL_RANGES.length; i++) {
        await rangeInputs.nth(i).click();
        await page.keyboard.press('Control+A');
        await page.keyboard.type(CONDITIONAL_RANGES[i]);
        await page.keyboard.press('Enter');
        await page.waitForTimeout(800);
      }
      await dialog.locator('[name="button-CLOSE"]').click();
      await page.waitForTimeout(1200);

      const coding = await conditionalRanges(page, 'CAST Idea ID');
      expect(coding.type).toBe('Conditional');
      expect(coding.ranges).toEqual(CONDITIONAL_RANGES);

      await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
      await waitPlotCanvas(page);
      const axes = await page.evaluate(() => {
        const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        return {x: sp.props.xColumnName, y: sp.props.yColumnName};
      });
      expect(axes.x).toBe('CAST Idea ID');
      expect(axes.y).toBe('Idea ID');

      await pickOnViewer(page, 'color', 'CAST Idea ID');
      const legend = await readLegend(page);
      expect(legend.present).toBe(true);
      expect(legend.colorLabels).toEqual(CONDITIONAL_RANGES);
      expect(errCount()).toBe(errBefore);
    } finally {
      // The demog view is identified by its columns rather than by a table
      // name — the loader does not name a table after its file.
      await page.evaluate(async () => {
        grok.shell.v.close();
        await new Promise((r) => setTimeout(r, 1500));
        for (const view of grok.shell.tableViews) {
          if (view.dataFrame.columns.contains('RACE')) {
            grok.shell.v = view;
            break;
          }
        }
        await new Promise((r) => setTimeout(r, 1500));
      });
    }
    const peak = await viewerColumns(page);
    expect(peak.color).toBe('RACE');
    expect(peak.markers).toBe('SEX');
  });

  await softStep('Legend visibility hides and restores the legend', async () => {
    const errBefore = errCount();
    const initial = await legendVisibility(page);
    const before = await readLegend(page);
    expect(before.present).toBe(true);
    expect(await legendHostCount(page)).toBe(1);

    try {
      await setLegendVisibility(page, 'Never');
      expect(await legendHostCount(page)).toBe(0);
      expect((await readLegend(page)).present).toBe(false);

      await setLegendVisibility(page, 'Always');
      expect(await legendHostCount(page)).toBe(1);
      const shown = await readLegend(page);
      expect(shown.present).toBe(true);
      expect(shown.all).toBe(before.all);
      expect(shown.colorLabels.sort()).toEqual([...before.colorLabels].sort());
      expect(shown.glyphLabels.sort()).toEqual([...before.glyphLabels].sort());

      await setLegendVisibility(page, initial);
      const restored = await readLegend(page);
      expect(restored.present).toBe(true);
      expect(restored.all).toBe(before.all);
      expect(restored.colorLabels.sort()).toEqual([...before.colorLabels].sort());
      expect(errCount()).toBe(errBefore);
    } finally {
      // Cleanup only: the value is forced back if the UI route left it anywhere
      // else, because the persistence tail must start from the peak setup.
      await page.evaluate((val: string) => {
        const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        if (sp && sp.props.legendVisibility !== val) sp.props.legendVisibility = val;
      }, initial);
      await page.waitForTimeout(1200);
    }
  });

  await softStep('Layout and project persistence at peak configuration — layout round-trip', async () => {
    const peak = await viewerColumns(page);
    expect(peak.color).toBe('RACE');
    expect(peak.markers).toBe('SEX');
    const peakLegend = await readLegend(page);
    expect(peakLegend.coloring).toBe(RACE_CATEGORIES.length);
    expect(peakLegend.extra).toBe(SEX_CATEGORIES.length);

    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1500));
      return layout.id as string;
    });
    try {
      const result = await page.evaluate(async (id: string) => {
        const tv = grok.shell.tv;
        tv.addViewer('Histogram');
        await new Promise((r) => setTimeout(r, 1200));
        const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        sp.props.colorColumnName = '';
        await new Promise((r) => setTimeout(r, 1000));
        const clearedColor = sp.props.colorColumnName;
        tv.loadLayout(await grok.dapi.layouts.find(id));
        await new Promise((r) => setTimeout(r, 4000));
        const restored = tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        return {
          clearedColor,
          hasScatter: tv.viewers.some((x: any) => x.type === 'Scatter plot'),
          hasHistogram: tv.viewers.some((x: any) => x.type === 'Histogram'),
          color: restored?.props.colorColumnName, markers: restored?.props.markersColumnName,
        };
      }, layoutId);

      // The perturbation really took effect before the re-apply.
      expect(result.clearedColor).toBe('');
      // The restored viewer set is the SAVED one, not the perturbed one.
      expect(result.hasScatter).toBe(true);
      expect(result.hasHistogram).toBe(false);
      expect(result.color).toBe('RACE');
      expect(result.markers).toBe('SEX');

      const legend = await readLegend(page);
      expect(legend.present).toBe(true);
      expect(legend.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
      expect(legend.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    } finally {
      await page.evaluate(async (id: string) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved) await grok.dapi.layouts.delete(saved);
        } catch (_) { /* best effort */ }
      }, layoutId);
    }
  });

  await softStep('Layout and project persistence at peak configuration — project save / Close All / reopen',
    async () => {
      const projName = 'zz-scatterplot-legend-probe-' + Date.now();
      let projectId: string | null = null;
      inProjectSaveWindow = true;
      try {
        const saved = await saveProjectViaUI(page, projName);
        projectId = saved.projectId;
        expect(projectId).toBeTruthy();

        const result = await page.evaluate(async (id: string) => {
          grok.shell.closeAll();
          await new Promise((r) => setTimeout(r, 1500));
          const full = await grok.dapi.projects.find(id);
          await full.open();
          let sp: any = null;
          for (let t = 0; t < 20 && !sp; t++) {
            await new Promise((r) => setTimeout(r, 1000));
            for (const view of grok.shell.tableViews)
              for (const vw of view.viewers)
                if (vw.type === 'Scatter plot') sp = vw;
          }
          const items = sp ? [...sp.root.querySelectorAll('[name="legend"] .d4-legend-item')] : [];
          const txt = (i: Element) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim();
          return {
            found: !!sp,
            color: sp?.props.colorColumnName, markers: sp?.props.markersColumnName,
            legendPresent: !!sp?.root.querySelector('[name="legend"]'),
            colorLabels: items.filter((i) => i.classList.contains('d4-legend-item-coloring')).map(txt),
            glyphLabels: items.filter((i) => i.querySelector('i[name="legend-icon-color-picker"]')).map(txt),
          };
        }, projectId);

        expect(result.found).toBe(true);
        expect(result.color).toBe('RACE');
        expect(result.markers).toBe('SEX');
        expect(result.legendPresent).toBe(true);
        expect(result.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
        expect(result.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
      } finally {
        inProjectSaveWindow = false;
        if (projectId) await deleteProjectWithCleanup(page, {projectId});
      }
    });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

/* ---
realizes: [scatterplot.cp.labels-tooltip, viewers.scatter-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const SETUP_X = 'WEIGHT';
const SETUP_Y = 'HEIGHT';
// The hover aim is computed from this row's X and Y values, but several demog
// rows share almost the same (WEIGHT, HEIGHT), so the row actually under the
// pointer is read back from the table and graded instead of assumed.
const REFERENCE_ROW = 10;
const CUSTOM_TOOLTIP_COLUMNS = ['AGE', 'SEX'];
const LABEL_COLUMN = 'AGE';
const SPGI_X = 'Whole blood assay 2 Date';
const SPGI_Y = 'Stereo Category';
const SPGI_LABEL_COLUMN = 'Id';

// Canvas readings scale with viewer size and dpr, so only these bounds are
// asserted, never a reading. Their ordering keeps the asserts honest:
// INK_SETTLE < OVERLAY_RESTORE_TOLERANCE < LABEL_INK_DELTA.
const INK_SETTLE = 60;
const LABEL_INK_DELTA = 500;
const OVERLAY_RESTORE_TOLERANCE = 150;
// How far, as a fraction of the viewport rectangle, the row the viewer reports
// may sit from the reference row's coordinates and still count as the aim.
const MARKER_AIM_TOLERANCE = 0.05;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

// No project save here, so no Dart NullError or "Stack trace" window is
// expected: nothing beyond the shared server's ambient noise is whitelisted,
// and the Infinity.ceil() class reaches the guards unconditionally.
const isBenignError = (text: string) => isAmbientError(text);

const CRASH_SIGNATURE = /Infinity\.ceil/i;

interface Rect {x: number; y: number; width: number; height: number}

/** Rect of the scatter-plot data canvas, resolved on the viewer root that is NOT
 * inside a dialog — a dialog can embed its own preview plot with a canvas of the
 * same name. */
const canvasRect = (page: Page): Promise<Rect> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

/** Pick a column through an on-viewer selector with trusted input only — the
 * shared helper also covers the long, space-bearing names of the chemical
 * dataset, and raises when the property does not end up holding the column. */
const pickOnViewer = (page: Page, role: string, column: string) =>
  v.pickColumnViaSelectorTrusted(page, {role, columnName: column});

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

/** Make a property editor reachable. A row inside a collapsed property-grid
 * category keeps its DOM node with an empty box, so the check is on the editor's
 * own rectangle; the category header is a toggle, so an attempt that collapsed
 * it is simply retried. */
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

/** Set a choice property: the value cell turns into a select once clicked. */
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
  const prop = rowName.replace(/^prop-/, '').replace(/-([a-z])/g, (_, c: string) => c.toUpperCase());
  for (let i = 0; i < 25; i++) {
    if (await readProp(page, prop) === value) break;
    await page.waitForTimeout(100);
  }
  await page.waitForTimeout(250);
}

const readProp = (page: Page, name: string) => page.evaluate((n: string) => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp.props[n];
}, name);

/**
 * Assign one of the viewer's list-valued properties — the label column list and
 * the custom tooltip column list. Their UI editor is the Select columns dialog,
 * whose checkbox list is a canvas-rendered grid: no DOM checkbox, All / None
 * ignore the search box, and no coordinate toggles a column. The column choice
 * is configuration here, not the graded behaviour — what is graded is the
 * tooltip text and the label overlay the viewer then produces.
 */
async function setListProp(page: Page, name: string, value: string | string[]): Promise<void> {
  await page.evaluate(({n, val}: {n: string; val: any}) => {
    const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    sp.props[n] = val;
  }, {n: name, val: value});
  await page.waitForTimeout(400);
}

/** Non-white opaque pixels on one of the plot's two canvases. `canvas` carries
 * the markers, axes and grid lines; `overlay` carries the marker labels. -1
 * signals a canvas that cannot be read, which the caller treats as a fault. */
const canvasInk = (page: Page, layer: 'canvas' | 'overlay') => page.evaluate((name: string) => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const c = root?.querySelector(`canvas[name="${name}"]`) as HTMLCanvasElement | null;
  const ctx = c?.getContext('2d');
  if (!c || !ctx) return -1;
  let data: Uint8ClampedArray;
  try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { return -1; }
  let n = 0;
  for (let k = 0; k < data.length; k += 16)
    if (data[k + 3] !== 0 && !(data[k] >= 250 && data[k + 1] >= 250 && data[k + 2] >= 250)) n++;
  return n;
}, layer);

/** Move the pointer clear of every viewer. The plot repaints the marker under
 * the cursor, and this file parks the pointer on markers on purpose to read the
 * tooltip, so a measurement taken without parking can catch that highlight. */
async function parkPointer(page: Page): Promise<void> {
  await page.mouse.move(4, 4);
  await page.waitForTimeout(250);
}

/** Read a canvas until two consecutive readings agree within the settle drift —
 * a repaint finishes after the change that triggered it has been applied. */
async function settledInk(page: Page, layer: 'canvas' | 'overlay'): Promise<number> {
  await parkPointer(page);
  let prev = await canvasInk(page, layer);
  let cur = prev;
  for (let i = 0; i < 12; i++) {
    await page.waitForTimeout(250);
    cur = await canvasInk(page, layer);
    if (cur >= 0 && Math.abs(cur - prev) <= INK_SETTLE) break;
    prev = cur;
  }
  return cur;
}

interface TooltipEntry {name: string; value: string}

/** The tooltip's contents as name/value pairs. The tooltip is a real DOM table:
 * each row carries the column name in its first cell and the formatted value in
 * a value cell. In the custom mode the name cell is absent, so the name comes
 * back empty and the values still read out in order. */
const tooltipEntries = (page: Page): Promise<TooltipEntry[]> => page.evaluate(() => {
  const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
  if (!tt) return [];
  return [...tt.querySelectorAll('tr')].map((tr) => {
    const value = (tr.querySelector('.d4-tooltip-text-value-cell')?.textContent ?? '').trim();
    const cells = [...tr.querySelectorAll('td')];
    const first = cells.length > 0 ? (cells[0].textContent ?? '').trim() : '';
    return {name: first === value ? '' : first, value};
  }).filter((e) => e.value.length > 0 || e.name.length > 0);
});

/** The tooltip's visible text, empty when no tooltip is up. The hover waits are
 * anchored on this rather than on a fixed pause: the tooltip is rebuilt a beat
 * after the pointer lands on a marker. */
const tooltipText = (page: Page): Promise<string> => page.evaluate(() => {
  const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
  if (!tt || tt.offsetParent === null) return '';
  return (tt.textContent ?? '').trim();
});

interface HoverResult {
  hitRow: number;
  entries: TooltipEntry[];
  aimOffsetX: number;
  aimOffsetY: number;
  viewportWidth: number;
  viewportHeight: number;
}

/**
 * Hover the marker of the reference row and read what came up. The pointer is
 * first parked elsewhere on the plot so that the move onto the marker is a real
 * change of position — the tooltip is rebuilt from that move. The aim goes
 * through the viewer's own world-to-screen mapping, which returns canvas-local
 * coordinates. The plot is dense, so the marker under that coordinate may belong
 * to a neighbouring row: the row the viewer reports is returned together with
 * its distance from the reference row, and the caller grades against that row.
 *
 * `expectTooltip: false` says the step configured the viewer to show none: an
 * absent tooltip cannot be waited for, so that case keeps a blind wait long
 * enough for one to have appeared.
 */
async function hoverReferenceMarker(
  page: Page, row: number, expectTooltip = true,
): Promise<HoverResult> {
  const rect = await canvasRect(page);
  await page.mouse.move(rect.x + rect.width * 0.03, rect.y + rect.height * 0.03);
  let parked = await tooltipText(page);
  for (let i = 0; i < 8 && parked.length > 0; i++) {
    await page.waitForTimeout(100);
    parked = await tooltipText(page);
  }
  const aim = await page.evaluate((i: number) => {
    const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
      .find((e) => !e.closest('.d4-dialog'))!;
    const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    const df = grok.shell.tv.dataFrame;
    const s = sp.worldToScreen(df.col(sp.props.xColumnName).get(i), df.col(sp.props.yColumnName).get(i));
    const b = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: b.x + s.x, y: b.y + s.y};
  }, row);
  await page.mouse.move(aim.x, aim.y, {steps: 8});
  if (expectTooltip) {
    for (let i = 0; i < 30; i++) {
      const t = await tooltipText(page);
      if (t.length > 0 && t !== parked) break;
      await page.waitForTimeout(100);
    }
    await page.waitForTimeout(200);
  } else
    await page.waitForTimeout(1600);
  const state = await page.evaluate((i: number) => {
    const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    const df = grok.shell.tv.dataFrame;
    const hit = df.mouseOverRowIdx as number;
    const xc = df.col(sp.props.xColumnName);
    const yc = df.col(sp.props.yColumnName);
    const vp = sp.viewport;
    return {
      hitRow: hit,
      aimOffsetX: hit >= 0 ? Math.abs(xc.get(hit) - xc.get(i)) : -1,
      aimOffsetY: hit >= 0 ? Math.abs(yc.get(hit) - yc.get(i)) : -1,
      viewportWidth: vp.width as number,
      viewportHeight: vp.height as number,
    };
  }, row);
  return {...state, entries: await tooltipEntries(page)};
}

const rowValues = (page: Page, row: number, columns: string[]) =>
  page.evaluate(({i, cols}: {i: number; cols: string[]}) => {
    const df = grok.shell.tv.dataFrame;
    const out: Record<string, string | null> = {};
    for (const c of cols) {
      const col = df.col(c);
      out[c] = col ? String(col.getString(i)) : null;
    }
    return out;
  }, {i: row, cols: columns});

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

async function clickCanvas(page: Page, p: Frac): Promise<void> {
  const r = await canvasRect(page);
  const pt = at(r, p);
  await page.mouse.click(pt.x, pt.y);
  await page.waitForTimeout(250);
}

const liveCount = (page: Page, kind: 'selection' | 'filter') =>
  page.evaluate((k: string) =>
    (k === 'selection' ? grok.shell.tv.dataFrame.selection.trueCount
      : grok.shell.tv.dataFrame.filter.trueCount) as number, kind);

/** Read a count once it stops moving: selection and filter are recomputed
 * asynchronously after a gesture. Used where the step is asserted to leave the
 * count as it is, so there is no known value to wait away from. */
async function settledCount(page: Page, kind: 'selection' | 'filter'): Promise<number> {
  let last = -1;
  let stable = 0;
  for (let i = 0; i < 30; i++) {
    const c = await liveCount(page, kind);
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

/** Read a count once it has left a known value and stopped moving. A count that
 * never leaves `from` is returned as it stands, so the caller's assert fails on
 * the real reading. */
async function settledCountAfterChange(
  page: Page, kind: 'selection' | 'filter', from: number, timeoutMs = 7000,
): Promise<number> {
  const deadline = Date.now() + timeoutMs;
  let prev = from;
  while (Date.now() < deadline) {
    const c = await liveCount(page, kind);
    if (c !== from) {
      if (c === prev) return c;
      prev = c;
    }
    await page.waitForTimeout(120);
  }
  return await liveCount(page, kind);
}

const rowCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);

/** Index of the filter card built for a column, found by the column name its
 * header carries — the cards are laid out in table order and this dataset's
 * panel is long, so position is not a usable handle. */
async function filterCardIndex(page: Page, column: string, timeoutMs = 45_000): Promise<number> {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    const idx = await page.evaluate((col: string) =>
      [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
        .findIndex((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col),
    column);
    if (idx >= 0) return idx;
    await page.waitForTimeout(300);
  }
  return -1;
}

interface CategoryState {filtered: number; survivors: string[]}

const categoryState = (page: Page, column: string): Promise<CategoryState> =>
  page.evaluate((col: string) => {
    const df = grok.shell.tv.dataFrame;
    const c = df.col(col);
    const counts: Record<string, number> = {};
    for (const x of c.categories) counts[x] = 0;
    for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[c.get(i)]++;
    return {
      filtered: df.filter.trueCount as number,
      survivors: Object.keys(counts).filter((k) => counts[k] > 0),
    };
  }, column);

/** Clear the filter cards. `from` is the filtered-row count the panel holds
 * before the click, so the release can be waited for instead of slept through. */
async function resetFilterPanel(page: Page, from?: number): Promise<void> {
  const reset = page.locator('[name="viewer-Filters"] [name="icon-arrow-rotate-left"]').first();
  await reset.scrollIntoViewIfNeeded();
  await reset.click();
  if (from === undefined) await page.waitForTimeout(800);
  else await settledCountAfterChange(page, 'filter', from, 8000);
}

/**
 * Leave exactly one category checked in a categorical filter card by clicking
 * its list. The list is painted on the card's canvas and carries no DOM text, so
 * the click cannot be aimed at a named category and the spec must not depend on
 * which one survives: a click anywhere on a category row narrows the filter to
 * that row's value. Several vertical positions are tried because the card's row
 * height follows its category count, and the filter is cleared between attempts
 * so each one starts from the whole table.
 */
async function narrowToOneCategory(
  page: Page, cardIndex: number, column: string, fullRows: number,
): Promise<CategoryState | null> {
  const card = page.locator('[name="viewer-Filters"] .d4-filter').nth(cardIndex);
  const canvas = card.locator('canvas').last();
  for (const fy of [0.5, 0.3, 0.7, 0.15, 0.85]) {
    await canvas.scrollIntoViewIfNeeded();
    const box = await canvas.boundingBox();
    if (!box || box.width === 0 || box.height === 0) break;
    await canvas.click({position: {x: box.width / 2, y: box.height * fy}});
    await settledCountAfterChange(page, 'filter', fullRows, 3500);
    const state = await categoryState(page, column);
    if (state.survivors.length === 1 && state.filtered > 0 && state.filtered < fullRows) return state;
    await resetFilterPanel(page, state.filtered === fullRows ? undefined : state.filtered);
  }
  return null;
}

test('Scatter Plot — Marker Labels and Tooltip', async ({page}: {page: Page}) => {
  test.setTimeout(1_500_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  const allMessages: string[] = [];
  page.on('pageerror', (e) => {
    allMessages.push(String(e));
    if (!isBenignError(String(e))) pageErrors.push(String(e));
  });
  page.on('console', (m) => {
    if (m.type() !== 'error') return;
    allMessages.push(m.text());
    if (!isBenignError(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;
  const crashCount = () => allMessages.filter((t) => CRASH_SIGNATURE.test(t)).length;

  await loginToDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
  await page.waitForTimeout(500);
  await pickOnViewer(page, 'x', SETUP_X);
  await pickOnViewer(page, 'y', SETUP_Y);
  expect(await readProp(page, 'xColumnName')).toBe(SETUP_X);
  expect(await readProp(page, 'yColumnName')).toBe(SETUP_Y);
  expect(await readProp(page, 'showTooltip')).toBe('inherit from table');
  expect(await readProp(page, 'labelColumnNames')).toEqual([]);

  await softStep('Tooltip inherited from the table', async () => {
    const errBefore = errCount();
    const hover = await hoverReferenceMarker(page, REFERENCE_ROW);
    expect(hover.hitRow).toBeGreaterThanOrEqual(0);
    expect(hover.aimOffsetX).toBeLessThanOrEqual(hover.viewportWidth * MARKER_AIM_TOLERANCE);
    expect(hover.aimOffsetY).toBeLessThanOrEqual(hover.viewportHeight * MARKER_AIM_TOLERANCE);

    const entries = hover.entries;
    expect(entries.length).toBeGreaterThan(1);
    const named = entries.filter((e) => e.name.length > 0);
    expect(named.length).toBe(entries.length);
    const table = await rowValues(page, hover.hitRow, named.map((e) => e.name));
    const unknown = named.filter((e) => table[e.name] === null).map((e) => e.name);
    const mismatched = named
      .filter((e) => table[e.name] !== null && table[e.name] !== e.value)
      .map((e) => `${e.name}: tooltip ${e.value} vs table ${table[e.name]}`);
    console.log(`inherited tooltip: row=${hover.hitRow} columns=${named.map((e) => e.name).join(',')}`);
    // Graded against whatever the table's tooltip column list happens to be:
    // every line must name a real column and repeat its value for the hit row.
    expect(unknown).toEqual([]);
    expect(mismatched).toEqual([]);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Custom tooltip column list', async () => {
    const errBefore = errCount();
    await setChoiceProp(page, 'prop-show-tooltip', 'prop-view-show-tooltip', 'tooltip', 'show custom tooltip');
    expect(await readProp(page, 'showTooltip')).toBe('show custom tooltip');
    await setListProp(page, 'rowTooltip', CUSTOM_TOOLTIP_COLUMNS.join('\n'));
    expect(await readProp(page, 'rowTooltip')).toBe(CUSTOM_TOOLTIP_COLUMNS.join('\n'));
    await setChoiceProp(page, 'prop-data-values', 'prop-view-data-values', 'tooltip', 'Do not add');
    expect(await readProp(page, 'dataValues')).toBe('Do not add');

    const custom = await hoverReferenceMarker(page, REFERENCE_ROW);
    expect(custom.hitRow).toBeGreaterThanOrEqual(0);
    const configured = await rowValues(page, custom.hitRow, CUSTOM_TOOLTIP_COLUMNS);
    console.log(`custom tooltip: row=${custom.hitRow} values=${custom.entries.map((e) => e.value).join('|')}`);
    // Nothing beyond the two configured values: the axis columns are excluded by
    // the Data Values setting.
    expect(custom.entries.map((e) => e.value))
      .toEqual(CUSTOM_TOOLTIP_COLUMNS.map((c) => configured[c]));

    await setChoiceProp(page, 'prop-show-tooltip', 'prop-view-show-tooltip', 'tooltip', 'do not show');
    expect(await readProp(page, 'showTooltip')).toBe('do not show');
    const silent = await hoverReferenceMarker(page, REFERENCE_ROW, false);
    expect(silent.hitRow).toBeGreaterThanOrEqual(0);
    expect(silent.entries).toEqual([]);

    await setChoiceProp(page, 'prop-show-tooltip', 'prop-view-show-tooltip', 'tooltip', 'inherit from table');
    await setListProp(page, 'rowTooltip', '');
    await setChoiceProp(page, 'prop-data-values', 'prop-view-data-values', 'tooltip', 'Merge');
    expect(await readProp(page, 'showTooltip')).toBe('inherit from table');
    expect(await readProp(page, 'rowTooltip')).toBe('');

    const restored = await hoverReferenceMarker(page, REFERENCE_ROW);
    expect(restored.hitRow).toBeGreaterThanOrEqual(0);
    expect(restored.entries.length).toBeGreaterThan(CUSTOM_TOOLTIP_COLUMNS.length);
    const namedAgain = restored.entries.filter((e) => e.name.length > 0);
    expect(namedAgain.length).toBe(restored.entries.length);
    const tableAgain = await rowValues(page, restored.hitRow, namedAgain.map((e) => e.name));
    expect(namedAgain.filter((e) => tableAgain[e.name] !== e.value)).toEqual([]);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Labels for the selected rows', async () => {
    const errBefore = errCount();
    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledCount(page, 'selection')).toBe(0);
    const baseline = await settledInk(page, 'overlay');
    // Comparing one fault value against another would let every claim below pass
    // without a measurement.
    expect(baseline).toBeGreaterThanOrEqual(0);

    await setListProp(page, 'labelColumnNames', [LABEL_COLUMN]);
    expect(await readProp(page, 'labelColumnNames')).toEqual([LABEL_COLUMN]);
    await setChoiceProp(page, 'prop-show-labels-for', 'prop-view-show-labels-for', 'labels', 'Selected');
    expect(await readProp(page, 'showLabelsFor')).toBe('Selected');
    // Configured but with nothing selected yet: this mode draws no label, so the
    // reading below measures the labels rather than the act of configuring them.
    const configured = await settledInk(page, 'overlay');
    expect(Math.abs(configured - baseline)).toBeLessThan(OVERLAY_RESTORE_TOLERANCE);

    await dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const selected = await settledCountAfterChange(page, 'selection', 0);
    expect(selected).toBeGreaterThan(0);
    const labelled = await settledInk(page, 'overlay');
    console.log(`label overlay ink: baseline=${baseline} configured=${configured} ` +
      `labelled=${labelled} selected=${selected}`);
    // The labels are painted on the overlay and contribute no DOM text, so the
    // label text is not asserted — only that the overlay rendering moved.
    expect(labelled - baseline).toBeGreaterThanOrEqual(LABEL_INK_DELTA);
    expect(errCount()).toBe(errBefore);

    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledCountAfterChange(page, 'selection', selected)).toBe(0);
    const cleared = await settledInk(page, 'overlay');
    console.log(`label overlay ink: cleared=${cleared}`);
    expect(Math.abs(cleared - baseline)).toBeLessThan(OVERLAY_RESTORE_TOLERANCE);

    await setListProp(page, 'labelColumnNames', []);
    await setChoiceProp(page, 'prop-show-labels-for', 'prop-view-show-labels-for', 'labels', 'All');
    expect(await readProp(page, 'labelColumnNames')).toEqual([]);
    expect(await readProp(page, 'showLabelsFor')).toBe('All');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Labels plus filtering on a datetime axis', async () => {
    await v.openTable(page, {path: spgiPath});
    const fullRows = await rowCount(page);
    expect(fullRows).toBeGreaterThan(0);
    await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
    await page.waitForTimeout(600);
    await pickOnViewer(page, 'x', SPGI_X);
    await pickOnViewer(page, 'y', SPGI_Y);
    expect(await readProp(page, 'xColumnName')).toBe(SPGI_X);
    expect(await readProp(page, 'yColumnName')).toBe(SPGI_Y);
    // The crash this scenario guards against comes from the label placement
    // arithmetic on a date axis, so the axis has to really be a date column.
    expect(await page.evaluate((c: string) =>
      String(grok.shell.tv.dataFrame.col(c).type), SPGI_X)).toBe('datetime');

    await setListProp(page, 'labelColumnNames', [SPGI_LABEL_COLUMN]);
    await setChoiceProp(page, 'prop-show-labels-for', 'prop-view-show-labels-for', 'labels', 'Selected');
    expect(await readProp(page, 'labelColumnNames')).toEqual([SPGI_LABEL_COLUMN]);
    expect(await readProp(page, 'showLabelsFor')).toBe('Selected');

    const errBefore = errCount();
    const crashBefore = crashCount();
    await v.openFilterPanel(page);
    const cardIndex = await filterCardIndex(page, SPGI_Y);
    expect(cardIndex).toBeGreaterThanOrEqual(0);
    const narrowed = await narrowToOneCategory(page, cardIndex, SPGI_Y, fullRows);
    expect(narrowed).not.toBeNull();
    // Which category the click leaves checked is up to the card's own ordering,
    // so it is read from the data and reported, never asserted; the graded fact
    // is that the table narrowed to a single category.
    console.log(`stereo category filter: kept=${narrowed!.survivors[0]} ` +
      `rows=${narrowed!.filtered}/${fullRows}`);
    expect(narrowed!.survivors.length).toBe(1);
    expect(narrowed!.filtered).toBeGreaterThan(0);
    expect(narrowed!.filtered).toBeLessThan(fullRows);

    await dragCanvas(page, {fx: 0.2, fy: 0.2}, {fx: 0.8, fy: 0.8}, ['Shift']);
    const selected = await settledCountAfterChange(page, 'selection', 0);
    expect(selected).toBeGreaterThan(0);
    await page.waitForTimeout(2000);
    expect(crashCount()).toBe(crashBefore);
    expect(errCount()).toBe(errBefore);

    await resetFilterPanel(page, narrowed!.filtered);
    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledCountAfterChange(page, 'selection', selected)).toBe(0);
    expect(await settledCount(page, 'filter')).toBe(fullRows);
    expect(crashCount()).toBe(crashBefore);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

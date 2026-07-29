/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const SETUP_X = 'WEIGHT';
const SETUP_Y = 'HEIGHT';
const WHISKER_X_MIN = 'AGE';
const WHISKER_X_MAX = 'WEIGHT';
const WHISKER_Y_MIN = 'HEIGHT';
const WHISKER_Y_MAX = 'WEIGHT';
const HISTOGRAM_BINS = 20;
const DEFAULT_HISTOGRAM_BINS = 10;
const TITLE_TEXT = 'Test Plot';
const DESCRIPTION_TEXT = 'Test description';
const LINES_ORDER_COLUMN = 'AGE';
const COLOR_COLUMN = 'RACE';
const LINES_BY_COLUMN = 'SEX';

// Canvas readings scale with viewer size and dpr, so only these bounds are
// asserted, never a reading. Their ordering is what keeps the asserts honest:
// settle < restore-ceiling < change-floor. The restore ceiling also bounds
// "the rendering did not change" (the connecting-lines split rules).
const CANVAS_SETTLE_TOLERANCE = 400;
const CANVAS_CHANGE_MIN = 2000;
const CANVAS_RESTORE_MAX = 1500;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

// Nothing beyond the shared server's ambient noise is whitelisted here: a Dart
// NullError or a "Stack trace" line while a visibility setting is toggled, a
// whisker column is bound or the context menu is opened is the regression signal.
const isBenignError = (text: string) => isAmbientError(text);

interface Rect {x: number; y: number; width: number; height: number}

/** Viewer root is disambiguated: a dialog carrying its own preview plot matches
 * the same name. */
const canvasRect = (page: Page): Promise<Rect> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

/** Poll `cond` until it holds. The caller keeps its own assert, so a poll that
 * runs out of time leaves the observed state to be graded as it is. */
async function waitUntil(cond: () => Promise<boolean>, timeoutMs: number, stepMs = 100): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  for (;;) {
    if (await cond()) return true;
    if (Date.now() >= deadline) return false;
    await new Promise((r) => setTimeout(r, stepMs));
  }
}

/** The plot repaints its mouse-over indication as the pointer crosses it, so
 * every rendering reading is taken from the same parked pointer position. The
 * longer settle serves the readings that have no known frame to leave. */
async function parkPointer(page: Page, settleMs = 300): Promise<void> {
  await page.mouse.move(2, 2);
  await page.waitForTimeout(settleMs);
}

/** Store a per-color histogram of the data canvas under `key`. False means the
 * canvas could not be read — a fault, not an equal reading. */
const captureCanvas = (page: Page, key: string) => page.evaluate((k: string) => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const c = root?.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
  const ctx = c?.getContext('2d');
  if (!c || !ctx) return false;
  let data: Uint8ClampedArray;
  try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { return false; }
  const colors = new Map<number, number>();
  for (let i = 0; i < data.length; i += 4) {
    const rgb = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
    colors.set(rgb, (colors.get(rgb) ?? 0) + 1);
  }
  const w = window as any;
  w.__spCanvasSnap = w.__spCanvasSnap || {};
  w.__spCanvasSnap[k] = colors;
  return true;
}, key);

/** Summed absolute per-color pixel delta against the reading stored under
 * `key`. -1 signals a fault. */
const diffCanvas = (page: Page, key: string) => page.evaluate((k: string) => {
  const w = window as any;
  const prev = w.__spCanvasSnap?.[k] as Map<number, number> | undefined;
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const c = root?.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
  const ctx = c?.getContext('2d');
  if (!prev || !c || !ctx) return -1;
  let data: Uint8ClampedArray;
  try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { return -1; }
  const colors = new Map<number, number>();
  for (let i = 0; i < data.length; i += 4) {
    const rgb = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
    colors.set(rgb, (colors.get(rgb) ?? 0) + 1);
  }
  let delta = 0;
  for (const [rgb, n] of colors) delta += Math.abs(n - (prev.get(rgb) ?? 0));
  for (const [rgb, n] of prev) if (!colors.has(rgb)) delta += n;
  return delta;
}, key);

/** Wait for the data canvas to stop repainting, then return its delta against
 * the reading stored under `key`. The settle gate compares two consecutive
 * readings of the current frame, so it stays independent of what the step
 * changed.
 *
 * `expectChange` is passed by the steps that claim the rendering moved: the poll
 * then starts at once and first waits for the delta against `key` to leave the
 * settle drift. The steps claiming the rendering did NOT move leave it off and
 * keep the blind settle, because an unchanged frame and one that has not
 * repainted yet read alike. */
async function settledCanvasDiff(page: Page, key: string, expectChange = false): Promise<number> {
  await parkPointer(page, expectChange ? 150 : 300);
  expect(await captureCanvas(page, '__settle')).toBe(true);
  if (expectChange) {
    for (let i = 0; i < 24; i++) {
      const moved = await diffCanvas(page, key);
      if (moved > CANVAS_SETTLE_TOLERANCE) break;
      await page.waitForTimeout(200);
    }
    await captureCanvas(page, '__settle');
  }
  for (let i = 0; i < 12; i++) {
    await page.waitForTimeout(expectChange ? 250 : 400);
    const drift = await diffCanvas(page, '__settle');
    await captureCanvas(page, '__settle');
    if (drift >= 0 && drift <= CANVAS_SETTLE_TOLERANCE) break;
  }
  const delta = await diffCanvas(page, key);
  expect(delta).toBeGreaterThanOrEqual(0);
  return delta;
}

async function captureBaseline(page: Page, key: string): Promise<void> {
  await settledCanvasDiff(page, '__settle');
  expect(await captureCanvas(page, key)).toBe(true);
}

/** Assign a column through the on-viewer selector. The shared helper verifies
 * the property itself, so the short commit settle is retried through a poll on
 * that property before the pick is repeated at the helper's own pace. */
async function pickOnViewer(page: Page, role: string, column: string): Promise<void> {
  try {
    await v.pickColumnViaSelectorTrusted(page, {role, columnName: column, commitSettleMs: 250});
    return;
  } catch (_) { /* fall through to the poll, then to a full retry */ }
  if (await waitUntil(async () => await readProp(page, `${role}ColumnName`) === column, 2500)) return;
  await v.pickColumnViaSelectorTrusted(page, {role, columnName: column});
}

// The helpers below serve PROPERTY-PANEL column rows only (whiskers, connecting
// lines). The shared pickColumnViaSelectorTrusted covers on-viewer selectors: it
// hovers the plot to reveal them and clicks their text. A panel row has neither
// — it is scrolled into view inside the settings panel and its popup opens from
// a mousedown on the nested combobox — so these paths stay local.
async function waitBackdrop(page: Page, timeout = 6000): Promise<boolean> {
  return await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout}).then(() => true).catch(() => false);
}

/** The popup grid is canvas-rendered, so the name goes in with real keys; the
 * first key is separated because the popup focuses its search input a tick
 * later. */
async function commitColumn(page: Page, column: string): Promise<void> {
  const text = column.toLowerCase();
  await page.keyboard.press(text[0]);
  await page.waitForTimeout(150);
  if (text.length > 1) await page.keyboard.type(text.slice(1));
  await page.waitForTimeout(200);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(200);
}

const settingsBuilt = (page: Page) =>
  page.evaluate(() => !!document.querySelector('[name="prop-category-axes"]'));

async function openSettings(page: Page): Promise<void> {
  for (let i = 0; i < 4; i++) {
    if (await settingsBuilt(page)) return;
    await v.openViewerGear(page, 'Scatter plot');
    await waitUntil(() => settingsBuilt(page), 2000);
  }
  throw new Error('the scatter plot settings panel did not build');
}

/** A row inside a collapsed category keeps its DOM node with an empty box, so
 * readiness is the editor's own rectangle, not the row's presence. The category
 * header is a toggle, hence the retry. */
async function revealPropEditor(page: Page, editorSelector: string, category: string): Promise<void> {
  const ready = () => page.evaluate((sel: string) => {
    const el = document.querySelector(sel) as HTMLElement | null;
    if (!el || !el.offsetParent) return false;
    const b = el.getBoundingClientRect();
    return b.width > 0 && b.height > 0;
  }, editorSelector);
  for (let i = 0; i < 8; i++) {
    if (await ready()) return;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();
    if (await waitUntil(ready, 1000)) return;
  }
  throw new Error(`property editor ${editorSelector} never became reachable`);
}

const readProp = (page: Page, name: string) => page.evaluate((n: string) => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp.props[n];
}, name);

/** Inline row opacity — the property grid's only greyed-out signal. */
const rowOpacity = (page: Page, rowName: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`) as HTMLElement | null)?.style.opacity ?? null, rowName);

const categoryCount = (page: Page, column: string) => page.evaluate((c: string) =>
  grok.shell.tv.dataFrame.col(c).categories.length as number, column);

const propIs = (page: Page, propName: string, value: unknown, timeoutMs: number) =>
  waitUntil(async () => await readProp(page, propName) === value, timeoutMs);

/** The row is the actuation; the viewer's own property confirms the click
 * landed before any rendering is read. */
async function setCheckboxProp(
  page: Page, rowName: string, category: string, propName: string, value: boolean,
): Promise<void> {
  await openSettings(page);
  const box = `[name="${rowName}"] input.property-grid-item-editor-checkbox`;
  await revealPropEditor(page, box, category);
  for (let i = 0; i < 3; i++) {
    if (await readProp(page, propName) === value) return;
    const locator = page.locator(box);
    await locator.scrollIntoViewIfNeeded();
    await locator.click();
    await propIs(page, propName, value, 2000);
  }
  if (await readProp(page, propName) !== value)
    throw new Error(`${propName} did not reach ${value} from the ${rowName} row`);
}

async function setNumericProp(
  page: Page, rowName: string, category: string, propName: string, value: number,
): Promise<void> {
  await openSettings(page);
  const selector = `[name="${rowName}"] input.property-grid-slider-textbox`;
  await revealPropEditor(page, selector, category);
  const locator = page.locator(selector).first();
  await locator.scrollIntoViewIfNeeded();
  await locator.click();
  await page.keyboard.press('Control+A');
  await page.keyboard.type(String(value));
  await page.keyboard.press('Enter');
  await propIs(page, propName, value, 2500);
  if (await readProp(page, propName) !== value)
    throw new Error(`${propName} did not reach ${value} from the ${rowName} row`);
}

/** The popup of a property-panel column row opens on a mousedown, never on a
 * plain click. */
async function openPanelColumnPopup(
  page: Page, rowName: string, comboName: string, category: string,
): Promise<void> {
  await openSettings(page);
  const combo = `[name="${rowName}"] [name="${comboName}"]`;
  await revealPropEditor(page, combo, category);
  await page.locator(combo).scrollIntoViewIfNeeded();
  await page.evaluate(({row, name}: {row: string; name: string}) => {
    const sel = document.querySelector(`[name="${row}"] [name="${name}"]`) as HTMLElement;
    const target = (sel.querySelector('.d4-column-selector-column') as HTMLElement) || sel;
    target.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
  }, {row: rowName, name: comboName});
  if (!await waitBackdrop(page)) throw new Error(`${comboName} popup did not open`);
}

async function pickPanelColumn(
  page: Page, rowName: string, comboName: string, category: string, propName: string, column: string,
): Promise<void> {
  await openPanelColumnPopup(page, rowName, comboName, category);
  await commitColumn(page, column);
  await propIs(page, propName, column, 2500);
  if (await readProp(page, propName) !== column)
    throw new Error(`${propName} did not reach ${column} from the ${rowName} row`);
}

/** The popup's first data row is the empty entry, one row height below the
 * popup's own header; picking it sets the property to the empty string. */
async function clearPanelColumn(
  page: Page, rowName: string, comboName: string, category: string, propName: string,
): Promise<void> {
  await openPanelColumnPopup(page, rowName, comboName, category);
  const row = await page.evaluate(() => {
    const b = document.querySelector('.d4-column-selector-backdrop')!.getBoundingClientRect();
    return {x: b.x + b.width / 2, y: b.y + 30};
  });
  await page.mouse.click(row.x, row.y);
  await waitUntil(async () => {
    const p = await readProp(page, propName);
    return p === '' || p === null;
  }, 2500);
  const left = await readProp(page, propName);
  if (left !== '' && left !== null)
    throw new Error(`${propName} was left as ${left} by the ${rowName} row`);
}

/** The value cell turns into a text editor only once clicked. */
async function setTextProp(
  page: Page, rowName: string, viewCell: string, propName: string, value: string,
): Promise<void> {
  await openSettings(page);
  await revealPropEditor(page, `[name="${viewCell}"]`, 'description');
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  const editor = page.locator(
    `[name="${rowName}"] input.property-grid-item-editor-textbox, ` +
    `[name="${rowName}"] input.property-grid-ellipsis-editor-input`).first();
  await waitUntil(async () => await editor.count() > 0 && await editor.isVisible(), 2000);
  await editor.click();
  await page.keyboard.press('Control+A');
  if (value === '') await page.keyboard.press('Delete');
  else await page.keyboard.type(value);
  await page.keyboard.press('Enter');
  await propIs(page, propName, value, 2500);
  if (await readProp(page, propName) !== value)
    throw new Error(`${propName} did not reach "${value}" from the ${rowName} row`);
}

/** State of the two on-viewer axis column selectors. The GROK-13533 guard turns
 * on the computed visibility: DOM presence, a live offsetParent and a non-zero
 * box all survive the setting being turned off, so each alone passes either way. */
const axisSelectorState = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const read = (role: string) => {
    const el = root.querySelector(`[name="div-column-combobox-${role}"]`) as HTMLElement | null;
    if (!el) return null;
    const b = el.getBoundingClientRect();
    return {
      visibility: getComputedStyle(el).visibility,
      inDom: true,
      hasOffsetParent: !!el.offsetParent,
      width: b.width,
      height: b.height,
    };
  };
  return {x: read('x'), y: read('y')};
});

/** Scoped to `.d4-menu-popup` — a bare `.d4-menu-item` query also matches the
 * application menubar. */
const menuEntryNames = (page: Page) => page.evaluate(() =>
  [...document.querySelectorAll('.d4-menu-popup [name]')]
    .map((e) => e.getAttribute('name')!).filter((n) => !!n));

const viewerAlive = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  return !!sp && !!root && document.body.contains(sp.root);
});

test('Scatter Plot — Secondary Settings Surface', async ({page}: {page: Page}) => {
  test.setTimeout(900_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
  await waitUntil(() => page.evaluate(() => {
    const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
      .find((e) => !e.closest('.d4-dialog'));
    const c = root?.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
    return !!c && c.getBoundingClientRect().width > 0;
  }), 15000);

  await pickOnViewer(page, 'x', SETUP_X);
  await pickOnViewer(page, 'y', SETUP_Y);
  await openSettings(page);

  await softStep('Scenario 1 — Axis histograms', async () => {
    const errBefore = errCount();
    await captureBaseline(page, 'hist-base');

    await setCheckboxProp(page, 'prop-show-x-histogram', 'axes', 'showXHistogram', true);
    expect(await settledCanvasDiff(page, 'hist-base', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    expect(await captureCanvas(page, 'hist-x')).toBe(true);
    await setCheckboxProp(page, 'prop-show-y-histogram', 'axes', 'showYHistogram', true);
    expect(await settledCanvasDiff(page, 'hist-x', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    expect(await captureCanvas(page, 'hist-xy')).toBe(true);
    await setNumericProp(page, 'prop-histogram-bins', 'axes', 'histogramBins', HISTOGRAM_BINS);
    expect(await settledCanvasDiff(page, 'hist-xy', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    await setCheckboxProp(page, 'prop-show-x-histogram', 'axes', 'showXHistogram', false);
    await setCheckboxProp(page, 'prop-show-y-histogram', 'axes', 'showYHistogram', false);
    await setNumericProp(page, 'prop-histogram-bins', 'axes', 'histogramBins', DEFAULT_HISTOGRAM_BINS);
    expect(await settledCanvasDiff(page, 'hist-base')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 — Grid lines, axes and selector visibility (GROK-13533)', async () => {
    const errBefore = errCount();
    await captureBaseline(page, 'vis-base');
    const before = await axisSelectorState(page);
    expect(before.x?.visibility).toBe('visible');
    expect(before.y?.visibility).toBe('visible');

    await setCheckboxProp(page, 'prop-show-vertical-grid-lines', 'x', 'showVerticalGridLines', false);
    await setCheckboxProp(page, 'prop-show-horizontal-grid-lines', 'y', 'showHorizontalGridLines', false);
    expect(await settledCanvasDiff(page, 'vis-base', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    expect(await captureCanvas(page, 'vis-nogrid')).toBe(true);
    await setCheckboxProp(page, 'prop-show-x-axis', 'x', 'showXAxis', false);
    await setCheckboxProp(page, 'prop-show-y-axis', 'y', 'showYAxis', false);
    expect(await settledCanvasDiff(page, 'vis-nogrid', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    await setCheckboxProp(page, 'prop-show-x-selector', 'x', 'showXSelector', false);
    await setCheckboxProp(page, 'prop-show-y-selector', 'y', 'showYSelector', false);
    await parkPointer(page);
    const hidden = await axisSelectorState(page);
    // The setting hides the selectors without removing them: computed
    // visibility is the only reading that separates the two states — the three
    // a presence-style check would use are asserted unchanged next to it.
    expect(hidden.x?.visibility).toBe('hidden');
    expect(hidden.y?.visibility).toBe('hidden');
    expect(hidden.x?.inDom).toBe(true);
    expect(hidden.y?.inDom).toBe(true);
    expect(hidden.x?.hasOffsetParent).toBe(true);
    expect(hidden.y?.hasOffsetParent).toBe(true);
    expect(hidden.x!.width).toBeGreaterThan(0);
    expect(hidden.x!.height).toBeGreaterThan(0);
    expect(hidden.y!.width).toBeGreaterThan(0);
    expect(hidden.y!.height).toBeGreaterThan(0);

    await setCheckboxProp(page, 'prop-show-x-selector', 'x', 'showXSelector', true);
    await setCheckboxProp(page, 'prop-show-y-selector', 'y', 'showYSelector', true);
    await setCheckboxProp(page, 'prop-show-x-axis', 'x', 'showXAxis', true);
    await setCheckboxProp(page, 'prop-show-y-axis', 'y', 'showYAxis', true);
    await setCheckboxProp(page, 'prop-show-vertical-grid-lines', 'x', 'showVerticalGridLines', true);
    await setCheckboxProp(page, 'prop-show-horizontal-grid-lines', 'y', 'showHorizontalGridLines', true);
    await parkPointer(page);
    const restored = await axisSelectorState(page);
    expect(restored.x?.visibility).toBe('visible');
    expect(restored.y?.visibility).toBe('visible');
    expect(await settledCanvasDiff(page, 'vis-base')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 3 — Whiskers', async () => {
    const errBefore = errCount();
    await captureBaseline(page, 'whisker-base');

    await pickPanelColumn(page, 'prop-x-whisker-min', 'div-column-combobox-x--whisker--min',
      'x', 'xWhiskerMinColumnName', WHISKER_X_MIN);
    await pickPanelColumn(page, 'prop-x-whisker-max', 'div-column-combobox-x--whisker--max',
      'x', 'xWhiskerMaxColumnName', WHISKER_X_MAX);
    await pickPanelColumn(page, 'prop-y-whisker-min', 'div-column-combobox-y--whisker--min',
      'y', 'yWhiskerMinColumnName', WHISKER_Y_MIN);
    await pickPanelColumn(page, 'prop-y-whisker-max', 'div-column-combobox-y--whisker--max',
      'y', 'yWhiskerMaxColumnName', WHISKER_Y_MAX);
    expect(await settledCanvasDiff(page, 'whisker-base', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);
    expect(errCount()).toBe(errBefore);

    await clearPanelColumn(page, 'prop-x-whisker-min', 'div-column-combobox-x--whisker--min',
      'x', 'xWhiskerMinColumnName');
    await clearPanelColumn(page, 'prop-x-whisker-max', 'div-column-combobox-x--whisker--max',
      'x', 'xWhiskerMaxColumnName');
    await clearPanelColumn(page, 'prop-y-whisker-min', 'div-column-combobox-y--whisker--min',
      'y', 'yWhiskerMinColumnName');
    await clearPanelColumn(page, 'prop-y-whisker-max', 'div-column-combobox-y--whisker--max',
      'y', 'yWhiskerMaxColumnName');
    expect(await settledCanvasDiff(page, 'whisker-base')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 4 — Context menu', async () => {
    const errBefore = errCount();
    const r = await canvasRect(page);
    await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
    await page.locator('.d4-menu-popup').last().waitFor({timeout: 10000});
    await waitUntil(async () => (await menuEntryNames(page)).includes('div-Properties...'), 3000);

    const entries = await menuEntryNames(page);
    expect(entries).toContain('div-Reset-View');
    expect(entries).toContain('div-Lasso-Tool');
    expect(entries).toContain('div-Tools');
    expect(entries).toContain('div-Properties...');

    await page.keyboard.press('Escape');
    await page.locator('.d4-menu-popup').first().waitFor({state: 'detached', timeout: 8000});
    expect(await page.locator('.d4-menu-popup').count()).toBe(0);
    expect(await viewerAlive(page)).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 5 — Title and description', async () => {
    const errBefore = errCount();
    // The title renders on the viewer's title bar, the description inside the
    // viewer element, so the two texts are read from their own scopes.
    const viewerText = () => page.evaluate(() => {
      const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
        .find((e) => !e.closest('.d4-dialog'))! as HTMLElement;
      const panel = root.closest('.panel-base') as HTMLElement;
      return {
        inViewer: (root.innerText ?? '').replace(/\s+/g, ' ').trim(),
        onViewer: (panel.innerText ?? '').replace(/\s+/g, ' ').trim(),
      };
    });

    const before = await viewerText();
    expect(before.onViewer).not.toContain(TITLE_TEXT);
    expect(before.inViewer).not.toContain(DESCRIPTION_TEXT);

    await setTextProp(page, 'prop-title', 'prop-view-title', 'title', TITLE_TEXT);
    expect((await viewerText()).onViewer).toContain(TITLE_TEXT);

    await setTextProp(page, 'prop-description', 'prop-view-description', 'description', DESCRIPTION_TEXT);
    const withBoth = await viewerText();
    expect(withBoth.onViewer).toContain(TITLE_TEXT);
    expect(withBoth.inViewer).toContain(DESCRIPTION_TEXT);

    await setTextProp(page, 'prop-title', 'prop-view-title', 'title', '');
    await setTextProp(page, 'prop-description', 'prop-view-description', 'description', '');
    const cleared = await viewerText();
    expect(cleared.onViewer).not.toContain(TITLE_TEXT);
    expect(cleared.inViewer).not.toContain(DESCRIPTION_TEXT);

    expect(errCount()).toBe(errBefore);
  });

  await softStep('Lines By overrides the color column when splitting connecting lines', async () => {
    const errBefore = errCount();
    // Lines By stays inert until a Lines Order column gives the lines an order.
    expect(await rowOpacity(page, 'prop-lines-by')).toBe('0.5');
    await captureBaseline(page, 'lines-base');

    await pickPanelColumn(page, 'prop-color', 'div-column-combobox-color',
      'color', 'colorColumnName', COLOR_COLUMN);
    // The two candidate split columns must really differ in how many series
    // they produce.
    const colorCategories = await categoryCount(page, COLOR_COLUMN);
    const linesByCategories = await categoryCount(page, LINES_BY_COLUMN);
    expect(linesByCategories).toBeGreaterThan(1);
    expect(colorCategories).toBeGreaterThan(linesByCategories);
    expect(await captureCanvas(page, 'lines-nolines')).toBe(true);

    await pickPanelColumn(page, 'prop-lines-order', 'div-column-combobox-lines--order',
      'data', 'linesOrderColumnName', LINES_ORDER_COLUMN);
    expect(await rowOpacity(page, 'prop-lines-by')).toBe('1');
    expect(await settledCanvasDiff(page, 'lines-nolines', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);
    // Lines split by the color column — the reading both Lines By settings are
    // judged against.
    expect(await captureCanvas(page, 'lines-color-split')).toBe(true);

    // Pointing Lines By at the color column keeps the split key the same, so the
    // picture must not move — this separates a changed split key from the mere
    // act of setting a property.
    await pickPanelColumn(page, 'prop-lines-by', 'div-column-combobox-lines--by',
      'data', 'linesByColumnName', COLOR_COLUMN);
    expect(await settledCanvasDiff(page, 'lines-color-split')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    // A column with fewer categories changes the split key, and the picture with
    // it — Lines By has taken over from the color column.
    await pickPanelColumn(page, 'prop-lines-by', 'div-column-combobox-lines--by',
      'data', 'linesByColumnName', LINES_BY_COLUMN);
    expect(await settledCanvasDiff(page, 'lines-color-split', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    await clearPanelColumn(page, 'prop-lines-by', 'div-column-combobox-lines--by',
      'data', 'linesByColumnName');
    expect(await settledCanvasDiff(page, 'lines-color-split')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    await clearPanelColumn(page, 'prop-lines-order', 'div-column-combobox-lines--order',
      'data', 'linesOrderColumnName');
    await clearPanelColumn(page, 'prop-color', 'div-column-combobox-color',
      'color', 'colorColumnName');
    expect(await settledCanvasDiff(page, 'lines-base')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    expect(errCount()).toBe(errBefore);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

/* ---
realizes: [scatterplot.cp.regression-and-formula-lines, viewers.scatter-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const SETUP_X = 'WEIGHT';
const SETUP_Y = 'HEIGHT';
const ALT_X = 'AGE';
const ALT_Y = 'WEIGHT';
const COLOR_COLUMN = 'RACE';
const DATETIME_X = 'STARTED';
const TIME_UNIT = 'year';
// The dialog accepts a formula only in the form `${column} = expression`, so a
// left-hand side that is not a single column is rejected by design.
const REJECTED_FORMULA = '${HEIGHT} - ${WEIGHT} = 1';
const EDITED_FORMULA = '${HEIGHT} = ${WEIGHT} + 1';
const MOVING_AVERAGE_WINDOW = 200;

// Canvas readings scale with viewer size and dpr, so only these bounds are
// asserted, never a reading. Their ordering keeps the asserts honest:
// INK_SETTLE < the two tolerances < the two change floors.
const INK_SETTLE = 60;
const REGRESSION_OVERLAY_DELTA = 500;
const MOVING_AVERAGE_DELTA = 300;
const OVERLAY_STABLE_TOLERANCE = 100;
const DATA_RESTORE_TOLERANCE = 150;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

// No project save here, so no Dart NullError or "Stack trace" window is
// expected: nothing beyond the shared server's ambient noise is whitelisted.
const isBenignError = (text: string) => isAmbientError(text);

interface Rect {x: number; y: number; width: number; height: number}

/** Rect of the scatter-plot data canvas, resolved on the viewer root that is NOT
 * inside a dialog — the Formula Lines dialog embeds its own preview plot with a
 * canvas of the same name. */
const canvasRect = (page: Page): Promise<Rect> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

/** Non-white opaque pixels on one of the plot's two canvases. `canvas` carries
 * the markers, axes and the moving-average line with its deviation band;
 * `overlay` carries the regression line and its equation box. -1 signals a
 * canvas that cannot be read, which the caller treats as a fault. */
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

/** Park the pointer clear of the plot: a hovered point paints its mouse-over
 * indicator into the frame and spoils the measurement. The longer settle serves
 * the readings that have no known value to leave. */
async function parkPointer(page: Page, settleMs = 500): Promise<void> {
  await page.mouse.move(4, 4);
  await page.waitForTimeout(settleMs);
}

/** Read a canvas until two consecutive readings agree within the settle drift —
 * a repaint finishes after the property that triggered it has already been set.
 *
 * `from` names a reading the layer is expected to leave: the poll then starts at
 * once and waits for the departure instead of sleeping. It is passed only where
 * the step claims the layer changed. A step claiming the layer did NOT change
 * omits it and keeps the blind settle, because an unchanged layer and one that
 * has not repainted yet read alike. */
async function settledInk(page: Page, layer: 'canvas' | 'overlay', from?: number): Promise<number> {
  await parkPointer(page, from === undefined ? 500 : 200);
  let prev = await canvasInk(page, layer);
  let cur = prev;
  if (from !== undefined) {
    for (let i = 0; i < 24 && !(cur >= 0 && Math.abs(cur - from) > INK_SETTLE); i++) {
      await page.waitForTimeout(150);
      cur = await canvasInk(page, layer);
    }
    prev = cur;
  }
  for (let i = 0; i < 12; i++) {
    await page.waitForTimeout(from === undefined ? 300 : 200);
    cur = await canvasInk(page, layer);
    if (cur >= 0 && Math.abs(cur - prev) <= INK_SETTLE) break;
    prev = cur;
  }
  return cur;
}

/** The data canvas carries the moving average, the overlay the regression line.
 * `fromData` is the data reading the step expects to leave; the overlay is read
 * blind throughout, since the claim made about it is that it stayed put. */
async function settledBoth(page: Page, fromData?: number): Promise<{data: number; overlay: number}> {
  const data = await settledInk(page, 'canvas', fromData);
  const overlay = await settledInk(page, 'overlay');
  return {data, overlay};
}

/** Move the real pointer over the plot. The on-viewer selectors carry
 * `visibility: hidden` until the pointer is over the viewer, so clearing a
 * column and reaching the time-unit control both start here. With a role given,
 * the wait is on that selector turning visible. */
async function hoverPlot(page: Page, role?: string): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  if (!role) {
    await page.waitForTimeout(400);
    return;
  }
  await waitUntil(() => page.evaluate((n: string) => {
    const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
      .find((e) => !e.closest('.d4-dialog'));
    const el = root?.querySelector(`[name="div-column-combobox-${n}"]`) as HTMLElement | null;
    return !!el && getComputedStyle(el).visibility !== 'hidden';
  }, role), 3000);
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

// The two helpers below serve the CLEAR path only. The shared
// pickColumnViaSelectorTrusted helper assigns a column and verifies it on the
// viewer — the opposite operation, with no empty-row coordinate.
/** Center point of an on-viewer column selector's clickable text. An unassigned
 * selector exposes only its caption, so the caption is the fallback target. */
const selectorPoint = (page: Page, role: string) =>
  page.evaluate((r: string) => {
    const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
      .find((e) => !e.closest('.d4-dialog'))!;
    const sel = root.querySelector(`[name="div-column-combobox-${r}"]`) as HTMLElement | null;
    if (!sel) return null;
    const candidates: Element[] = [];
    for (const q of ['.d4-column-selector-column', '.d4-column-selector-caption']) {
      const el = sel.querySelector(q);
      if (el) candidates.push(el);
    }
    candidates.push(sel);
    for (const el of candidates) {
      const b = el.getBoundingClientRect();
      if (b.width > 0 && b.height > 0) return {x: b.x + b.width / 2, y: b.y + b.height / 2};
    }
    return null;
  }, role);

async function waitBackdrop(page: Page, timeout = 5000): Promise<boolean> {
  return await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout}).then(() => true).catch(() => false);
}

/** Clear a column through the same on-viewer selector: the popup's first data
 * row is the empty entry, one row height below its own header. */
async function clearOnViewer(page: Page, role: string): Promise<void> {
  await hoverPlot(page, role);
  const pt = await selectorPoint(page, role);
  if (!pt) throw new Error(`on-viewer ${role} selector has no clickable text`);
  await page.mouse.click(pt.x, pt.y);
  if (!await waitBackdrop(page)) throw new Error(`${role} column popup did not open`);
  const row = await page.evaluate(() => {
    const b = document.querySelector('.d4-column-selector-backdrop')!.getBoundingClientRect();
    return {x: b.x + b.width / 2, y: b.y + 30};
  });
  await page.mouse.click(row.x, row.y);
  await waitUntil(async () => {
    const left = await readProp(page, `${role}ColumnName`);
    return left === '' || left === null;
  }, 2500);
}

const settingsBuilt = (page: Page) =>
  page.evaluate(() => !!document.querySelector('[name="prop-category-data"]'));

async function openSettings(page: Page): Promise<void> {
  for (let i = 0; i < 4; i++) {
    if (await settingsBuilt(page)) return;
    await v.openViewerGear(page, 'Scatter plot');
    await waitUntil(() => settingsBuilt(page), 2000);
  }
  throw new Error('the scatter plot settings panel did not build');
}

/** Make a property editor reachable. A row inside a collapsed property-grid
 * category keeps its DOM node with an empty box, so the check is on the editor's
 * own rectangle; the category header is a toggle, so an attempt that collapsed
 * it is simply retried. */
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

const propIs = (page: Page, propName: string, value: unknown, timeoutMs: number) =>
  waitUntil(async () => await readProp(page, propName) === value, timeoutMs);

/** Set a choice property: the value cell turns into a select once clicked. */
async function setChoiceProp(
  page: Page, rowName: string, viewCell: string, category: string, value: string, propName: string,
): Promise<void> {
  await openSettings(page);
  await revealPropEditor(page, `[name="${viewCell}"]`, category);
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  const select = page.locator(`[name="${rowName}"] select.property-grid-item-editor-spinner`);
  await waitUntil(async () => await select.count() > 0 && await select.isVisible(), 2000);
  await select.selectOption(value);
  await propIs(page, propName, value, 2500);
}

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

/** Set a numeric property. The editor is either a slider textbox or a plain
 * number input depending on whether the property declares bounds. */
async function setNumericProp(
  page: Page, rowName: string, category: string, propName: string, value: number,
): Promise<void> {
  await openSettings(page);
  const box = await page.evaluate((n: string) => {
    const row = document.querySelector(`[name="${n}"]`);
    if (!row) return null;
    return row.querySelector('input.property-grid-slider-textbox') ? 'slider' :
      row.querySelector('input') ? 'input' : null;
  }, rowName);
  if (!box) throw new Error(`no numeric editor inside ${rowName}`);
  const selector = box === 'slider' ?
    `[name="${rowName}"] input.property-grid-slider-textbox` : `[name="${rowName}"] input`;
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

const readProp = (page: Page, name: string) => page.evaluate((n: string) => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp.props[n];
}, name);

const viewerProps = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {
    x: sp.props.xColumnName, y: sp.props.yColumnName, color: sp.props.colorColumnName,
    xMap: sp.props.xMap, yAxisType: sp.props.yAxisType,
    showRegressionLine: sp.props.showRegressionLine,
    regressionPerCategory: sp.props.regressionPerCategory,
    showMovingAverageLine: sp.props.showMovingAverageLine,
    movingAverageWindow: sp.props.movingAverageWindow,
    movingAveragePerCategory: sp.props.movingAveragePerCategory,
    showMovingAverageDeviation: sp.props.showMovingAverageDeviation,
  };
});

/** The stored formula-line configuration — the only product-state signal the
 * dialog has: its list of lines is a canvas-rendered grid with no DOM text. */
const formulaLines = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  try { return JSON.parse(sp.props.formulaLines || '[]') as any[]; } catch (_) { return null; }
});

async function openPlotContextMenu(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await waitUntil(async () => (await menuLeafNames(page)).includes('div-Tools'), 3000);
}

const menuLeafNames = (page: Page) => page.evaluate(() =>
  [...document.querySelectorAll('.d4-menu-popup [name]')]
    .map((e) => e.getAttribute('name')!).filter((n) => !!n));

const FORMULA_DIALOG = '[name="dialog-Formula-Lines"]';

/** Values of the dialog's formula editors — the DOM signal that its list of
 * annotations gained or lost an entry. */
const formulaEditorValues = (page: Page) => page.evaluate((dlg: string) =>
  ([...document.querySelectorAll(`${dlg} textarea`)] as HTMLTextAreaElement[])
    .map((a) => a.value).join('|'), FORMULA_DIALOG);

async function openFormulaLinesDialog(page: Page): Promise<void> {
  await openPlotContextMenu(page);
  await page.locator('.d4-menu-popup [name="div-Tools"]').last().hover();
  await waitUntil(async () => (await menuLeafNames(page)).includes('div-Tools---Formula-Lines...'), 3000);
  expect(await menuLeafNames(page)).toContain('div-Tools---Formula-Lines...');
  await page.locator('.d4-menu-popup [name="div-Tools---Formula-Lines..."]').last().click();
  await page.locator(FORMULA_DIALOG).waitFor({timeout: 15000});
  const addNew = page.locator(`${FORMULA_DIALOG} [name="button-Add-new"]`);
  await waitUntil(async () => await addNew.count() > 0 && await addNew.isVisible(), 5000);
}

/** Close the dialog with OK and wait for it to leave the DOM — every reading of
 * the main plot happens after this, because the dialog carries its own plot. */
async function closeFormulaLinesDialog(page: Page): Promise<void> {
  await page.locator(`${FORMULA_DIALOG} [name="button-OK"]`).click();
  await page.locator(FORMULA_DIALOG).waitFor({state: 'detached', timeout: 15000});
  await waitUntil(async () => await page.locator('.d4-dialog').count() === 0, 3000);
}

/** Leave no dialog behind: a step that ends early can leave the Formula Lines
 * dialog on screen, and a modal dialog blocks the plot's own context menu. */
async function dismissStrayDialog(page: Page): Promise<void> {
  for (let i = 0; i < 3; i++) {
    if (await page.locator('.d4-dialog').count() === 0) return;
    await page.keyboard.press('Escape');
    await waitUntil(async () => await page.locator('.d4-dialog').count() === 0, 1500);
  }
}

/** Add an annotation from the dialog's ADD NEW menu — every entry is a
 * top-level leaf there. */
async function addNewAnnotation(page: Page, leaf: string): Promise<void> {
  const before = await formulaEditorValues(page);
  await page.locator(`${FORMULA_DIALOG} [name="button-Add-new"]`).click();
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await waitUntil(async () => (await menuLeafNames(page)).includes(leaf), 3000);
  expect(await menuLeafNames(page)).toContain(leaf);
  await page.locator(`.d4-menu-popup [name="${leaf}"]`).last().click();
  await waitUntil(async () => await formulaEditorValues(page) !== before, 4000);
}

/** Index of the formula editor inside the dialog: the textarea holding a column
 * reference. The editor carries no name attribute of its own. */
const formulaEditorIndex = (page: Page) => page.evaluate((dlg: string) => {
  const areas = [...document.querySelectorAll(`${dlg} textarea`)] as HTMLTextAreaElement[];
  return areas.findIndex((a) => /\$\{/.test(a.value));
}, FORMULA_DIALOG);

const formulaEditorState = (page: Page) => page.evaluate((dlg: string) => {
  const areas = [...document.querySelectorAll(`${dlg} textarea`)] as HTMLTextAreaElement[];
  const a = areas.find((t) => /\$\{/.test(t.value));
  if (!a) return {found: false, value: '', invalid: false};
  return {found: true, value: a.value, invalid: a.classList.contains('d4-forced-invalid')};
}, FORMULA_DIALOG);

/** Type a formula into the dialog's editor. `settleInvalid` names the validation
 * outcome the step goes on to grade — the wait is bounded, so an outcome that
 * never arrives is graded as it stands. */
async function setFormula(page: Page, formula: string, settleInvalid: boolean): Promise<void> {
  const idx = await formulaEditorIndex(page);
  expect(idx).toBeGreaterThanOrEqual(0);
  const editor = page.locator(`${FORMULA_DIALOG} textarea`).nth(idx);
  const wanted = formula.replace(/\s+/g, '');
  await editor.fill(formula);
  await waitUntil(async () => {
    const s = await formulaEditorState(page);
    return s.value.replace(/\s+/g, '') === wanted && s.invalid === settleInvalid;
  }, 3000);
  const state = await formulaEditorState(page);
  expect(state.value.replace(/\s+/g, '')).toBe(formula.replace(/\s+/g, ''));
}

/** Whether the dialog's OK button accepts a click — the dialog refuses to
 * commit a formula its validation rejected. */
const okEnabled = (page: Page) => page.evaluate((dlg: string) =>
  !(document.querySelector(`${dlg} [name="button-OK"]`)?.className ?? '').split(/\s+/).includes('disabled'),
FORMULA_DIALOG);

/** Delete every configured annotation. The stored configuration only updates
 * when the dialog is committed, so the count is taken before it opens and the
 * emptiness is verified after OK. */
async function deleteAllFormulaLines(page: Page): Promise<void> {
  for (let attempt = 0; attempt < 3; attempt++) {
    const lines = await formulaLines(page);
    if (lines !== null && lines.length === 0) return;
    await openFormulaLinesDialog(page);
    for (let i = 0; i < (lines?.length ?? 1); i++) {
      const row = await page.evaluate((dlg: string) => {
        const g = document.querySelector(`${dlg} [name="viewer-Grid"]`);
        if (!g) return null;
        const b = g.getBoundingClientRect();
        return {x: b.x + b.width / 2, y: b.y + 30};
      }, FORMULA_DIALOG);
      if (row) {
        await page.mouse.click(row.x, row.y);
        await page.waitForTimeout(300);
      }
      // Each configured item carries its own Delete button; the first one acts
      // on the item the grid currently points at.
      const before = await formulaEditorValues(page);
      await page.locator(`${FORMULA_DIALOG} [name="button-Delete"]`).first().click();
      await waitUntil(async () => await formulaEditorValues(page) !== before, 2500);
    }
    await closeFormulaLinesDialog(page);
  }
}

/** Page coordinates of a data point, aimed with the viewer's own world-to-screen
 * mapping, which returns canvas-local coordinates. */
const pointAt = (page: Page, row: number) => page.evaluate((i: number) => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const df = grok.shell.tv.dataFrame;
  const s = sp.worldToScreen(df.col(sp.props.xColumnName).get(i), df.col(sp.props.yColumnName).get(i));
  const b = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: b.x + s.x, y: b.y + s.y};
}, row);

const tooltipText = (page: Page) => page.evaluate(() =>
  ((document.querySelector('.d4-tooltip') as HTMLElement | null)?.innerText ?? '').trim());

/** Display style of the on-viewer time-unit selector — `none` for numeric and
 * categorical columns, laid out for a datetime column. */
const timeUnitDisplay = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const el = root.querySelector('[name="input-aggr-selector-x-map"]') as HTMLElement | null;
  return el ? getComputedStyle(el).display : null;
});

/** Choose a time unit in the on-viewer selector. The control is a native select
 * that ignores a scripted value assignment, so the option is picked through the
 * driver's own select interaction. */
async function chooseTimeUnit(page: Page, unit: string): Promise<boolean> {
  await hoverPlot(page);
  await page.locator('[name="viewer-Scatter-plot"] [name="input-aggr-selector-x-map"]').first()
    .selectOption(unit);
  await propIs(page, 'xMap', unit, 2500);
  return await readProp(page, 'xMap') === unit;
}

test('Scatter Plot — Regression Line, Formula Lines, Moving Average', async ({page}: {page: Page}) => {
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
  await waitUntil(() => page.evaluate(() => {
    const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
      .find((e) => !e.closest('.d4-dialog'));
    const c = root?.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
    return !!c && c.getBoundingClientRect().width > 0;
  }), 15000);
  await pickOnViewer(page, 'x', SETUP_X);
  await pickOnViewer(page, 'y', SETUP_Y);
  // Both per-category switches ship on, so they are normalized to off: a ladder
  // step that toggles a switch already in its target state grades nothing.
  await setCheckboxProp(page, 'prop-regression-per-category', 'lines', 'regressionPerCategory', false);
  await setCheckboxProp(page, 'prop-moving-average-per-category', 'lines',
    'movingAveragePerCategory', false);
  const setup = await viewerProps(page);
  expect(setup.x).toBe(SETUP_X);
  expect(setup.y).toBe(SETUP_Y);
  expect(setup.showRegressionLine).toBe(false);
  expect(setup.color === '' || setup.color === null).toBe(true);
  expect(await formulaLines(page)).toEqual([]);

  await softStep('Regression line on a logarithmic Y axis', async () => {
    const errBefore = errCount();
    await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'logarithmic', 'yAxisType');
    expect((await viewerProps(page)).yAxisType).toBe('logarithmic');

    // Baseline taken on the logarithmic axis itself, so the delta below can only
    // come from the regression line and its statistics box.
    const base = await settledInk(page, 'overlay');
    expect(base).toBeGreaterThan(0);

    await setCheckboxProp(page, 'prop-show-regression-line', 'lines', 'showRegressionLine', true);
    const drawn = await settledInk(page, 'overlay', base);
    console.log(`regression overlay ink: logBase=${base} logOn=${drawn}`);
    // The equation and statistics box is painted on this same overlay and
    // exposes no readable text and no accessor, so the coefficients are
    // deliberately not asserted.
    expect(Math.abs(drawn - base)).toBeGreaterThan(REGRESSION_OVERLAY_DELTA);
    expect(errCount()).toBe(errBefore);

    await setCheckboxProp(page, 'prop-show-regression-line', 'lines', 'showRegressionLine', false);
    await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'linear', 'yAxisType');
    const restored = await settledInk(page, 'overlay', drawn);
    console.log(`regression overlay ink: restored=${restored}`);
    const back = await viewerProps(page);
    expect(back.showRegressionLine).toBe(false);
    expect(back.yAxisType).toBe('linear');
    expect(Math.abs(restored - base)).toBeLessThan(Math.abs(restored - drawn));
  });

  await softStep('Regression Per Category and the axis time unit with the regression line on', async () => {
    const errBefore = errCount();
    await pickOnViewer(page, 'color', COLOR_COLUMN);
    expect((await viewerProps(page)).color).toBe(COLOR_COLUMN);
    await setCheckboxProp(page, 'prop-show-regression-line', 'lines', 'showRegressionLine', true);

    // The greying of the per-category row is governed by Show Regression Line
    // rather than by the number of colour categories, so no disabled-state claim
    // is made here — only a no-error floor over the toggle.
    await setCheckboxProp(page, 'prop-regression-per-category', 'lines', 'regressionPerCategory', true);
    expect(await settledInk(page, 'overlay')).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);

    await setCheckboxProp(page, 'prop-regression-per-category', 'lines', 'regressionPerCategory', false);
    expect(await settledInk(page, 'overlay')).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);

    await pickOnViewer(page, 'x', DATETIME_X);
    expect((await viewerProps(page)).x).toBe(DATETIME_X);
    // Moving a datetime column onto the axis switches the regression line off,
    // so it is turned back on before the time unit is chosen — the guard is
    // about a time unit picked WITH the line drawn.
    await setCheckboxProp(page, 'prop-show-regression-line', 'lines', 'showRegressionLine', true);
    expect(await timeUnitDisplay(page)).not.toBe('none');
    expect(await chooseTimeUnit(page, TIME_UNIT)).toBe(true);
    expect((await viewerProps(page)).xMap).toBe(TIME_UNIT);
    expect(await settledInk(page, 'canvas')).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);

    expect(await chooseTimeUnit(page, '')).toBe(true);
    await pickOnViewer(page, 'x', SETUP_X);
    await setCheckboxProp(page, 'prop-show-regression-line', 'lines', 'showRegressionLine', false);
    await clearOnViewer(page, 'color');
    const reverted = await viewerProps(page);
    expect(reverted.x).toBe(SETUP_X);
    expect(reverted.xMap).toBe('');
    expect(reverted.showRegressionLine).toBe(false);
    expect(reverted.regressionPerCategory).toBe(false);
    expect(reverted.color).toBe('');
  });

  await softStep('Formula lines dialog — add a line, edit it across the equals sign, reopen', async () => {
    await dismissStrayDialog(page);
    await openFormulaLinesDialog(page);
    await addNewAnnotation(page, 'div-Line');
    const prefilled = await formulaEditorState(page);
    expect(prefilled.found).toBe(true);
    expect(prefilled.value).toContain(SETUP_Y);
    expect(prefilled.value).toContain(SETUP_X);
    expect(prefilled.invalid).toBe(false);
    expect(await okEnabled(page)).toBe(true);

    // The rejected shape: the editor flags the input and holds the commit back,
    // while the form stays alive and editable.
    await setFormula(page, REJECTED_FORMULA, true);
    const rejected = await formulaEditorState(page);
    expect(rejected.invalid).toBe(true);
    expect(await okEnabled(page)).toBe(false);

    // Correcting the edit clears the flag — the editor does not stay stuck in
    // the state the rejected input put it in.
    await setFormula(page, EDITED_FORMULA, false);
    const accepted = await formulaEditorState(page);
    expect(accepted.invalid).toBe(false);
    expect(await okEnabled(page)).toBe(true);
    await closeFormulaLinesDialog(page);

    await openFormulaLinesDialog(page);
    const reopened = await formulaEditorState(page);
    expect(reopened.found).toBe(true);
    expect(reopened.value.replace(/\s+/g, '')).toBe(EDITED_FORMULA.replace(/\s+/g, ''));
    expect(reopened.invalid).toBe(false);

    const lines = await formulaLines(page);
    expect(lines).not.toBeNull();
    expect(lines!.length).toBe(1);
    expect(String(lines![0].formula).replace(/\s+/g, '')).toBe(EDITED_FORMULA.replace(/\s+/g, ''));
    // The line stays in place for the next scenario.
    await closeFormulaLinesDialog(page);
  });

  await softStep('Formula line across an axis-column change, and a band on a logarithmic axis', async () => {
    await dismissStrayDialog(page);
    const before = await formulaLines(page);
    expect(before).not.toBeNull();
    expect(before!.length).toBe(1);

    await pickOnViewer(page, 'x', ALT_X);
    await pickOnViewer(page, 'y', ALT_Y);
    const moved = await viewerProps(page);
    expect(moved.x).toBe(ALT_X);
    expect(moved.y).toBe(ALT_Y);
    expect(await formulaLines(page)).toEqual(before);

    await openFormulaLinesDialog(page);
    await addNewAnnotation(page, 'div-Band---Horizontal');
    await closeFormulaLinesDialog(page);
    const withBand = await formulaLines(page);
    expect(withBand!.length).toBe(2);
    expect(withBand!.some((l: any) => String(l.type).includes('band'))).toBe(true);

    const errBefore = errCount();
    await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'logarithmic', 'yAxisType');
    expect((await viewerProps(page)).yAxisType).toBe('logarithmic');
    expect(await settledInk(page, 'canvas')).toBeGreaterThan(0);
    await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'linear', 'yAxisType');
    expect((await viewerProps(page)).yAxisType).toBe('linear');
    expect(await settledInk(page, 'canvas')).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);

    await deleteAllFormulaLines(page);
    expect(await formulaLines(page)).toEqual([]);
    await pickOnViewer(page, 'x', SETUP_X);
    await pickOnViewer(page, 'y', SETUP_Y);
    const restored = await viewerProps(page);
    expect(restored.x).toBe(SETUP_X);
    expect(restored.y).toBe(SETUP_Y);
  });

  await softStep('Hover sweep with a formula line present', async () => {
    await dismissStrayDialog(page);
    await openFormulaLinesDialog(page);
    await addNewAnnotation(page, 'div-Line');
    await closeFormulaLinesDialog(page);
    expect((await formulaLines(page))!.length).toBe(1);

    const errBefore = errCount();
    let tooltips = 0;
    for (const row of [10, 40, 120, 300, 900]) {
      const pt = await pointAt(page, row);
      await page.mouse.move(pt.x, pt.y, {steps: 6});
      await waitUntil(async () => (await tooltipText(page)).length > 0, 1200);
      if ((await tooltipText(page)).length > 0) tooltips++;
      await page.mouse.move(pt.x + 4, pt.y + 4, {steps: 3});
      await page.waitForTimeout(250);
    }
    expect(tooltips).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);

    await deleteAllFormulaLines(page);
    expect(await formulaLines(page)).toEqual([]);
  });

  await softStep('Moving average line, window, per category and deviation', async () => {
    await dismissStrayDialog(page);
    const errBefore = errCount();
    await pickOnViewer(page, 'color', COLOR_COLUMN);
    // The regression line comes back on when an axis column moves back to a
    // numeric one, so this scenario puts it off itself: the moving-average
    // ladder is measured on a canvas the regression overlay does not touch.
    await setCheckboxProp(page, 'prop-show-regression-line', 'lines', 'showRegressionLine', false);
    const start = await viewerProps(page);
    expect(start.color).toBe(COLOR_COLUMN);
    expect(start.showRegressionLine).toBe(false);
    const baseline = await settledBoth(page);
    expect(baseline.data).toBeGreaterThan(0);
    // Both baselines have to be real readings: comparing one fault value against
    // another would let the "overlay unchanged" claim below pass unmeasured.
    expect(baseline.overlay).toBeGreaterThan(0);
    const defaultWindow = start.movingAverageWindow as number;
    expect(defaultWindow).toBeGreaterThan(0);

    await setCheckboxProp(page, 'prop-show-moving-average-line', 'lines', 'showMovingAverageLine', true);
    const lineOn = await settledBoth(page, baseline.data);
    expect(Math.abs(lineOn.data - baseline.data)).toBeGreaterThan(MOVING_AVERAGE_DELTA);

    await setNumericProp(page, 'prop-moving-average-window', 'lines',
      'movingAverageWindow', MOVING_AVERAGE_WINDOW);
    const widened = await settledBoth(page, lineOn.data);
    expect(Math.abs(widened.data - lineOn.data)).toBeGreaterThan(MOVING_AVERAGE_DELTA);

    await setCheckboxProp(page, 'prop-moving-average-per-category', 'lines',
      'movingAveragePerCategory', true);
    const perCategory = await settledBoth(page, widened.data);
    expect(Math.abs(perCategory.data - widened.data)).toBeGreaterThan(MOVING_AVERAGE_DELTA);

    await setCheckboxProp(page, 'prop-show-moving-average-deviation', 'lines',
      'showMovingAverageDeviation', true);
    const deviation = await settledBoth(page, perCategory.data);
    const steps = [
      Math.abs(lineOn.data - baseline.data),
      Math.abs(widened.data - lineOn.data),
      Math.abs(perCategory.data - widened.data),
    ];
    const deviationStep = Math.abs(deviation.data - perCategory.data);
    console.log(`moving average data ink: base=${baseline.data} line=${lineOn.data} ` +
      `window=${widened.data} perCat=${perCategory.data} deviation=${deviation.data}`);
    console.log(`moving average overlay ink: base=${baseline.overlay} line=${lineOn.overlay} ` +
      `window=${widened.overlay} perCat=${perCategory.overlay} deviation=${deviation.overlay}`);
    // The deviation band is the widest surface the ladder adds.
    expect(deviationStep).toBeGreaterThan(Math.max(...steps));

    // The moving average draws on the data canvas, unlike the regression line,
    // and that separation is the claim.
    for (const reading of [lineOn, widened, perCategory, deviation])
      expect(Math.abs(reading.overlay - baseline.overlay)).toBeLessThan(OVERLAY_STABLE_TOLERANCE);
    expect(errCount()).toBe(errBefore);

    await setCheckboxProp(page, 'prop-show-moving-average-deviation', 'lines',
      'showMovingAverageDeviation', false);
    await setCheckboxProp(page, 'prop-moving-average-per-category', 'lines',
      'movingAveragePerCategory', false);
    await setCheckboxProp(page, 'prop-show-moving-average-line', 'lines',
      'showMovingAverageLine', false);
    await setNumericProp(page, 'prop-moving-average-window', 'lines',
      'movingAverageWindow', defaultWindow);
    const restored = await settledBoth(page, deviation.data);
    console.log(`moving average data ink: restored=${restored.data}`);
    expect(Math.abs(restored.data - baseline.data)).toBeLessThan(DATA_RESTORE_TOLERANCE);

    await clearOnViewer(page, 'color');
    const final = await viewerProps(page);
    expect(final.color).toBe('');
    expect(final.showMovingAverageLine).toBe(false);
    expect(final.movingAveragePerCategory).toBe(false);
    expect(final.showMovingAverageDeviation).toBe(false);
    expect(final.movingAverageWindow).toBe(defaultWindow);
    expect(final.x).toBe(SETUP_X);
    expect(final.y).toBe(SETUP_Y);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

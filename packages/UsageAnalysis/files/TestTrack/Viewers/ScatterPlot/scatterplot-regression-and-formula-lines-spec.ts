/* ---
realizes: [scatterplot.cp.regression-and-formula-lines, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as sp from './scatterplot-shared';

declare const grok: any;

// Server lane: the Formula Lines dialog (steps 3-5) is contributed by the PowerPack package, which
// the local client does not load.
test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const SETUP_X = 'WEIGHT';
const SETUP_Y = 'HEIGHT';
const ALT_X = 'AGE';
const ALT_Y = 'WEIGHT';
const COLOR_COLUMN = 'RACE';
const DATETIME_X = 'STARTED';
const TIME_UNIT = 'year';

const REJECTED_FORMULA = '${HEIGHT} - ${WEIGHT} = 1';
const EDITED_FORMULA = '${HEIGHT} = ${WEIGHT} + 1';
const MOVING_AVERAGE_WINDOW = 200;

const INK_SETTLE = 60;
const REGRESSION_OVERLAY_DELTA = 500;
const MOVING_AVERAGE_DELTA = 300;
const OVERLAY_STABLE_TOLERANCE = 100;
const DATA_RESTORE_TOLERANCE = 150;

const settledInk = (page: Page, layer: 'canvas' | 'overlay', from?: number) =>
  sp.settledInk(page, layer, INK_SETTLE, from);

async function settledBoth(page: Page, fromData?: number): Promise<{data: number; overlay: number}> {
  const data = await settledInk(page, 'canvas', fromData);
  const overlay = await settledInk(page, 'overlay');
  return {data, overlay};
}

async function hoverPlot(page: Page): Promise<void> {
  const r = await sp.canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  await v.waitForViewerRendered(page, sp.SP_TYPE, 400);
}

const setYAxisType = (page: Page, value: 'logarithmic' | 'linear') =>
  sp.setChoiceProp(page, 'prop-y-axis-type', 'y-axis', value, 'yAxisType');

const setLinesCheckbox = (page: Page, rowName: string, value: boolean) =>
  sp.setCheckboxProp(page, rowName, 'lines', value);

const viewerProps = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {
    x: s.props.xColumnName, y: s.props.yColumnName, color: s.props.colorColumnName,
    xMap: s.props.xMap, yAxisType: s.props.yAxisType,
    showRegressionLine: s.props.showRegressionLine,
    regressionPerCategory: s.props.regressionPerCategory,
    showMovingAverageLine: s.props.showMovingAverageLine,
    movingAverageWindow: s.props.movingAverageWindow,
    movingAveragePerCategory: s.props.movingAveragePerCategory,
    showMovingAverageDeviation: s.props.showMovingAverageDeviation,
  };
});

const formulaLines = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  try { return JSON.parse(s.props.formulaLines || '[]') as any[]; } catch (_) { return null; }
});

const FORMULA_DIALOG = '[name="dialog-Formula-Lines"]';

const formulaEditorValues = (page: Page) => page.evaluate((dlg: string) =>
  ([...document.querySelectorAll(`${dlg} textarea`)] as HTMLTextAreaElement[])
    .map((a) => a.value).join('|'), FORMULA_DIALOG);

async function openFormulaLinesDialog(page: Page): Promise<void> {
  await sp.openPlotContextMenu(page);
  await sp.navigatePopup(page, ['div-Tools', 'div-Tools---Formula-Lines...']);
  await page.locator(FORMULA_DIALOG).waitFor({timeout: 15000});
  await page.locator(`${FORMULA_DIALOG} [name="button-Add-new"]`).waitFor({state: 'visible', timeout: 5000});
}

async function closeFormulaLinesDialog(page: Page): Promise<void> {
  await page.locator(`${FORMULA_DIALOG} [name="button-OK"]`).click();
  await page.locator(FORMULA_DIALOG).waitFor({state: 'detached', timeout: 15000});
  await v.pollValue(() => page.locator('.d4-dialog').count(), (n) => n === 0, 3000, 50);
}

async function dismissStrayDialog(page: Page): Promise<void> {
  for (let i = 0; i < 3 && await page.locator('.d4-dialog').count() > 0; i++) {
    await page.keyboard.press('Escape');
    await v.pollValue(() => page.locator('.d4-dialog').count(), (n) => n === 0, 1500, 50);
  }
}

async function addNewAnnotation(page: Page, leaf: string): Promise<void> {
  const before = await formulaEditorValues(page);
  await page.locator(`${FORMULA_DIALOG} [name="button-Add-new"]`).click();
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await sp.navigatePopup(page, [leaf]);
  await v.pollValue(() => formulaEditorValues(page), (cur) => cur !== before, 4000, 50);
}

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

async function setFormula(page: Page, formula: string, settleInvalid: boolean): Promise<void> {
  const idx = await formulaEditorIndex(page);
  expect(idx).toBeGreaterThanOrEqual(0);
  const editor = page.locator(`${FORMULA_DIALOG} textarea`).nth(idx);
  const wanted = formula.replace(/\s+/g, '');
  await editor.fill(formula);
  await v.pollValue(() => formulaEditorState(page),
    (s) => s.value.replace(/\s+/g, '') === wanted && s.invalid === settleInvalid, 3000, 50);
  const state = await formulaEditorState(page);
  expect(state.value.replace(/\s+/g, '')).toBe(wanted);
}

const okEnabled = (page: Page) => page.evaluate((dlg: string) =>
  !(document.querySelector(`${dlg} [name="button-OK"]`)?.className ?? '').split(/\s+/).includes('disabled'),
FORMULA_DIALOG);

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
        await v.waitForViewerRendered(page, sp.SP_TYPE, 300);
      }
      const before = await formulaEditorValues(page);
      await page.locator(`${FORMULA_DIALOG} [name="button-Delete"]`).first().click();
      await v.pollValue(() => formulaEditorValues(page), (cur) => cur !== before, 2500, 50);
    }
    await closeFormulaLinesDialog(page);
  }
}

const pointAt = (page: Page, row: number) => page.evaluate((i: number) => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const df = grok.shell.tv.dataFrame;
  const pt = s.worldToScreen(df.col(s.props.xColumnName).get(i), df.col(s.props.yColumnName).get(i));
  const b = s.root.querySelector('canvas[name="canvas"]').getBoundingClientRect();
  return {x: b.x + pt.x, y: b.y + pt.y};
}, row);

const tooltipText = (page: Page) => page.evaluate(() =>
  ((document.querySelector('.d4-tooltip') as HTMLElement | null)?.innerText ?? '').trim());

const timeUnitDisplay = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const el = s.root.querySelector('[name="input-aggr-selector-x-map"]') as HTMLElement | null;
  return el ? getComputedStyle(el).display : null;
});

async function chooseTimeUnit(page: Page, unit: string): Promise<boolean> {
  await hoverPlot(page);
  await page.locator('[name="viewer-Scatter-plot"] [name="input-aggr-selector-x-map"]').first()
    .selectOption(unit);
  return await sp.propIs(page, 'xMap', unit, 2500) === unit;
}

test('Scatter Plot — Regression Line, Formula Lines, Moving Average', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const errors = sp.trackErrors(page);
  const errCount = errors.count;
  const errorsSince = (n: number) => errors.all().slice(n);

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await sp.addScatterPlot(page);
  await sp.pickOnViewer(page, 'x', SETUP_X);
  await sp.pickOnViewer(page, 'y', SETUP_Y);

  await setLinesCheckbox(page, 'prop-regression-per-category', false);
  await setLinesCheckbox(page, 'prop-moving-average-per-category', false);
  const setup = await viewerProps(page);
  expect(setup.x).toBe(SETUP_X);
  expect(setup.y).toBe(SETUP_Y);
  expect(setup.showRegressionLine).toBe(false);
  expect(setup.color === '' || setup.color === null).toBe(true);
  expect(await formulaLines(page)).toEqual([]);

  await softStep('Regression line on a logarithmic Y axis', async () => {
    const errBefore = errCount();
    await setYAxisType(page, 'logarithmic');
    expect((await viewerProps(page)).yAxisType).toBe('logarithmic');

    const base = await settledInk(page, 'overlay');
    expect(base).toBeGreaterThan(0);

    await setLinesCheckbox(page, 'prop-show-regression-line', true);
    const drawn = await settledInk(page, 'overlay', base);
    console.log(`regression overlay ink: logBase=${base} logOn=${drawn}`);

    expect(Math.abs(drawn - base)).toBeGreaterThan(REGRESSION_OVERLAY_DELTA);
    expect(errorsSince(errBefore)).toEqual([]);

    await setLinesCheckbox(page, 'prop-show-regression-line', false);
    await setYAxisType(page, 'linear');
    const restored = await settledInk(page, 'overlay', drawn);
    console.log(`regression overlay ink: restored=${restored}`);
    const back = await viewerProps(page);
    expect(back.showRegressionLine).toBe(false);
    expect(back.yAxisType).toBe('linear');
    expect(Math.abs(restored - base)).toBeLessThan(Math.abs(restored - drawn));
  });

  await softStep('Regression Per Category and the axis time unit with the regression line on', async () => {
    const errBefore = errCount();
    await sp.pickOnViewer(page, 'color', COLOR_COLUMN);
    expect((await viewerProps(page)).color).toBe(COLOR_COLUMN);
    await setLinesCheckbox(page, 'prop-show-regression-line', true);

    await setLinesCheckbox(page, 'prop-regression-per-category', true);
    expect(await settledInk(page, 'overlay')).toBeGreaterThan(0);
    expect(errorsSince(errBefore)).toEqual([]);

    await setLinesCheckbox(page, 'prop-regression-per-category', false);
    expect(await settledInk(page, 'overlay')).toBeGreaterThan(0);
    expect(errorsSince(errBefore)).toEqual([]);

    await sp.pickOnViewer(page, 'x', DATETIME_X);
    expect((await viewerProps(page)).x).toBe(DATETIME_X);

    await setLinesCheckbox(page, 'prop-show-regression-line', true);
    expect(await timeUnitDisplay(page)).not.toBe('none');
    expect(await chooseTimeUnit(page, TIME_UNIT)).toBe(true);
    expect((await viewerProps(page)).xMap).toBe(TIME_UNIT);
    expect(await settledInk(page, 'canvas')).toBeGreaterThan(0);
    expect(errorsSince(errBefore)).toEqual([]);

    expect(await chooseTimeUnit(page, '')).toBe(true);
    await sp.pickOnViewer(page, 'x', SETUP_X);
    await setLinesCheckbox(page, 'prop-show-regression-line', false);
    await sp.clearOnViewer(page, 'color');
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

    await setFormula(page, REJECTED_FORMULA, true);
    const rejected = await formulaEditorState(page);
    expect(rejected.invalid).toBe(true);
    expect(await okEnabled(page)).toBe(false);

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

    await closeFormulaLinesDialog(page);
  });

  await softStep('Formula line across an axis-column change, and a band on a logarithmic axis', async () => {
    await dismissStrayDialog(page);
    const before = await formulaLines(page);
    expect(before).not.toBeNull();
    expect(before!.length).toBe(1);

    await sp.pickOnViewer(page, 'x', ALT_X);
    await sp.pickOnViewer(page, 'y', ALT_Y);
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
    await setYAxisType(page, 'logarithmic');
    expect((await viewerProps(page)).yAxisType).toBe('logarithmic');
    expect(await settledInk(page, 'canvas')).toBeGreaterThan(0);
    await setYAxisType(page, 'linear');
    expect((await viewerProps(page)).yAxisType).toBe('linear');
    expect(await settledInk(page, 'canvas')).toBeGreaterThan(0);
    expect(errorsSince(errBefore)).toEqual([]);

    await deleteAllFormulaLines(page);
    expect(await formulaLines(page)).toEqual([]);
    await sp.pickOnViewer(page, 'x', SETUP_X);
    await sp.pickOnViewer(page, 'y', SETUP_Y);
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
      const shown = await v.armEvent(page, 'grok.events.onTooltipShown', 1200);
      await page.mouse.move(pt.x, pt.y, {steps: 6});
      await shown();
      if ((await tooltipText(page)).length > 0) tooltips++;
      const shownAgain = await v.armEvent(page, 'grok.events.onTooltipShown', 250);
      await page.mouse.move(pt.x + 4, pt.y + 4, {steps: 3});
      await shownAgain();
    }
    expect(tooltips).toBeGreaterThan(0);
    expect(errorsSince(errBefore)).toEqual([]);

    await deleteAllFormulaLines(page);
    expect(await formulaLines(page)).toEqual([]);
  });

  await softStep('Moving average line, window, per category and deviation', async () => {
    await dismissStrayDialog(page);
    const errBefore = errCount();
    await sp.pickOnViewer(page, 'color', COLOR_COLUMN);

    await setLinesCheckbox(page, 'prop-show-regression-line', false);
    const start = await viewerProps(page);
    expect(start.color).toBe(COLOR_COLUMN);
    expect(start.showRegressionLine).toBe(false);
    const baseline = await settledBoth(page);
    expect(baseline.data).toBeGreaterThan(0);
    expect(baseline.overlay).toBeGreaterThan(0);
    const defaultWindow = start.movingAverageWindow as number;
    expect(defaultWindow).toBeGreaterThan(0);

    await setLinesCheckbox(page, 'prop-show-moving-average-line', true);
    const lineOn = await settledBoth(page, baseline.data);
    expect(Math.abs(lineOn.data - baseline.data)).toBeGreaterThan(MOVING_AVERAGE_DELTA);

    await sp.setNumericProp(page, 'prop-moving-average-window', 'lines', MOVING_AVERAGE_WINDOW);
    const widened = await settledBoth(page, lineOn.data);
    expect(Math.abs(widened.data - lineOn.data)).toBeGreaterThan(MOVING_AVERAGE_DELTA);

    await setLinesCheckbox(page, 'prop-moving-average-per-category', true);
    const perCategory = await settledBoth(page, widened.data);
    expect(Math.abs(perCategory.data - widened.data)).toBeGreaterThan(MOVING_AVERAGE_DELTA);

    await setLinesCheckbox(page, 'prop-show-moving-average-deviation', true);
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

    expect(deviationStep).toBeGreaterThan(Math.max(...steps));

    for (const reading of [lineOn, widened, perCategory, deviation])
      expect(Math.abs(reading.overlay - baseline.overlay)).toBeLessThan(OVERLAY_STABLE_TOLERANCE);
    expect(errorsSince(errBefore)).toEqual([]);

    await setLinesCheckbox(page, 'prop-show-moving-average-deviation', false);
    await setLinesCheckbox(page, 'prop-moving-average-per-category', false);
    await setLinesCheckbox(page, 'prop-show-moving-average-line', false);
    await sp.setNumericProp(page, 'prop-moving-average-window', 'lines', defaultWindow);
    const restored = await settledBoth(page, deviation.data);
    console.log(`moving average data ink: restored=${restored.data}`);
    expect(Math.abs(restored.data - baseline.data)).toBeLessThan(DATA_RESTORE_TOLERANCE);

    await sp.clearOnViewer(page, 'color');
    const final = await viewerProps(page);
    expect(final.color).toBe('');
    expect(final.showMovingAverageLine).toBe(false);
    expect(final.movingAveragePerCategory).toBe(false);
    expect(final.showMovingAverageDeviation).toBe(false);
    expect(final.movingAverageWindow).toBe(defaultWindow);
    expect(final.x).toBe(SETUP_X);
    expect(final.y).toBe(SETUP_Y);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});

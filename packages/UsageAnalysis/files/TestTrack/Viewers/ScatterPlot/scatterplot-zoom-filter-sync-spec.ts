/* ---
realizes: [scatterplot.cp.zoom-filter-sync, viewers.scatter-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const PROBE_X = 'ZZ_LOG_X_PROBE';
const PROBE_Y = 'ZZ_LOG_Y_PROBE';
const SPGI_X = 'CAST Idea ID';
const SPGI_Y = 'Idea ID';
const JITTER_X = 16;
const JITTER_Y = 17;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

// Nothing beyond the shared-server ambient noise is whitelisted here: a Dart
// NullError or a "Stack trace" line during a zoom, a reset or an axis-type
// switch is the regression signal.
const isBenignError = (text: string) => isAmbientError(text);

/** Viewer root is disambiguated: a dialog can embed its own preview plot under
 * the same name. */
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

/** The shared helper reads the property back and raises a named error when the
 * pick did not take, so callers do not wrap it in an assertion. */
const pickOnViewer = (page: Page, role: string, column: string) =>
  v.pickColumnViaSelectorTrusted(page, {role, columnName: column});

/** The property grid is rebuilt from scratch whenever the active table view
 * changes, hence the retry. */
async function openSettings(page: Page): Promise<void> {
  for (let i = 0; i < 4; i++) {
    const built = await page.evaluate(() => !!document.querySelector('[name="prop-category-data"]'));
    if (built) return;
    await v.openViewerGear(page, 'Scatter plot');
    await page.waitForTimeout(1000);
  }
  throw new Error('the scatter plot settings panel did not build');
}

/** A row inside a collapsed category keeps its DOM node with an empty box, so
 * readiness is the editor's own rectangle, not the row's presence. The category
 * header is a toggle, hence the retry. */
async function revealPropEditor(
  page: Page, editorSelector: string, category: string,
): Promise<void> {
  for (let i = 0; i < 8; i++) {
    const ready = await page.evaluate((sel: string) => {
      const el = document.querySelector(sel) as HTMLElement | null;
      if (!el || !el.offsetParent) return false;
      const b = el.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    }, editorSelector);
    if (ready) return;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();
    await page.waitForTimeout(800);
  }
  throw new Error(`property editor ${editorSelector} never became reachable`);
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

/** The viewer property behind a settings-panel row: `prop-jitter-size-y` is
 * `jitterSizeY`. */
const rowProp = (rowName: string) =>
  rowName.replace(/^prop-/, '').replace(/-([a-z])/g, (_, c: string) => c.toUpperCase());

const viewerProp = (page: Page, prop: string) => page.evaluate((p: string) => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp ? sp.props[p] : null;
}, prop);

/** The value cell turns into a select only once clicked. */
async function setChoiceProp(
  page: Page, rowName: string, viewCell: string, category: string, value: string,
): Promise<void> {
  await openSettings(page);
  await revealPropEditor(page, `[name="${viewCell}"]`, category);
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  await page.locator(`[name="${rowName}"] select.property-grid-item-editor-spinner`).selectOption(value);
  await waitFor(page, async () => await propCellText(page, viewCell) === value, 4000);
  await page.waitForTimeout(200);
}

async function setSliderProp(
  page: Page, rowName: string, category: string, value: number,
): Promise<void> {
  await openSettings(page);
  const box = `[name="${rowName}"] input.property-grid-slider-textbox`;
  await revealPropEditor(page, box, category);
  const locator = page.locator(box);
  await locator.scrollIntoViewIfNeeded();
  await locator.click();
  await page.keyboard.press('Control+A');
  await page.keyboard.type(String(value));
  await page.keyboard.press('Enter');
  await waitFor(page, async () => await viewerProp(page, rowProp(rowName)) === value, 4000);
  await page.waitForTimeout(200);
}

const propertyRowNames = (page: Page) => page.evaluate(() =>
  [...document.querySelectorAll('tr.property-grid-item[name]')]
    .map((e) => e.getAttribute('name') as string));

/** The settings-panel checkbox `checked` state mirrors the property. */
async function setFilterOutInvalid(page: Page, on: boolean): Promise<void> {
  await openSettings(page);
  await revealPropEditor(page, '[name="prop-view-filter-out-invalid"]', 'data');
  const box = page.locator('input[name="prop-view-filter-out-invalid"]');
  await box.scrollIntoViewIfNeeded();
  if (await box.isChecked() !== on) {
    await box.click();
    await waitFor(page, async () => await box.isChecked() === on, 3000);
  }
  expect(await box.isChecked()).toBe(on);
}

const propCellText = (page: Page, viewCell: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`)?.textContent ?? '').trim(), viewCell);

/** Inline row opacity — the only disabled signal the settings panel exposes. */
const rowOpacity = (page: Page, rowName: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`) as HTMLElement | null)?.style.opacity ?? null, rowName);

async function openPlotContextMenu(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await page.waitForTimeout(500);
}

/** Every level below the top one is laid out only while its parent is hovered,
 * so each hover waits for the entry it is meant to reveal; a choice leaf leaves
 * the popup on screen, so it is dismissed afterwards. */
async function clickContextMenuLeaf(page: Page, groups: string[], leaf: string): Promise<void> {
  await openPlotContextMenu(page);
  for (let i = 0; i < groups.length; i++) {
    await page.locator(`.d4-menu-popup [name="${groups[i]}"]`).last().hover();
    const next = i + 1 < groups.length ? groups[i + 1] : leaf;
    await page.locator(`.d4-menu-popup [name="${next}"]`).last().waitFor({timeout: 8000});
  }
  await page.locator(`.d4-menu-popup [name="${leaf}"]`).last().click();
  await page.waitForTimeout(1000);
  await page.keyboard.press('Escape');
  await page.waitForTimeout(500);
}

const X_AXIS_TYPE_MENU = ['div-Properties...', 'div-Properties...---X', 'div-Properties...---X---X-Axis-Type'];
const xAxisTypeLeaf = (choice: 'Linear' | 'Logarithmic') =>
  `div-Properties...---X---X-Axis-Type---${choice}`;

async function wheelZoomIn(page: Page, steps = 1): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  for (let i = 0; i < steps; i++) {
    await page.mouse.wheel(0, -300);
    await page.waitForTimeout(600);
  }
}

/** Reset View is a top-level entry of the plot's own context menu. */
async function resetView(page: Page): Promise<void> {
  await openPlotContextMenu(page);
  await page.locator('.d4-menu-popup [name="div-Reset-View"]').last().click();
  await page.waitForTimeout(800);
  await page.keyboard.press('Escape');
  await page.waitForTimeout(400);
}

const filterCount = (page: Page) =>
  page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount as number);

/** Filtered row count where a NEW value is expected. The filter update trails
 * the gesture, so the settle gate is the count LEAVING the value the caller
 * already knows and only then holding still — agreeing samples on the
 * pre-gesture value are not a settled reading. A count that never moves is
 * returned as it is, so the caller's assertion fails on the real value. */
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

/** Filtered row count where the expectation is that it does NOT move. Here
 * "did not change" cannot be told apart from "has not changed yet", so
 * sampling starts only after a time floor that outlasts the filter update. */
async function settledFilterCountUnchanged(page: Page, floorMs = 2000): Promise<number> {
  await page.waitForTimeout(floorMs);
  let last = -1;
  let stable = 0;
  for (let i = 0; i < 40; i++) {
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

const rowCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);

interface Rect {x: number; y: number; width: number; height: number}

/** The viewer's own viewport rectangle. `props.viewport` is null on this viewer
 * and is not a substitute. */
const viewport = (page: Page): Promise<Rect> => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const vp = sp.viewport;
  return {x: vp.x, y: vp.y, width: vp.width, height: vp.height};
});

/** Two viewport rectangles are the same view when every side agrees to within a
 * small fraction of the reference rectangle — a refit recomputes the rect and is
 * not required to be bit-exact. */
const VIEWPORT_TOLERANCE = 0.02;

function expectSameViewport(actual: Rect, reference: Rect): void {
  expect(Math.abs(actual.width - reference.width)).toBeLessThan(reference.width * VIEWPORT_TOLERANCE);
  expect(Math.abs(actual.height - reference.height)).toBeLessThan(reference.height * VIEWPORT_TOLERANCE);
  expect(Math.abs(actual.x - reference.x)).toBeLessThan(reference.width * VIEWPORT_TOLERANCE);
  expect(Math.abs(actual.y - reference.y)).toBeLessThan(reference.height * VIEWPORT_TOLERANCE);
}

/** Restrict a numeric column to a band through the Filter Panel's range filter,
 * with the band taken from the column's own extremes. */
const applyRangeFilter = (page: Page, column: string, loFraction: number, hiFraction: number) =>
  page.evaluate(async ({col, lo, hi}: {col: string; lo: number; hi: number}) => {
    const tv = grok.shell.tv;
    const stats = tv.dataFrame.col(col).stats;
    const min = stats.min as number;
    const max = stats.max as number;
    const from = min + (max - min) * lo;
    const to = min + (max - min) * hi;
    tv.getFiltersGroup().updateOrAdd({
      type: (window as any).DG.FILTER_TYPE.HISTOGRAM, column: col, min: from, max: to,
    });
    await new Promise((r) => setTimeout(r, 1500));
    return {min, max, from, to};
  }, {col: column, lo: loFraction, hi: hiFraction});

const viewerProps = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {
    x: sp.props.xColumnName, y: sp.props.yColumnName,
    xAxisType: sp.props.xAxisType, yAxisType: sp.props.yAxisType,
    zoomAndFilter: sp.props.zoomAndFilter,
    jitter: sp.props.jitterSize, jitterY: sp.props.jitterSizeY,
  };
});

/** Status-bar filtered-rows indicator. The element is absent from the DOM while
 * the table is unfiltered, so presence is reported alongside the text. */
const filteredIndicator = (page: Page) => page.evaluate(() => {
  const el = document.querySelector('[name="span-filtered"]') as HTMLElement | null;
  return {present: !!el, text: el ? (el.innerText ?? '').trim() : null};
});

async function closeFilterPanel(page: Page): Promise<void> {
  await page.evaluate(() => {
    const f = grok.shell.tv.viewers.find((x: any) => x.type === 'Filters');
    if (f) f.close();
  });
  await page.waitForTimeout(1000);
}

test('Scatter Plot — Zoom and Filter Synchronization', async ({page}: {page: Page}) => {
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
  const fullRowCount = await rowCount(page);
  expect(fullRowCount).toBeGreaterThan(0);

  await softStep('Zoom drives the table filter, reset restores it', async () => {
    const errBefore = errCount();

    // The mode a freshly added viewer carries is read from the settings panel,
    // not from the property object.
    await openSettings(page);
    await revealPropEditor(page, '[name="prop-view-zoom-and-filter"]', 'data');
    if (await propCellText(page, 'prop-view-zoom-and-filter') !== 'filter by zoom')
      await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data', 'filter by zoom');
    expect(await propCellText(page, 'prop-view-zoom-and-filter')).toBe('filter by zoom');

    expect(await settledFilterCountUnchanged(page, 500)).toBe(fullRowCount);

    await wheelZoomIn(page);
    const firstStep = await settledFilterCountAfterChange(page, fullRowCount);
    expect(firstStep).toBeLessThan(fullRowCount);

    await wheelZoomIn(page);
    const secondStep = await settledFilterCountAfterChange(page, firstStep);
    expect(secondStep).toBeLessThan(firstStep);

    await resetView(page);
    expect(await settledFilterCountAfterChange(page, secondStep)).toBe(fullRowCount);

    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Pack and Zoom keeps the zoom filtering resettable', async () => {
    const errBefore = errCount();

    await wheelZoomIn(page, 2);
    const zoomed = await settledFilterCountAfterChange(page, fullRowCount);
    expect(zoomed).toBeLessThan(fullRowCount);

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data',
      'pack and zoom by filter');
    expect(await propCellText(page, 'prop-view-zoom-and-filter')).toBe('pack and zoom by filter');

    // Switching the mode may itself lift the zoom filter, so the reset is
    // graded against whatever the count is when it starts.
    const beforeReset = await filterCount(page);
    await resetView(page);
    const afterReset = beforeReset === fullRowCount
      ? await settledFilterCountUnchanged(page)
      : await settledFilterCountAfterChange(page, beforeReset);
    expect(afterReset).toBe(fullRowCount);
    expect(errCount()).toBe(errBefore);

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data', 'filter by zoom');
    await resetView(page);
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);
    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
  });

  await softStep('External filtering drives the viewport in zoom by filter mode', async () => {
    const errBefore = errCount();

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data',
      'zoom by filter');
    expect((await viewerProps(page)).zoomAndFilter).toBe('zoom by filter');

    const baseline = await viewport(page);
    expect(baseline.width).toBeGreaterThan(0);
    expect(baseline.height).toBeGreaterThan(0);
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);

    await v.openFilterPanel(page);
    let narrowedBand = fullRowCount;
    try {
      // WEIGHT is the column on the X axis, so the surviving rows span a
      // narrower range on that axis and the viewport has to follow.
      const band = await applyRangeFilter(page, 'WEIGHT', 0.4, 0.6);
      expect(band.from).toBeGreaterThan(band.min);
      expect(band.to).toBeLessThan(band.max);

      const narrowed = await settledFilterCountAfterChange(page, fullRowCount);
      narrowedBand = narrowed;
      expect(narrowed).toBeGreaterThan(0);
      expect(narrowed).toBeLessThan(fullRowCount);

      const fitted = await viewport(page);
      expect(fitted.width).toBeGreaterThan(0);
      // Well past the band within which two rectangles count as the same view.
      expect(fitted.width).toBeLessThan(baseline.width * (1 - VIEWPORT_TOLERANCE));
      expect(errCount()).toBe(errBefore);
    } finally {
      await v.resetFilters(page);
    }

    expect(await settledFilterCountAfterChange(page, narrowedBand)).toBe(fullRowCount);
    expectSameViewport(await viewport(page), baseline);

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data',
      'filter by zoom');
    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
    await closeFilterPanel(page);
  });

  await softStep('Filter Panel reset clears the scatter plot\'s contribution', async () => {
    const errBefore = errCount();
    await v.openFilterPanel(page);

    await wheelZoomIn(page, 2);
    const zoomed = await settledFilterCountAfterChange(page, fullRowCount);
    expect(zoomed).toBeLessThan(fullRowCount);

    const reported = await filteredIndicator(page);
    expect(reported.present).toBe(true);
    expect(reported.text).toContain(String(zoomed));

    await page.locator('[name="viewer-Filters"] [name="icon-arrow-rotate-left"]').click();
    expect(await settledFilterCountAfterChange(page, zoomed)).toBe(fullRowCount);

    expect((await filteredIndicator(page)).present).toBe(false);
    expect(errCount()).toBe(errBefore);

    await closeFilterPanel(page);
    await resetView(page);
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);
  });

  await softStep('Axis type switch on a datetime axis keeps the applied filter', async () => {
    const errBefore = errCount();
    await pickOnViewer(page, 'x', 'STARTED');
    expect((await viewerProps(page)).x).toBe('STARTED');

    await v.openFilterPanel(page);
    await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
    const baseline = await settledFilterCountAfterChange(page, fullRowCount);
    expect(baseline).toBeLessThan(fullRowCount);
    expect(baseline).toBeGreaterThan(0);

    // While a datetime column sits on the axis the settings-panel row is greyed
    // out, so the switch is driven from the context menu instead.
    await openSettings(page);
    await revealPropEditor(page, '[name="prop-view-x-axis-type"]', 'x');
    expect(await rowOpacity(page, 'prop-x-axis-type')).toBe('0.5');

    expect((await viewerProps(page)).xAxisType).toBe('linear');
    await clickContextMenuLeaf(page, X_AXIS_TYPE_MENU, xAxisTypeLeaf('Logarithmic'));

    // An unchanged filter over an axis that never switched would grade nothing,
    // so the switch is asserted before the count.
    expect((await viewerProps(page)).xAxisType).toBe('logarithmic');

    expect(await settledFilterCountUnchanged(page)).toBe(baseline);
    expect(errCount()).toBe(errBefore);

    await clickContextMenuLeaf(page, X_AXIS_TYPE_MENU, xAxisTypeLeaf('Linear'));
    expect((await viewerProps(page)).xAxisType).toBe('linear');
    await v.resetFilters(page);
    await closeFilterPanel(page);
    await pickOnViewer(page, 'x', 'WEIGHT');
    await resetView(page);
    const reverted = await viewerProps(page);
    expect(reverted.x).toBe('WEIGHT');
    expect(await settledFilterCountAfterChange(page, baseline)).toBe(fullRowCount);
  });

  await softStep('Filter Out Invalid removes the rows a logarithmic axis cannot draw', async () => {
    const errBefore = errCount();

    // Every numeric column of demog is strictly positive, so the probe pair is
    // built here: an X free of both non-positive and missing values (Filter Out
    // Invalid drops a row undrawable on either axis), and a Y carrying
    // non-positive values, counted from the data.
    const fixture = await page.evaluate(({xName, yName}) => {
      const df = grok.shell.tv.dataFrame;
      for (const n of [xName, yName]) if (df.col(n)) df.columns.remove(n);
      const xc = df.columns.addNewFloat(xName);
      const yc = df.columns.addNewFloat(yName);
      for (let i = 0; i < df.rowCount; i++) {
        xc.set(i, i + 1);
        yc.set(i, i % 4 === 0 ? -(i + 1) : (i + 1));
      }
      let nonPositive = 0;
      let xUndrawable = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (!(yc.get(i) > 0)) nonPositive++;
        if (xc.isNone(i) || !isFinite(xc.get(i))) xUndrawable++;
      }
      return {nonPositive, xUndrawable, rowCount: df.rowCount as number};
    }, {xName: PROBE_X, yName: PROBE_Y});
    await page.waitForTimeout(1500);
    expect(fixture.rowCount).toBe(fullRowCount);
    expect(fixture.xUndrawable).toBe(0);
    expect(fixture.nonPositive).toBeGreaterThan(0);
    expect(fixture.nonPositive).toBeLessThan(fullRowCount);

    try {
      await pickOnViewer(page, 'x', PROBE_X);
      await pickOnViewer(page, 'y', PROBE_Y);
      const axes = await viewerProps(page);
      expect(axes.x).toBe(PROBE_X);
      expect(axes.y).toBe(PROBE_Y);

      await openSettings(page);
      await revealPropEditor(page, '[name="prop-view-filter-out-invalid"]', 'data');
      expect(await propertyRowNames(page)).toContain('prop-filter-out-invalid');

      await setFilterOutInvalid(page, false);
      expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);

      await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'logarithmic');
      expect((await viewerProps(page)).yAxisType).toBe('logarithmic');
      expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);

      await setFilterOutInvalid(page, true);
      expect(await settledFilterCountAfterChange(page, fullRowCount))
        .toBe(fullRowCount - fixture.nonPositive);

      await setFilterOutInvalid(page, false);
      expect(await settledFilterCountAfterChange(page, fullRowCount - fixture.nonPositive))
        .toBe(fullRowCount);

      // Back on linear the same values are drawable, so the setting removes
      // nothing — the drop above came from the axis type.
      await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'linear');
      await setFilterOutInvalid(page, true);
      expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);
      expect(errCount()).toBe(errBefore);
    } finally {
      await setFilterOutInvalid(page, false);
      await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'linear');
      await pickOnViewer(page, 'x', 'WEIGHT');
      await pickOnViewer(page, 'y', 'HEIGHT');
      await page.evaluate((names: string[]) => {
        const df = grok.shell.tv.dataFrame;
        for (const n of names) if (df.col(n)) df.columns.remove(n);
      }, [PROBE_X, PROBE_Y]);
      await page.waitForTimeout(1500);
    }
    const reverted = await viewerProps(page);
    expect(reverted.x).toBe('WEIGHT');
    expect(reverted.y).toBe('HEIGHT');
    expect(reverted.yAxisType).toBe('linear');
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);
  });

  await softStep('Large jitter with a logarithmic axis does not filter rows', async () => {
    const errBefore = errCount();
    await page.evaluate(async (path: string) => {
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise((resolve) => {
        const s = df.onSemanticTypeDetected.subscribe(() => { s.unsubscribe(); resolve(null); });
        setTimeout(resolve, 6000);
      });
    }, spgiPath);
    await page.waitForFunction((col: string) =>
      !!grok.shell.tv?.dataFrame?.col(col), SPGI_X, {timeout: 30000});
    await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
    await waitPlotCanvas(page);

    // The auto-picked columns are set explicitly all the same; the Y column
    // carries strictly positive identifiers, which is what makes a logarithmic
    // Y axis meaningful rather than degenerate.
    await pickOnViewer(page, 'x', SPGI_X);
    await pickOnViewer(page, 'y', SPGI_Y);
    const axes = await viewerProps(page);
    expect(axes.x).toBe(SPGI_X);
    expect(axes.y).toBe(SPGI_Y);
    const yMin = await page.evaluate((c: string) =>
      grok.shell.tv.dataFrame.col(c).stats.min as number, SPGI_Y);
    expect(yMin).toBeGreaterThan(0);

    const spgiRowCount = await rowCount(page);
    expect(spgiRowCount).toBeGreaterThan(0);
    expect(await settledFilterCountUnchanged(page, 500)).toBe(spgiRowCount);

    await setSliderProp(page, 'prop-jitter-size', 'marker', JITTER_X);
    await setSliderProp(page, 'prop-jitter-size-y', 'marker', JITTER_Y);
    const jittered = await viewerProps(page);
    expect(jittered.jitter).toBe(JITTER_X);
    expect(jittered.jitterY).toBe(JITTER_Y);

    await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'logarithmic');
    expect((await viewerProps(page)).yAxisType).toBe('logarithmic');

    expect(await settledFilterCountUnchanged(page)).toBe(spgiRowCount);
    expect(errCount()).toBe(errBefore);

    await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y', 'linear');
    await setSliderProp(page, 'prop-jitter-size', 'marker', 0);
    await setSliderProp(page, 'prop-jitter-size-y', 'marker', 0);
    const back = await viewerProps(page);
    expect(back.yAxisType).toBe('linear');
    expect(back.jitter).toBe(0);
    expect(back.jitterY).toBe(0);
    expect(await settledFilterCountUnchanged(page, 500)).toBe(spgiRowCount);
    await page.evaluate(() => grok.shell.v.close());
    await page.waitForTimeout(1500);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

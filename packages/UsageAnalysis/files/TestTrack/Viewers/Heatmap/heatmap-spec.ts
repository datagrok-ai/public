/* ---
realizes: [viewers.heat-map, viewers.filters.histogram]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {knownOpenBug} from '../../helpers/known-open-bug';

test.use(specTestOptions);

const VIEWER_NAME = 'Heat-map';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const VIEWER_TYPE = 'Heat map';
const datasetPath = 'System:DemoFiles/demog.csv';

const CONTENT_CANVAS = 'canvas[name="canvas"]';

const shownValue = (page: Page, prop: string) => v.propertyGridValue(page, prop);
const category = (page: Page, cat: string, probe: string) =>
  v.ensurePropertyCategory(page, VIEWER_NAME, cat, probe);
const snapshot = (page: Page) => v.snapshotCanvasColors(page, VIEWER_TYPE, CONTENT_CANVAS);
const diff = (page: Page) => v.diffCanvasColors(page, VIEWER_TYPE, CONTENT_CANVAS);
const settle = (page: Page) => v.waitForCanvasQuiet(page, VIEWER_TYPE, {canvasSelector: CONTENT_CANVAS});
const repaint = (page: Page, minDelta: number, timeoutMs?: number) =>
  v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta, canvasSelector: CONTENT_CANVAS, timeoutMs});

async function sliderSpan(page: Page, slider: 'x' | 'y'): Promise<number> {
  return page.evaluate(({sel, s}) => {
    const root = document.querySelector(sel) as HTMLElement;
    const svg = root.querySelector(`[name="${s}-slider"]`) as SVGElement | null;
    const read = (name: string) => {
      const el = svg?.querySelector(`[name="${name}"]`) as SVGCircleElement | null;
      return el ? {x: Number(el.getAttribute('cx')), y: Number(el.getAttribute('cy'))} : null;
    };
    const min = read('min-handle');
    const max = read('max-handle');
    if (!min || !max) return -1;
    return s === 'x' ? max.x - min.x : max.y - min.y;
  }, {sel: VIEWER, s: slider});
}

test('Heat map', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Add Heat map from the Viewers toolbox', async () => {
    await page.locator('[name="icon-heat-map"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});
    await expect.poll(async () =>
      (await v.countCanvasPixels(page, VIEWER_TYPE, {canvasSelector: CONTENT_CANVAS})).total,
    {timeout: 30_000}).toBeGreaterThan(1000);

    await v.openViewerProperties(page, VIEWER_NAME);
    await category(page, 'misc', 'is-heatmap');
    expect(await shownValue(page, 'is-heatmap')).toBe('true');
  });

  await softStep('Heatmap Colors off falls back to plain cells', async () => {

    await settle(page);
    await category(page, 'misc', 'heatmap-colors');
    expect(await shownValue(page, 'heatmap-colors')).toBe('true');

    await snapshot(page);
    expect(await v.togglePropertyGridCheckbox(page, 'heatmap-colors')).toBe(false);

    await knownOpenBug('GROK-20619', async () => {
      await repaint(page, 1000, 4000);
    });

    expect(await v.togglePropertyGridCheckbox(page, 'heatmap-colors')).toBe(true);
  });

  await softStep('Global Color Scaling recolours the cells', async () => {
    await settle(page);
    await category(page, 'misc', 'global-color-scaling');
    const was = await shownValue(page, 'global-color-scaling');
    await snapshot(page);
    expect(String(await v.togglePropertyGridCheckbox(page, 'global-color-scaling'))).not.toBe(was);
    await repaint(page, 500);

    expect(String(await v.togglePropertyGridCheckbox(page, 'global-color-scaling'))).toBe(was);
  });

  await softStep('Col Labels Orientation rotates the header', async () => {
    await category(page, 'style', 'col-labels-orientation');
    expect(await shownValue(page, 'col-labels-orientation')).toBe('Auto');

    for (const orientation of ['Vert', 'Horz', 'Auto']) {
      await v.selectPropertyGridChoice(page, 'col-labels-orientation', orientation);
      await v.waitForPropertyValue(page, 'col-labels-orientation', orientation);
    }
  });

  await softStep('Max Heatmap Columns limits the columns on screen', async () => {
    await settle(page);
    await category(page, 'columns', 'max-heatmap-columns');
    expect(await shownValue(page, 'max-heatmap-columns')).toBe('100');

    await snapshot(page);
    await v.setPropertyGridValue(page, 'max-heatmap-columns', '3');
    expect(await shownValue(page, 'max-heatmap-columns')).toBe('3');
    await repaint(page, 1000);

    await v.setPropertyGridValue(page, 'max-heatmap-columns', '100');
    await v.waitForPropertyValue(page, 'max-heatmap-columns', '100');
  });

  await softStep('Show Heatmap Scrollbars hides the range sliders', async () => {
    await category(page, 'misc', 'show-heatmap-scrollbars');
    expect(await v.togglePropertyGridCheckbox(page, 'show-heatmap-scrollbars')).toBe(false);
    await expect(page.locator(`${VIEWER} [name="x-slider"]`).first()).toBeHidden();

    expect(await v.togglePropertyGridCheckbox(page, 'show-heatmap-scrollbars')).toBe(true);
    await expect(page.locator(`${VIEWER} [name="x-slider"]`).first()).toBeVisible();
  });

  await softStep('Is Heatmap off turns the viewer into a plain grid', async () => {
    await settle(page);
    await category(page, 'misc', 'is-heatmap');
    await snapshot(page);
    expect(await v.togglePropertyGridCheckbox(page, 'is-heatmap')).toBe(false);
    await repaint(page, 1000);

    await settle(page);
    await category(page, 'misc', 'is-heatmap');
    await snapshot(page);
    expect(await v.togglePropertyGridCheckbox(page, 'is-heatmap')).toBe(true);
    await repaint(page, 1000);
  });

  await softStep('Row Height is not offered in heatmap mode', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);

    await knownOpenBug('GROK-20619', async () => {
      expect(await page.locator('.property-grid tr[name="prop-row-height"]').count()).toBe(0);
    });
  });

  await softStep('Clicking a cell makes its row current', async () => {
    const box = (await page.locator(VIEWER).boundingBox())!;
    const currentRow = () =>
      page.evaluate(() => (window as any).grok.shell.t.currentRowIdx as number);
    const before = await currentRow();
    await page.mouse.click(box.x + box.width * 0.5, box.y + box.height * 0.5);
    await expect.poll(currentRow, {timeout: 8000}).not.toBe(before);
    expect(await currentRow()).toBeGreaterThanOrEqual(0);
  });

  await softStep('Alt+drag zooms into an area and the slider follows', async () => {
    const box = (await page.locator(VIEWER).boundingBox())!;
    const before = await sliderSpan(page, 'y');
    expect(before).toBeGreaterThan(0);
    await settle(page);
    await snapshot(page);

    await page.keyboard.down('Alt');
    await page.mouse.move(box.x + box.width * 0.3, box.y + box.height * 0.3);
    await page.mouse.down();
    await page.mouse.move(box.x + box.width * 0.6, box.y + box.height * 0.6, {steps: 20});
    await page.mouse.up();
    await page.keyboard.up('Alt');

    await expect.poll(() => sliderSpan(page, 'y'), {timeout: 10_000}).toBeLessThan(before);
    await repaint(page, 500);

    await page.locator(`${VIEWER} [name="y-slider"]`).first().dblclick();
    await expect.poll(() => sliderSpan(page, 'y'), {timeout: 10_000}).toBeGreaterThan(before);
  });

  await softStep('Filtering the table redraws the heat map', async () => {
    await settle(page);
    await snapshot(page);
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['M']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeLessThan(total);
    await repaint(page, 500);

    await v.resetFilters(page);
  });

  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});

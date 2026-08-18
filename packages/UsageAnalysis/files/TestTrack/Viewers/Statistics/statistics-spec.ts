/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = 'Statistics';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const VIEWER_TYPE = 'Statistics';
const datasetPath = 'System:DemoFiles/demog.csv';

const CONTENT_CANVAS = 'canvas[name="canvas"]';
const snapshot = (page: Page) => v.snapshotCanvasColors(page, VIEWER_TYPE, CONTENT_CANVAS);
const diff = (page: Page) => v.diffCanvasColors(page, VIEWER_TYPE, CONTENT_CANVAS);
const content = (page: Page) =>
  v.countCanvasPixels(page, VIEWER_TYPE, {canvasSelector: CONTENT_CANVAS});

async function menuItemState(page: Page, item: string): Promise<'on' | 'off'> {
  const cls = await page.locator(`[name="${item}"]:visible`).last()
    .locator('.d4-menu-item-check i').getAttribute('class') ?? '';
  return /fa-check|fa-dot-circle/.test(cls) ? 'on' : 'off';
}

test('Statistics viewer', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Add, close and re-add the viewer', async () => {
    await page.locator('[name="icon-statistics"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});
    const drawn = () => expect.poll(async () => (await content(page)).total,
      {timeout: 30_000}).toBeGreaterThan(1000);
    await drawn();

    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);

    await page.locator('[name="icon-statistics"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});
    await drawn();
  });

  await softStep('Statistics > sum adds a column to the viewer', async () => {
    const box = (await page.locator(VIEWER).first().boundingBox())!;
    await snapshot(page);

    await page.mouse.click(box.x + box.width * 0.5, box.y + box.height * 0.4, {button: 'right'});
    await page.locator('.d4-menu-popup:visible').first().waitFor({timeout: 10_000});
    await page.locator('[name="div-Statistics"]:visible').last().hover();
    await page.locator('[name="div-Statistics---sum"]:visible').last().waitFor({timeout: 10_000});

    expect(await menuItemState(page, 'div-Statistics---sum')).toBe('off');
    await page.locator('[name="div-Statistics---sum"]:visible').last().click();

    await v.waitForCanvasChange(page, VIEWER_TYPE,
      {minDelta: 500, canvasSelector: CONTENT_CANVAS});
  });

  await softStep('Row Source Filtered follows the table filter', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);
    await v.ensurePropertyCategory(page, VIEWER_NAME, 'data', 'row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'Filtered');
    expect(await v.propertyGridValue(page, 'row-source')).toBe('Filtered');

    await v.waitForCanvasQuiet(page, VIEWER_TYPE, {canvasSelector: CONTENT_CANVAS});
    await snapshot(page);
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['M']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeGreaterThan(0);
    expect(filteredCount).toBeLessThan(total);

    await v.waitForCanvasChange(page, VIEWER_TYPE,
      {minDelta: 500, canvasSelector: CONTENT_CANVAS});

    await v.resetFilters(page);
  });

  await softStep('Row Source Selected follows the selection', async () => {
    await v.ensurePropertyCategory(page, VIEWER_NAME, 'data', 'row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'Selected');
    expect(await v.propertyGridValue(page, 'row-source')).toBe('Selected');

    await v.waitForCanvasQuiet(page, VIEWER_TYPE, {canvasSelector: CONTENT_CANVAS});
    await snapshot(page);
    await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      df.selection.init((i: number) => i % 2 === 0);
    });
    await v.waitForCanvasChange(page, VIEWER_TYPE,
      {minDelta: 500, canvasSelector: CONTENT_CANVAS});

    await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
    await v.ensurePropertyCategory(page, VIEWER_NAME, 'data', 'row-source');
    await v.selectPropertyGridChoice(page, 'row-source', 'All');
    expect(await v.propertyGridValue(page, 'row-source')).toBe('All');
  });

  await softStep('The Columns property lists every column of the table', async () => {
    await v.ensurePropertyCategory(page, VIEWER_NAME, 'data', 'column-names');
    const columns = await v.propertyGridValue(page, 'column-names');
    const total = await page.evaluate(() => (window as any).grok.shell.t.columns.length);
    expect(columns).toContain(`/ ${total}`);
  });

  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});

/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = 'Tree-map';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const VIEWER_TYPE = 'Tree map';
const datasetPath = 'System:DemoFiles/demog.csv';

const splitSelect = (page: Page, level: number) =>
  page.locator(`${VIEWER} select.d4-column-selector-tree-map`).nth(level);

const shownValue = (page: Page, prop: string) => v.propertyGridValue(page, prop);
const category = (page: Page, cat: string, probe: string) =>
  v.ensurePropertyCategory(page, VIEWER_NAME, cat, probe);

async function categoryCounts(page: Page, column: string, filteredOnly = false): Promise<Record<string, number>> {
  return page.evaluate(({col, filtered}) => {
    const df = (window as any).grok.shell.t;
    const c = df.col(col);
    const counts: Record<string, number> = {};
    for (let i = 0; i < df.rowCount; i++) {
      if (filtered && !df.filter.get(i)) continue;
      const value = String(c.get(i));
      counts[value] = (counts[value] ?? 0) + 1;
    }
    return counts;
  }, {col: column, filtered: filteredOnly});
}

async function hoverLeaf(page: Page): Promise<{name: string; rows: number; x: number; y: number}> {
  const box = (await page.locator(VIEWER).boundingBox())!;
  for (const [dx, dy] of [[0.3, 0.3], [0.7, 0.3], [0.5, 0.7], [0.2, 0.6], [0.85, 0.7]]) {
    const x = box.x + box.width * dx;
    const y = box.y + box.height * dy;
    await page.mouse.move(x - 20, y - 20);
    await page.mouse.move(x, y, {steps: 5});
    await page.waitForTimeout(400);
    const text = await page.evaluate(() => {
      const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
      return t && getComputedStyle(t).display !== 'none' ? t.innerText : '';
    });
    const match = text.match(/^(.+?)\s*\n\s*(\d+)\s+rows/);
    if (match) return {name: match[1].trim(), rows: Number(match[2]), x, y};
  }
  throw new Error('no leaf tooltip appeared anywhere on the tree map');
}

test('Tree map', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Add Tree map from the Viewers toolbox', async () => {
    await page.locator('[name="icon-tree-map"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});
    await expect.poll(async () => (await v.countCanvasPixels(page, VIEWER_TYPE)).total,
      {timeout: 30_000}).toBeGreaterThan(1000);

    expect(await splitSelect(page, 0).inputValue()).not.toBe('');

    expect(await splitSelect(page, 1).inputValue()).toBe('');
  });

  await softStep('Split by RACE draws one leaf per race', async () => {
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await splitSelect(page, 0).selectOption('RACE');

    expect(await splitSelect(page, 0).inputValue()).toBe('RACE');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500});

    const leaf = await hoverLeaf(page);
    const counts = await categoryCounts(page, 'RACE');
    expect(Object.keys(counts)).toContain(leaf.name);
    expect(leaf.rows).toBe(counts[leaf.name]);
  });

  await softStep('Adding SEX as a second level nests the leaves', async () => {
    const levels = page.locator(`${VIEWER} select.d4-column-selector-tree-map`);
    const levelsBefore = await levels.count();
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    await splitSelect(page, 1).selectOption('SEX');

    await expect(levels).toHaveCount(levelsBefore + 1);
    expect(await splitSelect(page, 1).inputValue()).toBe('SEX');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500});

    const leaf = await hoverLeaf(page);
    const raceCounts = await categoryCounts(page, 'RACE');
    expect(Math.max(...Object.values(raceCounts))).toBeGreaterThan(leaf.rows);

    await splitSelect(page, 1).selectOption('');
    await expect(levels).toHaveCount(levelsBefore);
  });

  await softStep('Colour by AGE and switch the aggregation', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);
    await category(page, 'data', 'color');
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'AGE', viewerType: VIEWER_TYPE,
      propName: 'colorColumnName', scopeSelector: '.property-grid',
    });
    expect(await shownValue(page, 'color')).toBe('AGE');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 300});

    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.selectPropertyGridChoice(page, 'color-aggr-type', 'max');
    expect(await shownValue(page, 'color-aggr-type')).toBe('max');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 300});
  });

  await softStep('Size by WEIGHT rescales the rectangles', async () => {
    await category(page, 'data', 'size');
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    await v.pickColumnViaSelectorTrusted(page, {
      role: 'size', columnName: 'WEIGHT', viewerType: VIEWER_TYPE,
      propName: 'sizeColumnName', scopeSelector: '.property-grid',
    });
    expect(await shownValue(page, 'size')).toBe('WEIGHT');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 300});

    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.selectPropertyGridChoice(page, 'size-aggr-type', 'max');
    expect(await shownValue(page, 'size-aggr-type')).toBe('max');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 300});
  });

  await softStep('Clicking a leaf selects exactly that group', async () => {
    const selected = () =>
      page.evaluate(() => (window as any).grok.shell.t.selection.trueCount as number);
    await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
    await expect.poll(selected, {timeout: 5000}).toBe(0);

    const leaf = await hoverLeaf(page);
    await page.mouse.click(leaf.x, leaf.y);
    await expect.poll(selected, {timeout: 8000}).toBe(leaf.rows);

    await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
  });

  await softStep('Filtering the table reshapes the map', async () => {
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['M']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeLessThan(total);
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 300});

    const leaf = await hoverLeaf(page);
    const counts = await categoryCounts(page, 'RACE', true);
    expect(leaf.rows).toBe(counts[leaf.name]);

    await v.resetFilters(page);
  });

  await softStep('Show Column Selection Panel hides the on-viewer selectors', async () => {
    await category(page, 'misc', 'show-column-selection-panel');

    await v.setPropertyGridCheckbox(page, 'show-column-selection-panel', false);
    await expect(splitSelect(page, 0)).toBeHidden();

    await v.setPropertyGridCheckbox(page, 'show-column-selection-panel', true);
    await expect(splitSelect(page, 0)).toBeVisible();
  });

  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});

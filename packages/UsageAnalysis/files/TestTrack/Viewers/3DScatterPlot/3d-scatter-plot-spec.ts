/* ---
realizes: []
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = '3d-scatter-plot';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const VIEWER_TYPE = '3d scatter plot';
const datasetPath = 'System:DemoFiles/demog.csv';

const signature = (page: Page) => v.viewerSignature(page, VIEWER_NAME);
const repaints = (page: Page, before: string) => v.waitForViewerRepaint(page, VIEWER_NAME, before);
const shownValue = (page: Page, prop: string) => v.propertyGridValue(page, prop);
const category = (page: Page, cat: string, probe: string) =>
  v.ensurePropertyCategory(page, VIEWER_NAME, cat, probe);

async function selectorText(page: Page, role: string): Promise<string> {
  return (await page.locator(`${VIEWER} [name="div-column-combobox-${role}"]`).first().innerText())
    .replace(/\s+/g, ' ').trim();
}

async function plotCentre(page: Page): Promise<{x: number; y: number; box: any}> {
  const box = (await page.locator(VIEWER).boundingBox())!;
  return {x: box.x + box.width / 2, y: box.y + box.height / 2, box};
}

test('3D scatter plot', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Add 3D scatter plot from the Viewers toolbox', async () => {
    await page.locator('[name="icon-3d-scatter-plot"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});

    await expect.poll(() => selectorText(page, 'x'), {timeout: 30_000}).toBe('X: AGE');
    expect(await selectorText(page, 'y')).toBe('Y: HEIGHT');
    expect(await selectorText(page, 'z')).toBe('Z: WEIGHT');
  });

  await softStep('Reassign X and Z with the on-viewer selectors', async () => {
    const before = await signature(page);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'x', columnName: 'WEIGHT', viewerType: VIEWER_TYPE, propName: 'xColumnName',
    });
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'z', columnName: 'AGE', viewerType: VIEWER_TYPE, propName: 'zColumnName',
    });
    expect(await selectorText(page, 'x')).toBe('X: WEIGHT');
    expect(await selectorText(page, 'z')).toBe('Z: AGE');
    await repaints(page, before);

    await v.pickColumnViaSelectorTrusted(page, {
      role: 'x', columnName: 'AGE', viewerType: VIEWER_TYPE, propName: 'xColumnName',
    });
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'z', columnName: 'WEIGHT', viewerType: VIEWER_TYPE, propName: 'zColumnName',
    });
  });

  await softStep('Color by SEX shows a categorical legend', async () => {
    const before = await signature(page);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'SEX', viewerType: VIEWER_TYPE, propName: 'colorColumnName',
    });

    await expect.poll(async () => (await v.readLegend(page, VIEWER_TYPE)).labels.sort(),
      {timeout: 10_000}).toEqual(['F', 'M']);
    expect((await v.readLegend(page, VIEWER_TYPE)).legendRendered).toBe(true);
    await repaints(page, before);
  });

  await softStep('Color by AGE switches the legend to a gradient', async () => {
    const before = await signature(page);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'AGE', viewerType: VIEWER_TYPE, propName: 'colorColumnName',
    });

    await expect.poll(async () => (await v.readLegend(page, VIEWER_TYPE)).labels,
      {timeout: 10_000}).not.toEqual(['F', 'M']);
    await repaints(page, before);
  });

  await softStep('Marker type redraws the markers', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);
    await category(page, 'marker', 'marker-type');

    const shapes: Record<string, string> = {};
    for (const shape of ['box', 'sphere', 'cylinder']) {
      await v.selectPropertyGridChoice(page, 'marker-type', shape);
      expect(await shownValue(page, 'marker-type')).toBe(shape);
      shapes[shape] = await signature(page);
    }
    expect(shapes['box']).not.toBe(shapes['sphere']);
    expect(shapes['cylinder']).not.toBe(shapes['sphere']);
  });

  await softStep('Marker opacity redraws the markers', async () => {
    await category(page, 'marker', 'marker-opacity');
    const before = await signature(page);
    await v.setPropertyGridValue(page, 'marker-opacity', '25');
    expect(await shownValue(page, 'marker-opacity')).toBe('25');
    await repaints(page, before);
    await v.setPropertyGridValue(page, 'marker-opacity', '100');
  });

  await softStep('Show Axes hides and restores the axes', async () => {
    await category(page, 'misc', 'show-axes');
    const before = await signature(page);
    expect(await v.togglePropertyGridCheckbox(page, 'show-axes')).toBe(false);
    const withoutAxes = await repaints(page, before);

    expect(await v.togglePropertyGridCheckbox(page, 'show-axes')).toBe(true);
    await repaints(page, withoutAxes);
  });

  await softStep('X axis type switches to logarithmic', async () => {
    await category(page, 'x-axis', 'x-axis-type');
    const before = await signature(page);
    await v.selectPropertyGridChoice(page, 'x-axis-type', 'logarithmic');
    expect(await shownValue(page, 'x-axis-type')).toBe('logarithmic');
    await repaints(page, before);
    await v.selectPropertyGridChoice(page, 'x-axis-type', 'linear');
    expect(await shownValue(page, 'x-axis-type')).toBe('linear');
  });

  await softStep('Drag rotates the scene and Reset View restores it', async () => {
    const {x, y, box} = await plotCentre(page);
    const before = await signature(page);

    await page.mouse.move(x, y);
    await page.mouse.down();
    await page.mouse.move(box.x + box.width * 0.75, box.y + box.height * 0.3, {steps: 20});
    await page.mouse.up();
    const rotated = await repaints(page, before);

    await page.mouse.click(x, y, {button: 'right'});
    await page.locator('.d4-menu-popup').first().waitFor({timeout: 5000});
    await page.locator('.d4-menu-item')
      .filter({has: page.locator('.d4-menu-item-label', {hasText: /^Reset View$/})}).first().click();
    await expect(page.locator('.d4-menu-popup')).toHaveCount(0);
    await repaints(page, rotated);
  });

  await softStep('Mouse wheel zooms the scene', async () => {
    const {x, y} = await plotCentre(page);
    await page.mouse.move(x, y);
    const before = await signature(page);

    await page.mouse.wheel(0, -600);
    const zoomedIn = await repaints(page, before);

    await page.mouse.wheel(0, 600);
    await repaints(page, zoomedIn);
  });

  await softStep('Click makes a row current, Shift+click selects it', async () => {
    const {x, y, box} = await plotCentre(page);
    const startRow = await page.evaluate(() => (window as any).grok.shell.t.currentRowIdx);

    const currentRow = () =>
      page.evaluate(() => (window as any).grok.shell.t.currentRowIdx as number);
    let current = startRow;
    for (const [dx, dy] of [[0.5, 0.5], [0.45, 0.55], [0.55, 0.45], [0.5, 0.6]]) {
      await page.mouse.click(box.x + box.width * dx, box.y + box.height * dy);

      current = await v.pollValue(currentRow, (n) => n !== startRow && n >= 0, 1500);
      if (current !== startRow && current >= 0) break;
    }
    expect(current).toBeGreaterThanOrEqual(0);
    expect(current).not.toBe(startRow);

    const selected = () =>
      page.evaluate(() => (window as any).grok.shell.t.selection.trueCount as number);
    const selectedBefore = await selected();
    await page.keyboard.down('Shift');
    await page.mouse.click(x, y);
    await page.keyboard.up('Shift');
    await expect.poll(selected, {timeout: 8000}).toBeGreaterThan(selectedBefore);

    await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
  });

  await softStep('Show Filtered Out Points repaints the filtered-away rows', async () => {
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['F']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeLessThan(total);

    await category(page, 'data', 'show-filtered-out-points');
    const filteredOnly = await signature(page);
    expect(await v.togglePropertyGridCheckbox(page, 'show-filtered-out-points')).toBe(true);
    await repaints(page, filteredOnly);

    expect(await v.togglePropertyGridCheckbox(page, 'show-filtered-out-points')).toBe(false);
    await v.resetFilters(page);
  });

  await softStep('Hovering a bar chart bin highlights the matching 3D points', async () => {
    await page.locator('[name="icon-bar-chart"]').first().click();
    await page.locator('[name="viewer-Bar-chart"]').first().waitFor({timeout: 30_000});
    await expect.poll(async () =>
      (await v.countCanvasPixels(page, 'Bar chart')).total, {timeout: 30_000}).toBeGreaterThan(1000);

    const bar = (await page.locator('[name="viewer-Bar-chart"]').boundingBox())!;
    const idle = await signature(page);
    await page.mouse.move(bar.x + bar.width * 0.3, bar.y + bar.height * 0.6);
    await repaints(page, idle);

    await page.mouse.move(bar.x + bar.width * 0.9, bar.y + bar.height * 0.05);
    await category(page, 'selection', 'show-mouse-over-row-group');
    expect(await v.togglePropertyGridCheckbox(page, 'show-mouse-over-row-group')).toBe(false);
    expect(await shownValue(page, 'show-mouse-over-row-group')).toBe('false');

    await category(page, 'selection', 'show-mouse-over-row-group');
    expect(await v.togglePropertyGridCheckbox(page, 'show-mouse-over-row-group')).toBe(true);
    await v.clickViewerTitlebarIcon(page, 'Bar-chart', 'Close');
    await expect(page.locator('[name="viewer-Bar-chart"]')).toHaveCount(0);
  });

  await softStep('Legend position moves the legend', async () => {
    // an earlier step left Color = AGE (numeric), which renders NO legend element at all —
    // there is nothing to move, so the repaint could never happen. Colour by a CATEGORICAL
    // column first so a legend exists to reposition.
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'SEX', viewerType: VIEWER_TYPE, propName: 'colorColumnName',
    });
    await expect.poll(async () => (await v.readLegend(page, VIEWER_TYPE)).legendRendered,
      {timeout: 10_000}).toBe(true);

    await category(page, 'legend', 'legend-position');
    await v.selectPropertyGridChoice(page, 'legend-visibility', 'Always');
    const before = await signature(page);
    await v.selectPropertyGridChoice(page, 'legend-position', 'Left');
    expect(await shownValue(page, 'legend-position')).toBe('Left');
    await repaints(page, before);
    await v.selectPropertyGridChoice(page, 'legend-position', 'Auto');
  });

  await v.cleanupShell(page);

  v.finishSpec();
});

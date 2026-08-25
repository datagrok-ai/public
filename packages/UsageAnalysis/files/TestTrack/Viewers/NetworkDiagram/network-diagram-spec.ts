/* ---
realizes: [viewers.network-diagram, entities.viewer.action.close-viewer]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = 'Network-diagram';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const VIEWER_TYPE = 'Network diagram';
const datasetPath = 'System:DemoFiles/demog.csv';

const shownValue = (page: Page, prop: string, cat?: string) => v.propertyGridValue(page, prop, cat);
const category = (page: Page, cat: string, probe: string) =>
  v.ensurePropertyCategory(page, VIEWER_NAME, cat, probe);

const selectionCount = (page: Page) =>
  page.evaluate(() => (window as any).grok.shell.t.selection.trueCount as number);

const clearSelection = async (page: Page) => {
  await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
  await expect.poll(() => selectionCount(page), {timeout: 5000}).toBe(0);
};

async function selectionSettles(page: Page): Promise<number> {
  let last = -1;
  for (let i = 0; i < 8; i++) {
    await page.waitForTimeout(150);
    const now = await selectionCount(page);
    if (now === last && now > 0) return now;
    last = now;
  }
  return last;
}

async function selectorText(page: Page, role: string): Promise<string> {
  return (await page.locator(`${VIEWER} [name="div-column-combobox-${role}"]`).first().innerText())
    .replace(/\s+/g, ' ').trim();
}

async function nodePositions(page: Page): Promise<{x: number; y: number; n: number}[]> {
  return page.evaluate((sel) => {
    const root = document.querySelector(sel) as HTMLElement;
    const cv = root.querySelector('canvas') as HTMLCanvasElement;
    const r = cv.getBoundingClientRect();
    const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
    const sx = cv.width / r.width;
    const sy = cv.height / r.height;
    const buckets = new Map<string, {n: number; x: number; y: number}>();
    for (let y = 0; y < cv.height; y += 3) {
      for (let x = 0; x < cv.width; x += 3) {
        const i = (y * cv.width + x) * 4;
        const rr = img[i], gg = img[i + 1], bb = img[i + 2];
        if (img[i + 3] === 0) continue;
        const mx = Math.max(rr, gg, bb), mn = Math.min(rr, gg, bb);
        if ((mx > 245 && mn > 245) || mx - mn < 25) continue;
        const key = `${Math.floor(x / 40)}:${Math.floor(y / 40)}`;
        const b = buckets.get(key) ?? {n: 0, x: 0, y: 0};
        b.n++; b.x += x; b.y += y;
        buckets.set(key, b);
      }
    }
    return [...buckets.values()].filter((b) => b.n > 20)
      .sort((a, b) => b.n - a.n).slice(0, 8)
      .map((b) => ({n: b.n, x: r.x + (b.x / b.n) / sx, y: r.y + (b.y / b.n) / sy}));
  }, VIEWER);
}

test('Network diagram', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Add Network diagram from the Viewers toolbox', async () => {
    await page.locator('[name="icon-network-diagram"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});

    expect(await selectorText(page, 'node1')).toBe('SEX');
    expect(await selectorText(page, 'node2')).toBe('CONTROL');

    await expect.poll(async () => (await v.countCanvasPixels(page, VIEWER_TYPE)).total,
      {timeout: 30_000}).toBeGreaterThan(1000);
  });

  await softStep('Suspend simulation freezes the layout', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);
    await category(page, 'misc', 'suspend-simulation');
    expect(await v.togglePropertyGridCheckbox(page, 'suspend-simulation', 'misc')).toBe(true);

    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
  });

  await softStep('Switch Node 1 to RACE', async () => {
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'node1', columnName: 'RACE', viewerType: VIEWER_TYPE, propName: 'node1ColumnName',
    });
    expect(await selectorText(page, 'node1')).toBe('RACE');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500});
  });

  await softStep('Colour nodes by SEX and size them by AGE', async () => {
    await category(page, 'data', 'node1-color');
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    await v.pickColumnViaSelectorTrusted(page, {
      role: 'node1--color', columnName: 'SEX', viewerType: VIEWER_TYPE,
      propName: 'node1ColorColumnName', scopeSelector: '.property-grid',
    });
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'node1--size', columnName: 'AGE', viewerType: VIEWER_TYPE,
      propName: 'node1SizeColumnName', scopeSelector: '.property-grid',
    });

    expect(await shownValue(page, 'node1-color')).toBe('SEX');
    expect(await shownValue(page, 'node1-size')).toBe('AGE');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500});
  });

  await softStep('Colour edges by AGE', async () => {
    await category(page, 'data', 'edge-color');
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'edge--color', columnName: 'AGE', viewerType: VIEWER_TYPE,
      propName: 'edgeColorColumnName', scopeSelector: '.property-grid',
    });
    expect(await shownValue(page, 'edge-color')).toBe('AGE');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 300});
  });

  await softStep('Clicking a node selects the rows behind it', async () => {
    await clearSelection(page);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    const nodes = await nodePositions(page);
    expect(nodes.length).toBeGreaterThan(0);

    let selected = 0;
    for (const node of nodes.slice(0, 6)) {
      await page.mouse.click(node.x, node.y);
      selected = await selectionSettles(page);
      if (selected > 0) break;
    }
    expect(selected).toBeGreaterThan(0);
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});
  });

  await softStep('Click-selection switches off stop the clicks from selecting', async () => {
    await category(page, 'misc', 'select-rows-on-click');
    expect(await v.togglePropertyGridCheckbox(page, 'select-rows-on-click', 'misc')).toBe(false);
    expect(await v.togglePropertyGridCheckbox(page, 'select-edges-on-click', 'misc')).toBe(false);
    await clearSelection(page);

    const nodes = await nodePositions(page);
    for (const node of nodes.slice(0, 4)) {
      await page.mouse.click(node.x, node.y);
      await selectionSettles(page);
    }
    expect(await selectionCount(page)).toBe(0);

    await category(page, 'misc', 'select-rows-on-click');
    expect(await v.togglePropertyGridCheckbox(page, 'select-rows-on-click', 'misc')).toBe(true);
    expect(await v.togglePropertyGridCheckbox(page, 'select-edges-on-click', 'misc')).toBe(true);
  });

  await softStep('Show Column Selectors hides the on-viewer selectors', async () => {
    await category(page, 'misc', 'show-column-selectors');
    expect(await v.togglePropertyGridCheckbox(page, 'show-column-selectors', 'misc')).toBe(false);
    await expect(page.locator(`${VIEWER} [name="div-column-combobox-node1"]`)).toBeHidden();

    await category(page, 'misc', 'show-column-selectors');
    expect(await v.togglePropertyGridCheckbox(page, 'show-column-selectors', 'misc')).toBe(true);
    await expect(page.locator(`${VIEWER} [name="div-column-combobox-node1"]`).first()).toBeVisible();
  });

  await softStep('Show Arrows draws directions on the edges', async () => {
    await category(page, 'misc', 'show-arrows');
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.selectPropertyGridChoice(page, 'show-arrows', 'to', 'misc');
    await v.waitForPropertyValue(page, 'show-arrows', 'to', 'misc');

    await v.waitForCanvasChange(page, VIEWER_TYPE, {timeoutMs: 10_000});

    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'node1', columnName: 'SEX', viewerType: VIEWER_TYPE, propName: 'node1ColumnName',
    });
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500, timeoutMs: 20_000});
    expect(await shownValue(page, 'show-arrows', 'misc')).toBe('to');
  });

  await softStep('Show Filtered Out Nodes brings the filtered-away nodes back', async () => {

    if (await selectorText(page, 'node1') !== 'SEX')
      await v.pickColumnViaSelectorTrusted(page, {
        role: 'node1', columnName: 'SEX', viewerType: VIEWER_TYPE, propName: 'node1ColumnName',
      });

    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['F']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeGreaterThan(0);
    expect(filteredCount).toBeLessThan(total);

    await category(page, 'misc', 'show-filtered-out-nodes');

    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    const hidden = (await v.countCanvasPixels(page, VIEWER_TYPE)).total;
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    expect(await v.togglePropertyGridCheckbox(page, 'show-filtered-out-nodes', 'misc')).toBe(true);
    expect(await shownValue(page, 'show-filtered-out-nodes', 'misc')).toBe('true');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {timeoutMs: 10_000});

    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    expect((await v.countCanvasPixels(page, VIEWER_TYPE)).total).toBeGreaterThan(hidden);

    expect(await v.togglePropertyGridCheckbox(page, 'show-filtered-out-nodes', 'misc')).toBe(false);
    await v.resetFilters(page);
  });

  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});

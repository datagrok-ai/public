/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = 'Word-cloud';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const VIEWER_TYPE = 'Word cloud';
const datasetPath = 'System:DemoFiles/demog.csv';

const shownValue = (page: Page, prop: string) => v.propertyGridValue(page, prop);
const category = (page: Page, cat: string, probe: string) =>
  v.ensurePropertyCategory(page, VIEWER_NAME, cat, probe);

async function renderedWords(page: Page): Promise<{name: string; value: number}[]> {
  return page.evaluate(() => {
    const wc = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Word cloud') as any;
    const data = wc.chart?.getOption()?.series?.[0]?.data ?? [];
    return data.map((d: any) => ({name: String(d.name), value: Number(d.value)}));
  });
}

const selectionCount = (page: Page) =>
  page.evaluate(() => (window as any).grok.shell.t.selection.trueCount as number);

const clearSelection = async (page: Page) => {
  await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
  await expect.poll(() => selectionCount(page), {timeout: 5000}).toBe(0);
};

async function wordPositions(page: Page): Promise<{x: number; y: number}[]> {
  return page.evaluate((sel) => {
    const root = document.querySelector(sel) as HTMLElement;
    const cv = root.querySelector('canvas') as HTMLCanvasElement;
    const r = cv.getBoundingClientRect();
    const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
    const sx = cv.width / r.width, sy = cv.height / r.height;
    const buckets = new Map<string, {n: number; x: number; y: number}>();
    for (let y = 0; y < cv.height; y += 2) {
      for (let x = 0; x < cv.width; x += 2) {
        const i = (y * cv.width + x) * 4;
        if (img[i + 3] === 0) continue;
        if (img[i] > 240 && img[i + 1] > 240 && img[i + 2] > 240) continue;
        const key = `${Math.floor(x / 50)}:${Math.floor(y / 30)}`;
        const b = buckets.get(key) ?? {n: 0, x: 0, y: 0};
        b.n++; b.x += x; b.y += y;
        buckets.set(key, b);
      }
    }
    return [...buckets.values()].filter((b) => b.n > 15).sort((a, b) => b.n - a.n).slice(0, 6)
      .map((b) => ({x: r.x + (b.x / b.n) / sx, y: r.y + (b.y / b.n) / sy}));
  }, VIEWER);
}

async function columnCounts(page: Page, column: string): Promise<Record<string, number>> {
  return page.evaluate((col) => {
    const c = (window as any).grok.shell.t.col(col);
    const counts: Record<string, number> = {};
    for (let i = 0; i < c.length; i++) {
      const value = String(c.get(i));
      counts[value] = (counts[value] ?? 0) + 1;
    }
    return counts;
  }, column);
}

test('Word cloud', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Add Word cloud from the Viewers toolbox', async () => {
    await page.locator('[name="icon-Word-cloud"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});

    await expect.poll(async () => (await v.countCanvasPixels(page, VIEWER_TYPE)).total,
      {timeout: 30_000}).toBeGreaterThan(500);

    await expect.poll(async () => (await renderedWords(page)).length,
      {timeout: 3000, intervals: [100, 150]}).toBeGreaterThan(0);
    const words = await renderedWords(page);
    const expected = await columnCounts(page, await shownValue(page, 'column') || 'SEX');
    expect(words.length).toBeGreaterThan(0);
    for (const w of words)
      expect(w.value).toBe(expected[w.name]);
  });

  await softStep('Column RACE draws one word per race with its row count', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);
    await category(page, 'data', 'column');
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    await v.pickColumnViaSelectorTrusted(page, {
      role: 'column', columnName: 'RACE', viewerType: VIEWER_TYPE,
      propName: 'columnColumnName', scopeSelector: '.property-grid',
    });
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500});

    const expected = await columnCounts(page, 'RACE');

    await expect.poll(async () => (await renderedWords(page)).length,
      {timeout: 3000, intervals: [100, 150]}).toBe(Object.keys(expected).length);
    const words = await renderedWords(page);
    expect(words.length).toBe(Object.keys(expected).length);
    for (const w of words)
      expect(w.value).toBe(expected[w.name]);
  });

  await softStep('A high-cardinality column shows an error instead of a cloud', async () => {
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'column', columnName: 'USUBJID', viewerType: VIEWER_TYPE,
      propName: 'columnColumnName', scopeSelector: '.property-grid',
    });

    const error = page.locator(`${VIEWER} .d4-viewer-error`).first();
    await expect(error).toBeVisible({timeout: 15_000});
    expect(await error.innerText()).toContain('500 or fewer unique categories');

    await v.pickColumnViaSelectorTrusted(page, {
      role: 'column', columnName: 'RACE', viewerType: VIEWER_TYPE,
      propName: 'columnColumnName', scopeSelector: '.property-grid',
    });
    await expect(page.locator(`${VIEWER} .d4-viewer-error`)).toHaveCount(0, {timeout: 15_000});
    await expect.poll(async () => (await renderedWords(page)).length,
      {timeout: 15_000}).toBeGreaterThan(0);
  });

  await softStep('Shape rearranges the words', async () => {
    await category(page, 'misc', 'shape');
    const shapes: Record<string, number> = {};
    for (const shape of ['star', 'diamond', 'circle']) {
      await v.waitForCanvasQuiet(page, VIEWER_TYPE);
      await v.snapshotCanvasColors(page, VIEWER_TYPE);
      await v.selectPropertyGridChoice(page, 'shape', shape);
      expect(await shownValue(page, 'shape')).toBe(shape);
      shapes[shape] = await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});
    }
    for (const shape of Object.keys(shapes))
      expect(shapes[shape]).toBeGreaterThan(100);
  });

  await softStep('Equal min and max text size flattens the word sizes', async () => {
    await category(page, 'misc', 'min-text-size');
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.setPropertyGridValue(page, 'min-text-size', '20');
    await v.setPropertyGridValue(page, 'max-text-size', '20');
    expect(await shownValue(page, 'min-text-size')).toBe('20');
    expect(await shownValue(page, 'max-text-size')).toBe('20');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});

    await v.setPropertyGridValue(page, 'min-text-size', '12');
    await v.setPropertyGridValue(page, 'max-text-size', '48');
  });

  await softStep('Rotation range changes the word angles', async () => {
    await category(page, 'misc', 'min-rotation-degree');
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.setPropertyGridValue(page, 'min-rotation-degree', '0');
    await v.setPropertyGridValue(page, 'max-rotation-degree', '0');
    expect(await shownValue(page, 'min-rotation-degree')).toBe('0');
    expect(await shownValue(page, 'max-rotation-degree')).toBe('0');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});
  });

  await softStep('Bold and font family restyle the words', async () => {
    await category(page, 'misc', 'bold');

    const boldWas = await shownValue(page, 'bold');
    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    expect(String(await v.togglePropertyGridCheckbox(page, 'bold'))).not.toBe(boldWas);
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});

    await v.waitForCanvasQuiet(page, VIEWER_TYPE);
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.selectPropertyGridChoice(page, 'font-family', 'monospace');
    expect(await shownValue(page, 'font-family')).toBe('monospace');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});

    expect(String(await v.togglePropertyGridCheckbox(page, 'bold'))).toBe(boldWas);
    await v.selectPropertyGridChoice(page, 'font-family', 'sans-serif');
  });

  await softStep('Hovering a word shows its row count', async () => {
    const words = await wordPositions(page);
    expect(words.length).toBeGreaterThan(0);

    let tooltip = {visible: false, text: ''};
    for (const w of words) {
      await page.mouse.move(w.x - 20, w.y - 20);
      await page.mouse.move(w.x, w.y, {steps: 6});

      tooltip = await v.pollValue(() => page.evaluate(() => {
        const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
        return {visible: !!t && getComputedStyle(t).display !== 'none', text: t?.innerText ?? ''};
      }), (tip) => tip.visible, 1200);
      if (tooltip.visible) break;
    }
    expect(tooltip.visible).toBe(true);
    expect(tooltip.text).toMatch(/\d+\s+rows/);

    const shownCount = parseInt(tooltip.text.match(/(\d+)\s+rows/)![1], 10);
    expect((await renderedWords(page)).map((w) => w.value)).toContain(shownCount);
  });

  await softStep('Clicking a word selects exactly that word\'s rows', async () => {
    await clearSelection(page);

    const words = await renderedWords(page);
    const counts = words.map((w) => w.value);
    let selected = 0;
    for (const position of await wordPositions(page)) {
      await page.mouse.click(position.x, position.y);

      selected = await v.pollValue(() => selectionCount(page), (n) => n > 0);
      if (selected > 0) break;
    }
    expect(selected).toBeGreaterThan(0);
    expect(counts).toContain(selected);
  });

  await softStep('Filtering the table leaves the cloud rendering', async () => {
    const before = await renderedWords(page);
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['M']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeGreaterThan(0);
    expect(filteredCount).toBeLessThan(total);

    await expect.poll(async () => (await renderedWords(page)).length,
      {timeout: 3000, intervals: [100, 150]}).toBe(before.length);
    const after = await renderedWords(page);
    expect(after.map((w) => w.name).sort()).toEqual(before.map((w) => w.name).sort());
    expect(await page.locator(`${VIEWER} .d4-viewer-error`).count()).toBe(0);
    expect((await v.countCanvasPixels(page, VIEWER_TYPE)).total).toBeGreaterThan(500);

    await v.resetFilters(page);
  });

  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.closeAllAndWait(page);

  v.finishSpec();
});

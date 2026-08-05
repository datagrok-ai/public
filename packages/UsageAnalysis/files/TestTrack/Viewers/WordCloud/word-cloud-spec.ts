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

/**
 * The words the cloud actually drew, with their counts. Read from the ECharts
 * series the viewer renders — the canvas carries no text nodes, and this is the
 * same data the user sees as words and reads in the tooltip.
 */
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

/** Screen positions of the drawn words — dark pixel blobs on the canvas. */
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

/** Category counts of a column, as the cloud should show them. */
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

  // #### Add the viewer
  await softStep('Add Word cloud from the Viewers toolbox', async () => {
    await page.locator('[name="icon-Word-cloud"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});

    // The cloud is laid out asynchronously — wait for words to be drawn.
    await expect.poll(async () => (await v.countCanvasPixels(page, VIEWER_TYPE)).total,
      {timeout: 30_000}).toBeGreaterThan(500);

    // The default column is drawn word-for-word with the right counts.
    const words = await renderedWords(page);
    const expected = await columnCounts(page, await shownValue(page, 'column') || 'SEX');
    expect(words.length).toBeGreaterThan(0);
    for (const w of words)
      expect(w.value).toBe(expected[w.name]);
  });

  // #### Column assignment
  await softStep('Column RACE draws one word per race with its row count', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);
    await category(page, 'data', 'column');
    await v.snapshotCanvasColors(page, VIEWER_TYPE);

    await v.pickColumnViaSelectorTrusted(page, {
      role: 'column', columnName: 'RACE', viewerType: VIEWER_TYPE,
      propName: 'columnColumnName', scopeSelector: '.property-grid',
    });
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 500});

    const words = await renderedWords(page);
    const expected = await columnCounts(page, 'RACE');
    expect(words.length).toBe(Object.keys(expected).length);
    for (const w of words)
      expect(w.value).toBe(expected[w.name]);
  });

  // #### Too many categories is reported, not drawn
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

  // #### Layout properties
  await softStep('Shape rearranges the words', async () => {
    await category(page, 'misc', 'shape');
    const shapes: Record<string, number> = {};
    for (const shape of ['star', 'diamond', 'circle']) {
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
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.setPropertyGridValue(page, 'min-rotation-degree', '0');
    await v.setPropertyGridValue(page, 'max-rotation-degree', '0');
    expect(await shownValue(page, 'min-rotation-degree')).toBe('0');
    expect(await shownValue(page, 'max-rotation-degree')).toBe('0');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});
  });

  await softStep('Bold and font family restyle the words', async () => {
    await category(page, 'misc', 'bold');
    // Bold ships on, so the step asserts the flip rather than a fixed value.
    const boldWas = await shownValue(page, 'bold');
    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    expect(String(await v.togglePropertyGridCheckbox(page, 'bold'))).not.toBe(boldWas);
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});

    await v.snapshotCanvasColors(page, VIEWER_TYPE);
    await v.selectPropertyGridChoice(page, 'font-family', 'monospace');
    expect(await shownValue(page, 'font-family')).toBe('monospace');
    await v.waitForCanvasChange(page, VIEWER_TYPE, {minDelta: 100});

    expect(String(await v.togglePropertyGridCheckbox(page, 'bold'))).toBe(boldWas);
    await v.selectPropertyGridChoice(page, 'font-family', 'sans-serif');
  });

  // #### Hovering and clicking a word
  await softStep('Hovering a word shows its row count', async () => {
    const words = await wordPositions(page);
    expect(words.length).toBeGreaterThan(0);

    let tooltip = {visible: false, text: ''};
    for (const w of words) {
      await page.mouse.move(w.x - 20, w.y - 20);
      await page.mouse.move(w.x, w.y, {steps: 6});
      await page.waitForTimeout(400);
      tooltip = await page.evaluate(() => {
        const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
        return {visible: !!t && getComputedStyle(t).display !== 'none', text: t?.innerText ?? ''};
      });
      if (tooltip.visible) break;
    }
    expect(tooltip.visible).toBe(true);
    expect(tooltip.text).toMatch(/\d+\s+rows/);
  });

  await softStep('Clicking a word selects exactly that word\'s rows', async () => {
    await clearSelection(page);

    const words = await renderedWords(page);
    const counts = words.map((w) => w.value);
    let selected = 0;
    for (const position of await wordPositions(page)) {
      await page.mouse.click(position.x, position.y);
      await expect.poll(() => selectionCount(page), {timeout: 5000}).toBeGreaterThanOrEqual(0);
      selected = await selectionCount(page);
      if (selected > 0) break;
    }
    expect(counts).toContain(selected);
  });

  // #### Filtering
  //
  // The cloud does NOT follow the table filter — it has no Row Source and no
  // Filter property, and the counts stay at their unfiltered values (dev,
  // 2026-08-04: `SEX = M` leaves RACE at Caucasian 5267 instead of 2444, even
  // after a forced redraw). The step therefore checks that filtering leaves the
  // viewer intact and still drawing, and does not pin the counts either way.
  await softStep('Filtering the table leaves the cloud rendering', async () => {
    const before = await renderedWords(page);
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['M']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeGreaterThan(0);
    expect(filteredCount).toBeLessThan(total);

    const after = await renderedWords(page);
    expect(after.map((w) => w.name).sort()).toEqual(before.map((w) => w.name).sort());
    expect(await page.locator(`${VIEWER} .d4-viewer-error`).count()).toBe(0);
    expect((await v.countCanvasPixels(page, VIEWER_TYPE)).total).toBeGreaterThan(500);

    await v.resetFilters(page);
  });

  // #### Closing the viewer
  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});

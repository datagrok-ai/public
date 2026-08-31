/* ---
realizes: [filters.cp.compose-with-viewer-filtering, filters.int.and-combination]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {expectHeaderCounter, expectHeaderCounterNow, trueCount} from '../../helpers/filter-panel';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const FULL = 5850;
const PANEL_CATEGORY = 'Asian';
const SPLIT_COLUMN = 'SEX';
const RACE_CARD_CANVAS = '.d4-filter [name="canvas"]';

async function raceSelectedCategories(page: Page): Promise<string[] | null> {
  return page.evaluate(() => {
    const states = grok.shell.tv.getFiltersGroup().getStates('RACE', 'categorical') as any[];
    if (!states || states.length === 0) return null;
    const s = states[0];
    return Array.isArray(s.selected) ? s.selected.map((c: any) => String(c)) : null;
  });
}

async function seedPanelCriterion(page: Page, category: string): Promise<number> {
  await page.evaluate(() => { grok.shell.tv.dataFrame.selection.setAll(false); });
  await v.resetFilters(page, {clearScatterFilter: true});
  const {filteredCount} = await v.applyCategoricalFilter(page, 'RACE', [category]);
  return filteredCount;
}

async function addViewer(page: Page, type: string): Promise<void> {
  await page.evaluate(async (t: string) => {
    grok.shell.tv.addViewer(t);
    await new Promise((r) => setTimeout(r, 1500));
  }, type);
}

async function viewerCanvasRect(page: Page, type: string):
    Promise<{x: number; y: number; w: number; h: number} | null> {
  return page.evaluate((t: string) => {
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === t);
    if (!vw) return null;
    const cv = vw.root.querySelector('canvas[name="canvas"]') || vw.root.querySelector('canvas');
    if (!cv) return null;
    const r = cv.getBoundingClientRect();
    return {x: r.x, y: r.y, w: r.width, h: r.height};
  }, type);
}

async function closeViewerByTitlebar(page: Page, viewerName: string, type: string): Promise<void> {
  await v.clickViewerTitlebarIcon(page, viewerName, 'Close');
  await expect.poll(async () => page.evaluate((t: string) =>
    Array.from(grok.shell.tv.viewers).filter((x: any) => x.type === t).length, type),
  {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(0);
}

async function driveViewerContextMenuLeaf(page: Page, type: string, path: string[]): Promise<boolean> {
  const rect = await viewerCanvasRect(page, type);
  if (!rect) return false;
  await page.mouse.click(rect.x + rect.w / 2, rect.y + rect.h / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 15_000});
  return page.evaluate(async (segments: string[]) => {
    const findLabel = (text: string) => {
      const popup = [...document.querySelectorAll('.d4-menu-popup')].pop();
      return [...(popup?.querySelectorAll('.d4-menu-item-label') ?? [])]
        .find((l) => (l.textContent ?? '').trim() === text) as HTMLElement | undefined;
    };
    const waitForLabel = async (text: string) => {
      const deadline = Date.now() + 5000;
      let l = findLabel(text);
      while (!l && Date.now() < deadline) {
        await new Promise((res) => setTimeout(res, 100));
        l = findLabel(text);
      }
      return l;
    };
    for (let i = 0; i < segments.length; i++) {
      const label = await waitForLabel(segments[i]);
      if (!label) return false;
      const item = (label.closest('.d4-menu-item') ?? label.closest('.d4-menu-group')) as HTMLElement | null;
      if (!item) return false;
      if (i === segments.length - 1) { item.click(); break; }
      const r = item.getBoundingClientRect();
      const at = {clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, bubbles: true};
      item.dispatchEvent(new MouseEvent('mouseenter', at));
      item.dispatchEvent(new MouseEvent('mousemove', at));
      await new Promise((res) => setTimeout(res, 400));
    }
    await new Promise((res) => setTimeout(res, 500));
    return true;
  }, path);
}

async function dragSliderHandle(
  page: Page, rootSelector: string, handleName: string, dx: number, dy: number,
): Promise<boolean> {
  const box = await page.evaluate(({sel, hn}: {sel: string; hn: string}) => {
    const root = document.querySelector(sel);
    if (!root) return null;
    const h = root.querySelector(`[name="${hn}"]`);
    if (!h) return null;
    const r = h.getBoundingClientRect();
    if (r.width === 0 && r.height === 0) return null;
    return {x: r.x + r.width / 2, y: r.y + r.height / 2};
  }, {sel: rootSelector, hn: handleName});
  if (!box) return false;
  await page.mouse.move(box.x, box.y);
  await page.mouse.down();
  await page.mouse.move(box.x + dx / 3, box.y + dy / 3, {steps: 4});
  await page.mouse.move(box.x + dx, box.y + dy, {steps: 10});
  await page.mouse.up();
  await page.waitForTimeout(900);
  return true;
}

async function armOnClickFilter(page: Page, type: string): Promise<void> {
  expect(await driveViewerContextMenuLeaf(page, type, ['On Click', 'Filter']),
    `the ${type} context menu has no On Click > Filter leaf`).toBe(true);
  await expect.poll(async () => page.evaluate((t: string) =>
    grok.shell.tv.viewers.find((x: any) => x.type === t)?.props?.onClick, type),
  {timeout: 10_000, intervals: [300, 600, 1200]}).toBe('Filter');
}

async function clickViewerAt(page: Page, type: string, fx: number, fy: number): Promise<boolean> {
  const rect = await viewerCanvasRect(page, type);
  if (!rect) return false;
  await page.mouse.click(rect.x + rect.w * fx, rect.y + rect.h * fy);
  await page.waitForTimeout(900);
  return true;
}

async function resetSliderByDoubleClick(page: Page, rootSelector: string): Promise<boolean> {
  const at = await page.evaluate((sel: string) => {
    const g = document.querySelector(sel);
    if (!g) return null;
    const r = g.getBoundingClientRect();
    if (r.width === 0 && r.height === 0) return null;
    return {x: r.x + r.width / 2, y: r.y + r.height / 2};
  }, rootSelector);
  if (!at) return false;
  await page.mouse.dblclick(at.x, at.y);
  await page.waitForTimeout(900);
  return true;
}

async function pollUntil<T>(
  read: () => Promise<T>, accept: (v: T) => boolean, timeoutMs = 6000,
): Promise<T> {
  const deadline = Date.now() + timeoutMs;
  for (;;) {
    const value = await read();
    if (accept(value) || Date.now() >= deadline) return value;
    await new Promise((r) => setTimeout(r, 200));
  }
}

const HISTOGRAM_CANVAS = '[name="viewer-Histogram"] canvas[name="canvas"]';
const HISTOGRAM_SLIDER_INSET = 25;
const PC_SLIDER = '[name="viewer-PC-Plot"] .d4-pc-plot-filter-slider svg[type="range-slider"]';

async function findHistogramMaxHandleX(
  page: Page, rect: {x: number; y: number; w: number; h: number}, y: number,
): Promise<number | null> {
  const readCursor = () => page.evaluate((sel: string) =>
    (document.querySelector(sel) as HTMLElement | null)?.style.cursor ?? '', HISTOGRAM_CANVAS);
  await page.mouse.move(rect.x + rect.w / 2, rect.y + rect.h / 2);
  for (let pass = 0; pass < 2; pass++) {
    for (let x = rect.x + rect.w - 2; x > rect.x; x -= 3) {
      await page.mouse.move(x, y);
      if (await pollUntil(readCursor, (c) => c === 'ew-resize', 200) === 'ew-resize') return x;
    }
  }
  return null;
}

async function dragHistogramHandle(page: Page, fromX: number, y: number, toX: number): Promise<void> {
  await page.mouse.move(fromX, y);
  await page.mouse.down();
  const steps = 12;
  for (let i = 1; i <= steps; i++)
    await page.mouse.move(fromX + (toX - fromX) * i / steps, y);
  await page.mouse.up();
}

const TRELLIS_CELL = '[name="viewer-Trellis-plot"] .d4-trellis-plot-cell';
const TRELLIS_CELL_CLICK_FY = 0.85;

async function trellisCurrentCells(page: Page): Promise<number> {
  return page.evaluate(() => document.querySelectorAll('.d4-trellis-cell-current').length);
}

async function clickTrellisCellFiltering(page: Page): Promise<boolean> {
  const base = await trueCount(page);
  const cells = page.locator(TRELLIS_CELL);
  const count = await cells.count();
  expect(count, 'the Trellis plot exposes no .d4-trellis-plot-cell to click').toBeGreaterThan(0);
  let delivered = false;
  for (let i = 0; i < count; i++) {
    const box = await cells.nth(i).boundingBox();
    if (!box) continue;
    await page.mouse.click(box.x + box.width * 0.5, box.y + box.height * TRELLIS_CELL_CLICK_FY);
    if (await pollUntil(() => trellisCurrentCells(page), (n) => n > 0, 3000) === 0) continue;
    delivered = true;
    const now = await pollUntil(() => trueCount(page), (n) => n !== base, 6000);
    if (now > 0 && now < base) return true;
    await page.keyboard.press('Escape');
    const reverted = await pollUntil(() => trueCount(page), (n) => n === base, 6000);
    expect(reverted, 'Escape did not put the filtered row count back to the value the Trellis cell ' +
      'was clicked from, so the next cell would be measured against a baseline that is no longer the ' +
      'pre-click count').toBe(base);
  }
  expect(delivered,
    'no Trellis cell click ever made a cell current (.d4-trellis-cell-current) — the click is ' +
    'landing on the chart drawn inside the cell instead of on the grid').toBe(true);
  return false;
}

const VIEWER_CHANNELS: Array<{
  type: string; domName: string;
  arm: (p: Page) => Promise<void>;
  drive: (p: Page) => Promise<boolean>;
  revert: (p: Page) => Promise<boolean>;
}> = [
  {
    type: 'Histogram', domName: 'Histogram',
    arm: async (p) => {
      await p.evaluate(async () => {
        const h = grok.shell.tv.viewers.find((x: any) => x.type === 'Histogram');
        h.props.valueColumnName = 'AGE';
        h.props.filteringEnabled = true;
        await new Promise((r) => setTimeout(r, 900));
      });
      expect(await p.evaluate(() =>
        grok.shell.tv.viewers.find((x: any) => x.type === 'Histogram').props.filteringEnabled),
      'the histogram is not routing its range slider to the dataframe filter').toBe(true);
    },
    drive: async (p) => {
      const rect = await viewerCanvasRect(p, 'Histogram');
      if (!rect) return false;
      const y = rect.y + rect.h - HISTOGRAM_SLIDER_INSET;
      const handleX = await findHistogramMaxHandleX(p, rect, y);
      if (handleX === null) return false;
      await dragHistogramHandle(p, handleX, y, handleX - rect.w * 0.42);
      return true;
    },
    revert: async (p) => {
      const rect = await viewerCanvasRect(p, 'Histogram');
      if (!rect) return false;
      const y = rect.y + rect.h - HISTOGRAM_SLIDER_INSET;
      const handleX = await findHistogramMaxHandleX(p, rect, y);
      if (handleX === null) return false;
      await dragHistogramHandle(p, handleX, y, rect.x + rect.w - 1);
      return true;
    },
  },
  {
    type: 'PC Plot', domName: 'PC-Plot',
    arm: async (p) => {
      await p.evaluate(async () => {
        const pc = grok.shell.tv.viewers.find((x: any) => x.type === 'PC Plot');
        pc.props.showFilters = true;
        await new Promise((r) => setTimeout(r, 900));
      });
      expect(await p.evaluate(() =>
        grok.shell.tv.viewers.find((x: any) => x.type === 'PC Plot').props.showFilters),
      'the PC Plot is not showing its per-axis filter sliders — there is no handle to drag').toBe(true);
    },
    drive: async (p) => {
      await p.locator('[name="viewer-PC-Plot"]').first().hover();
      await p.waitForTimeout(600);
      const height = await p.evaluate((sel: string) =>
        document.querySelector(sel)?.getBoundingClientRect().height ?? 0, PC_SLIDER);
      if (height === 0) return false;
      return dragSliderHandle(p, PC_SLIDER, 'max-handle', 0, height * 0.45);
    },
    revert: async (p) => resetSliderByDoubleClick(p, PC_SLIDER),
  },
  {
    type: 'Trellis plot', domName: 'Trellis-plot',
    arm: async (p) => armOnClickFilter(p, 'Trellis plot'),
    drive: clickTrellisCellFiltering,
    revert: async (p) => {
      await p.keyboard.press('Escape');
      return true;
    },
  },
  {
    type: 'Pie chart', domName: 'Pie-chart',
    arm: async (p) => armOnClickFilter(p, 'Pie chart'),
    drive: async (p) => clickViewerAt(p, 'Pie chart', 0.62, 0.38),
    revert: async (p) => clickViewerAt(p, 'Pie chart', 0.02, 0.02),
  },
];

async function closeFilterPanel(page: Page): Promise<void> {
  await page.locator('[name="viewer-Filters"]').first().hover();
  let clicked = false;
  for (const icon of ['icon-times', 'Close']) {
    try {
      await v.clickViewerTitlebarIcon(page, 'Filters', icon);
      clicked = true;
      break;
    } catch (_) { /* try the other title-bar close control */ }
  }
  expect(clicked, 'the Filters panel exposes no title-bar close control').toBe(true);
  await expect.poll(async () => page.locator('[name="viewer-Filters"]').count(),
    {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(0);
}

test('Filters — Composition of Filter Panel with Viewer-Driven Filtering', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  const total = await trueCount(page);
  expect(total).toBe(FULL);
  const truePanel = await seedPanelCriterion(page, PANEL_CATEGORY);
  expect(truePanel).toBeLessThan(FULL);
  expect(truePanel).toBeGreaterThan(0);
  await expectHeaderCounter(page, '1',
    'the seeded RACE criterion is the only thing filtering, so the header counter must read 1');

  let trueAfterZoom = -1;
  let rowCount = -1;
  let secondCategoryCount = -1;
  let selectionAfterAdd = -1;
  let inverted = -1;

  try {
    await softStep('Scenario 1 Step 4 Scatter zoom narrows below panel count, counter still 1', async () => {
      await addViewer(page, 'Scatter plot');
      const mode = await page.evaluate(() =>
        grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot').props.zoomAndFilter);
      expect(mode).toBe('filter by zoom');

      const rect = await viewerCanvasRect(page, 'Scatter plot');
      expect(rect, 'the Scatter Plot has no canvas to zoom on').not.toBeNull();
      const before = await trueCount(page);
      expect(before).toBe(truePanel);

      const x1 = rect!.x + rect!.w * 0.30, y1 = rect!.y + rect!.h * 0.30;
      const x2 = rect!.x + rect!.w * 0.70, y2 = rect!.y + rect!.h * 0.70;
      await page.mouse.move(x1, y1);
      await page.mouse.down();
      await page.mouse.move((x1 + x2) / 2, (y1 + y2) / 2, {steps: 6});
      await page.mouse.move(x2, y2, {steps: 8});
      await page.mouse.up();
      await expect.poll(async () => trueCount(page),
        {timeout: 15_000, intervals: [400, 800, 1500]}).toBeLessThan(truePanel);
      trueAfterZoom = await trueCount(page);

      expect(trueAfterZoom).toBeLessThan(truePanel);
      expect(trueAfterZoom).toBeGreaterThan(0);
      await expectHeaderCounterNow(page, '1',
        'the scatter zoom is not a panel card, so the header counter must still read 1');
      expect(await raceSelectedCategories(page)).toEqual([PANEL_CATEGORY]);
    });

    await softStep('Scenario 1 Step 5 Change panel criterion, zoom still biting (github-2642 guard)', async () => {
      const zoomedFullCaucasian = await page.evaluate(() => {
        const raceCol = grok.shell.tv.dataFrame.col('RACE');
        let n = 0;
        for (let i = 0; i < grok.shell.tv.dataFrame.rowCount; i++) if (raceCol.get(i) === 'Caucasian') n++;
        return n;
      });
      const {filteredCount: changed} = await v.applyCategoricalFilter(page, 'RACE', ['Caucasian'], 1000);
      expect(changed).not.toBe(trueAfterZoom);
      expect(changed).toBeGreaterThan(0);
      expect(changed).toBeLessThan(zoomedFullCaucasian);
      await expectHeaderCounterNow(page, '1',
        'changing the criterion on the one card must leave the header counter at 1');

      const {filteredCount: restored} = await v.applyCategoricalFilter(page, 'RACE', [PANEL_CATEGORY], 900);
      expect(restored).toBeLessThan(truePanel);
    });

    await softStep('Scenario 1 Step 8 Close Scatter Plot → trueCount returns to panel count exactly', async () => {
      const beforeClose = await trueCount(page);
      expect(beforeClose).toBeLessThan(truePanel);
      await closeViewerByTitlebar(page, 'Scatter-plot', 'Scatter plot');
      await expect.poll(async () => trueCount(page),
        {timeout: 10_000, intervals: [400, 800, 1500]}).toBe(truePanel);
      expect(await trueCount(page)).toBe(truePanel);
      await expectHeaderCounterNow(page, '1',
        'closing the scatter plot leaves the panel criterion alone filtering, so the counter must read 1');
      expect(await raceSelectedCategories(page)).toEqual([PANEL_CATEGORY]);
    });

    await softStep('Scenario 1 Step 7 Bar-chart click narrows below panel count, counter still 1', async () => {
      await addViewer(page, 'Bar chart');
      await page.evaluate(async (split: string) => {
        const bc = grok.shell.tv.viewers.find((x: any) => x.type === 'Bar chart');
        bc.props.splitColumnName = split;
        await new Promise((r) => setTimeout(r, 1000));
      }, SPLIT_COLUMN);
      await armOnClickFilter(page, 'Bar chart');

      const rect = await viewerCanvasRect(page, 'Bar chart');
      expect(rect, 'the Bar Chart has no canvas to click').not.toBeNull();
      const before = await trueCount(page);
      expect(before,
        'the bar chart is not starting from the panel-only base — a previous viewer is still filtering')
        .toBe(truePanel);

      const probeBars = async () => page.evaluate(({cat, split}: {cat: string; split: string}) => {
        const df = grok.shell.tv.dataFrame;
        const splitCol = df.col(split);
        const raceCol = df.col('RACE');
        const counts: Record<string, number> = {};
        for (const c of splitCol.categories) counts[c] = 0;
        for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[splitCol.get(i)]++;
        const survivors = Object.keys(counts).filter((c) => counts[c] > 0);
        let expected = -1;
        if (survivors.length === 1) {
          expected = 0;
          for (let i = 0; i < df.rowCount; i++)
            if (raceCol.get(i) === cat && splitCol.get(i) === survivors[0]) expected++;
        }
        return {survivors, expected, filtered: df.filter.trueCount};
      }, {cat: PANEL_CATEGORY, split: SPLIT_COLUMN});

      const positions = [
        {x: rect!.x + rect!.w * 0.5, y: rect!.y + rect!.h * 0.4},
        {x: rect!.x + rect!.w * 0.5, y: rect!.y + rect!.h * 0.6},
        {x: rect!.x + rect!.w * 0.35, y: rect!.y + rect!.h * 0.5},
        {x: rect!.x + rect!.w * 0.65, y: rect!.y + rect!.h * 0.5},
        {x: rect!.x + rect!.w * 0.5, y: rect!.y + rect!.h * 0.75},
      ];
      let landed: {survivors: string[]; expected: number; filtered: number} | null = null;
      for (const pos of positions) {
        await page.mouse.click(pos.x, pos.y);
        await page.waitForTimeout(900);
        const probe = await probeBars();
        if (probe.survivors.length === 1 && probe.filtered > 0 && probe.filtered < truePanel) {
          landed = probe;
          break;
        }
      }
      expect(landed,
        `no bar-chart click left exactly one ${SPLIT_COLUMN} category standing below the panel count ` +
        `(last read: ${JSON.stringify(await probeBars())})`).not.toBeNull();
      expect(landed!.filtered,
        `the surviving rows are not exactly the rows that are both ${PANEL_CATEGORY} and ` +
        `${landed!.survivors[0]}, so the bar-chart click and the panel criterion did not intersect`)
        .toBe(landed!.expected);
      expect(landed!.filtered).toBeLessThan(truePanel);
      expect(landed!.filtered).toBeGreaterThan(0);
      await expectHeaderCounterNow(page, '1',
        'the bar-chart click is not a panel card, so the header counter must still read 1');
    });

    await softStep('Scenario 1 Step 9 Close Bar Chart → trueCount returns to panel count exactly', async () => {
      const beforeClose = await trueCount(page);
      expect(beforeClose,
        'the bar-chart narrowing is not in place, so closing it cannot be shown to release anything')
        .toBeLessThan(truePanel);
      await closeViewerByTitlebar(page, 'Bar-chart', 'Bar chart');
      await expect.poll(async () => trueCount(page),
        {timeout: 10_000, intervals: [400, 800, 1500]}).toBe(truePanel);
      expect(await trueCount(page)).toBe(truePanel);
      await expectHeaderCounterNow(page, '1',
        'closing the bar chart leaves the panel criterion alone filtering, so the counter must read 1');
      expect(await raceSelectedCategories(page)).toEqual([PANEL_CATEGORY]);
    });

    for (const viewer of VIEWER_CHANNELS) {
      await softStep(`Scenario 4 ${viewer.type} composes with the panel criterion, then reverts`, async () => {
        const before = await trueCount(page);
        expect(before,
          `${viewer.type} is not starting from the panel-only base — something else is still filtering`)
          .toBe(truePanel);

        await addViewer(page, viewer.type);
        expect(await page.evaluate((t: string) =>
          Array.from(grok.shell.tv.viewers).some((x: any) => x.type === t), viewer.type),
        `the ${viewer.type} viewer was not added`).toBe(true);

        try {
          await viewer.arm(page);
          expect(await viewer.drive(page),
            `the ${viewer.type} filtering gesture could not be ISSUED — its target (canvas rect, ` +
            'slider handle or cell box) was not on screen').toBe(true);

          await expect.poll(async () => trueCount(page),
            {timeout: 15_000, intervals: [400, 800, 1500]}).toBeLessThan(truePanel);
          const composed = await trueCount(page);
          expect(composed).toBeLessThan(truePanel);
          expect(composed, `the ${viewer.type} gesture filtered every row away`).toBeGreaterThan(0);
          await expectHeaderCounterNow(page, '1',
            `the ${viewer.type} gesture is not a panel card, so the header counter must still read 1`);
          expect(await raceSelectedCategories(page)).toEqual([PANEL_CATEGORY]);

          expect(await viewer.revert(page),
            `the ${viewer.type} revert gesture could not be ISSUED — its target was not on screen`)
            .toBe(true);
          await expect.poll(async () => trueCount(page),
            {timeout: 15_000, intervals: [400, 800, 1500]}).toBe(truePanel);
          expect(await trueCount(page)).toBe(truePanel);
          expect(await raceSelectedCategories(page)).toEqual([PANEL_CATEGORY]);
        } finally {
          await closeViewerByTitlebar(page, viewer.domName, viewer.type);
        }
        await expect.poll(async () => trueCount(page),
          {timeout: 10_000, intervals: [400, 800, 1500]}).toBe(truePanel);
        await expectHeaderCounterNow(page, '1',
          `closing the ${viewer.type} leaves the panel criterion alone filtering, so the counter reads 1`);
      });
    }

    await softStep('Scenario 2 Step 5 Layout round-trip with combined filtering → no error, count restored', async () => {
      let savedLayout: any = null;
      try {
        await seedPanelCriterion(page, PANEL_CATEGORY);
        await addViewer(page, 'Scatter plot');
        const rect = await viewerCanvasRect(page, 'Scatter plot');
        expect(rect).not.toBeNull();
        await page.mouse.move(rect!.x + rect!.w * 0.30, rect!.y + rect!.h * 0.30);
        await page.mouse.down();
        await page.mouse.move(rect!.x + rect!.w * 0.70, rect!.y + rect!.h * 0.70, {steps: 8});
        await page.mouse.up();
        await expect.poll(async () => trueCount(page),
          {timeout: 15_000, intervals: [400, 800, 1500]}).toBeLessThan(truePanel);
        const zoomed = await trueCount(page);
        const beforeSave = await page.evaluate(() => {
          const tv = grok.shell.tv;
          const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
          return {
            filtered: tv.dataFrame.filter.trueCount,
            zoomAndFilter: sp?.props?.zoomAndFilter ?? null,
            xMin: sp?.props?.xMin ?? null, xMax: sp?.props?.xMax ?? null,
            yMin: sp?.props?.yMin ?? null, yMax: sp?.props?.yMax ?? null,
            spFilter: sp ? String(sp.props.filter ?? '') : null,
          };
        });

        const pageErrors: string[] = [];
        const onPageError = (e: Error) => pageErrors.push(e.message);
        page.on('pageerror', onPageError);
        await page.evaluate(() => {
          const w = window as any;
          w.__errBalloons = [];
          const seen = new WeakSet<Element>();
          const record = (el: Element) => {
            if (seen.has(el)) return;
            seen.add(el);
            w.__errBalloons.push((el.textContent ?? '').trim());
          };
          const scan = (n: Node) => {
            if (!(n instanceof Element)) return;
            if (n.matches('.d4-balloon.error')) record(n);
            for (const el of Array.from(n.querySelectorAll('.d4-balloon.error'))) record(el);
          };
          w.__errBalloonObs = new MutationObserver((records: MutationRecord[]) => {
            for (const r of records) {
              if (r.type === 'attributes') {
                const t = r.target as Element;
                if (t.matches('.d4-balloon.error')) record(t);
              } else
                for (const n of Array.from(r.addedNodes)) scan(n);
            }
          });
          w.__errBalloonObs.observe(document.body,
            {childList: true, subtree: true, attributes: true, attributeFilter: ['class']});
        });

        const beforeIds = await page.evaluate(async () => {
          const me = String(grok.shell.user.id);
          const ls = (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)) ?? [];
          return ls.filter((l: any) => !l.author || !l.author.id || String(l.author.id) === me).map((l: any) => String(l.id));
        });

        expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Save to Gallery'])).toBe(true);
        let fresh: string[] = [];
        await expect.poll(async () => {
          fresh = (await page.evaluate(async (prev: string[]) => {
            const me = String(grok.shell.user.id);
            const ls = (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)) ?? [];
            return ls.filter((l: any) => !l.author || !l.author.id || String(l.author.id) === me)
              .map((l: any) => String(l.id)).filter((id: string) => !prev.includes(id));
          }, beforeIds));
          return fresh.length;
        }, {timeout: 25_000, intervals: [500, 1000, 2000, 3000]}).toBeGreaterThanOrEqual(1);
        expect(fresh.length, `expected exactly 1 new layout, got ${fresh.length}`).toBe(1);
        savedLayout = fresh[0];

        await closeFilterPanel(page);

        await page.evaluate(async (layoutId: string) => {
          const saved = await grok.dapi.layouts.find(layoutId);
          const applied = new Promise<void>((resolve) => {
            const sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); resolve(); });
            setTimeout(resolve, 8000);
          });
          grok.shell.tv.loadLayout(saved);
          await applied;
        }, savedLayout);

        await expect.poll(async () => page.locator('[name="viewer-Filters"] .d4-filter').count(),
          {timeout: 20_000, intervals: [500, 1000, 2000, 3000]}).toBeGreaterThanOrEqual(1);
        await expect.poll(async () => raceSelectedCategories(page),
          {timeout: 15_000, intervals: [500, 1000, 2000]}).toEqual([PANEL_CATEGORY]);
        const afterLayout = await page.evaluate(() => {
          const tv = grok.shell.tv;
          const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
          return {
            viewers: Array.from(tv.viewers).map((x: any) => x.type),
            filtered: tv.dataFrame.filter.trueCount,
            zoomAndFilter: sp?.props?.zoomAndFilter ?? null,
            xMin: sp?.props?.xMin ?? null, xMax: sp?.props?.xMax ?? null,
            yMin: sp?.props?.yMin ?? null, yMax: sp?.props?.yMax ?? null,
            spFilter: sp ? String(sp.props.filter ?? '') : null,
          };
        });
        expect(afterLayout.viewers,
          'the re-applied layout did not bring the Scatter Plot back, so the row count settling at ' +
          'the panel-only value below says nothing about a dropped zoom — it would be the count of a ' +
          `view with no scatter plot in it at all; viewers: [${afterLayout.viewers.join(', ')}]`)
          .toContain('Scatter plot');
        expect(afterLayout.zoomAndFilter,
          'the restored Scatter Plot is no longer routing its zoom to the dataframe filter, so the ' +
          'panel-only row count below would be explained by the viewer being disarmed rather than by ' +
          'the zoom itself being dropped from the layout').toBe('filter by zoom');
        expect(zoomed, 'the pre-save zoom was not narrower than the panel criterion alone, so the ' +
          'round-trip bound below could not tell a restored zoom from a lost one').toBeLessThan(truePanel);
        await expect.poll(async () => trueCount(page), {
          timeout: 30_000,
          intervals: [500, 1000, 2000, 3000],
          message: 'the re-applied layout did not settle at the panel-only row count: ' +
            `before=${JSON.stringify(beforeSave)} after=${JSON.stringify(afterLayout)}`,
        }).toBe(truePanel);
        const restored = await trueCount(page);
        expect(restored,
          'the layout round-trip no longer settles at the panel-only row count. Below it means the ' +
          'scatter-plot zoom is now restored with the layout (scope_reductions SR-06 is obsolete — ' +
          'tighten this to toBeLessThan(truePanel)); above it means the panel criterion itself was ' +
          `lost. before=${JSON.stringify(beforeSave)} after=${JSON.stringify(afterLayout)}`)
          .toBe(truePanel);

        const errBalloons = await page.evaluate(async () => {
          const w = window as any;
          let last = (w.__errBalloons ?? []).length;
          let since = Date.now();
          while (Date.now() - since < 2500) {
            await new Promise((r) => setTimeout(r, 200));
            const n = (w.__errBalloons ?? []).length;
            if (n !== last) { last = n; since = Date.now(); }
          }
          w.__errBalloonObs?.disconnect();
          return (w.__errBalloons ?? []) as string[];
        });
        page.off('pageerror', onPageError);
        expect(errBalloons,
          `error balloons appeared during the layout round-trip — GROK-18281: ${errBalloons.join(' | ')}`)
          .toEqual([]);
        expect(pageErrors, `layout re-apply raised page errors — GROK-18281: ${pageErrors.join('; ')}`).toEqual([]);
      } finally {
        await page.evaluate(() => { try { (window as any).__errBalloonObs?.disconnect(); } catch (_) {} });
        if (savedLayout) {
          await page.evaluate(async (layoutId: string) => {
            try { const s = await grok.dapi.layouts.find(layoutId); await grok.dapi.layouts.delete(s); } catch (_) {}
          }, savedLayout);
        }
        await page.evaluate(() => { try { grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot')?.close(); } catch (_) {} });
      }
    });

    await softStep('Scenario 3 Step 2 Select the visible rows → selection equals the panel-filtered set', async () => {
      const seeded = await seedPanelCriterion(page, PANEL_CATEGORY);
      expect(seeded, 'the panel-only base was not re-established before the selection block')
        .toBe(truePanel);
      const stray = await page.evaluate(() => Array.from(grok.shell.tv.viewers)
        .map((x: any) => x.type).filter((t: string) => t !== 'Grid' && t !== 'Filters'));
      expect(stray, `viewers left over from an earlier block are still filtering: ${stray.join(', ')}`).toEqual([]);

      const sel = await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        df.selection.copyFrom(df.filter);
        await new Promise((r) => setTimeout(r, 400));
        return {afterVisible: df.selection.trueCount, rowCount: df.rowCount};
      });
      rowCount = sel.rowCount;
      expect(rowCount, 'demog reports no rows, so every count below would be vacuous').toBeGreaterThan(0);
      expect(sel.afterVisible, 'selecting the visible rows did not select the panel-filtered set')
        .toBe(truePanel);
    });

    await softStep('Scenario 3 Step 3 Add the second category to the selection', async () => {
      const sel = await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        const raceCol = df.col('RACE');
        const kept: boolean[] = [];
        for (let i = 0; i < df.rowCount; i++) kept.push(df.selection.get(i));
        let caucasian = 0;
        for (let i = 0; i < df.rowCount; i++) if (raceCol.get(i) === 'Caucasian') caucasian++;
        df.selection.init((i: number) => kept[i] || raceCol.get(i) === 'Caucasian');
        await new Promise((r) => setTimeout(r, 400));
        return {afterAdd: df.selection.trueCount, caucasian, filter: df.filter.trueCount};
      });
      secondCategoryCount = sel.caucasian;
      selectionAfterAdd = sel.afterAdd;
      expect(secondCategoryCount, 'the second category holds no rows, so growing the selection by it ' +
        'could be satisfied by a command that did nothing').toBeGreaterThan(0);
      expect(sel.afterAdd, 'adding the hidden category did not grow the selection by that category')
        .toBe(truePanel + secondCategoryCount);
      expect(sel.filter, 'adding rows to the selection moved the filter as well').toBe(truePanel);
    });

    await softStep('Scenario 3 Step 4 Select > Invert leaves exactly the complement', async () => {
      expect(await v.driveTopMenuLeaf(page, ['Select', 'Invert']),
        'the Select > Invert menu leaf was not actuated').toBe(true);
      await expect.poll(async () => page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount),
        {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(rowCount - selectionAfterAdd);
      inverted = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
      expect(inverted, 'the inversion did not leave exactly the complement').toBe(rowCount - selectionAfterAdd);
      expect(inverted).toBeGreaterThan(0);
      expect(inverted).toBeLessThan(rowCount);
      const filterBefore = await trueCount(page);
      expect(filterBefore).toBe(truePanel);
      expect(filterBefore, 'the filter already equals the inverted selection, so the Step 5 equality ' +
        'would be satisfied by a command that does nothing').not.toBe(inverted);
      expect(await raceSelectedCategories(page)).toEqual([PANEL_CATEGORY]);
    });

    await softStep('Scenario 3 Step 5 Selection to filter → count matches inverted selection, criterion kept, card repaints', async () => {
      const panel = await page.evaluate((cs: string) => {
        const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
        return {
          captions: cards.map((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() ?? ''),
          hasCanvas: !!document.querySelector('[name="viewer-Filters"]')?.querySelector(cs),
        };
      }, RACE_CARD_CANVAS);
      expect(panel.captions[0],
        `the first filter card is not RACE (cards: ${panel.captions.join(', ')})`).toBe('RACE');
      expect(panel.hasCanvas, `the RACE card has no ${RACE_CARD_CANVAS} to snapshot`).toBe(true);

      const look = await page.evaluate(() => {
        const tv = grok.shell.tv;
        const f = tv.viewers.find((x: any) => x.type === 'Filters');
        f.props.showCurrentRow = false;
        f.props.showMouseOverRow = false;
        f.props.showMouseOverGroupRow = false;
        tv.dataFrame.currentRowIdx = -1;
        tv.dataFrame.mouseOverRowIdx = -1;
        return {current: f.props.showCurrentRow, mouseOver: f.props.showMouseOverRow};
      });
      expect(look, 'the card row-marker rendering could not be switched off — the guard below would ' +
        'not tell a command repaint apart from a pointer repaint').toEqual({current: false, mouseOver: false});
      const park = await page.evaluate(() => {
        const g = document.querySelector('.d4-table-view [name="viewer-Grid"]') ?? document.body;
        const r = g.getBoundingClientRect();
        return {x: Math.round(r.x + r.width / 2), y: Math.round(r.y + r.height / 2)};
      });
      await page.mouse.move(park.x, park.y);

      await v.waitForCanvasQuiet(page, 'Filters', {canvasSelector: RACE_CARD_CANVAS, stableReads: 3});
      expect(await v.snapshotCanvasColors(page, 'Filters', RACE_CARD_CANVAS),
        'the RACE card canvas could not be read — the repaint guard below would be vacuous').toBe(true);
      await page.waitForTimeout(1500);
      const idle = await v.diffCanvasColors(page, 'Filters', RACE_CARD_CANVAS);
      expect(idle.deltaPx,
        'the RACE card canvas kept changing while nothing acted on it (or could not be read: -1) — ' +
        'a post-command delta would not be attributable to the command')
        .toBe(0);

      expect(await v.driveTopMenuLeaf(page, ['Select', 'Selection to Filter']),
        'the Select > Selection to Filter menu leaf was not actuated').toBe(true);
      await page.mouse.move(park.x, park.y);
      await page.waitForTimeout(1500);

      const outcome = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const raceCol = df.col('RACE');
        const surviving: Record<string, number> = {};
        for (const c of raceCol.categories) surviving[c] = 0;
        for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) surviving[raceCol.get(i)]++;
        const filteredOut = Object.keys(surviving).filter((c) => surviving[c] === 0);
        return {filterCount: df.filter.trueCount, surviving, filteredOut};
      });

      expect(outcome.filterCount,
        `surviving = ${JSON.stringify(outcome.surviving)}`).toBe(inverted);
      expect(outcome.filteredOut,
        `surviving = ${JSON.stringify(outcome.surviving)}`).toContain(PANEL_CATEGORY);
      const raceSelected = await raceSelectedCategories(page);
      expect(raceSelected, 'the RACE card exposes no criterion after selection-to-filter').not.toBeNull();
      expect(raceSelected!,
        'the command erased the RACE card criterion instead of only copying the selection onto the filter')
        .toEqual([PANEL_CATEGORY]);

      await v.waitForCanvasQuiet(page, 'Filters', {canvasSelector: RACE_CARD_CANVAS, stableReads: 3});
      const {deltaPx} = await v.diffCanvasColors(page, 'Filters', RACE_CARD_CANVAS);
      expect(deltaPx,
        'the RACE card canvas is identical to its pre-command frame (or could not be read: -1) — the ' +
        'card did not repaint after the filter was replaced from outside the panel (GROK-16713)')
        .toBeGreaterThan(0);
    });
  } finally {
    await page.evaluate(() => {
      try {
        grok.shell.tv?.dataFrame?.resetFilter();
        grok.shell.tv?.dataFrame?.selection?.setAll(false);
      } catch (_) {}
    });
    await v.cleanupShell(page);
  }

  v.finishSpec();
});

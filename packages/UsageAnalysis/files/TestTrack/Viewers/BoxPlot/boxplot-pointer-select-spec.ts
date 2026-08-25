/* ---
realizes: [boxplot.cp.pointer-select-highlight]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {armBalloonRecorderProved, expectNoBalloonSinceArmed} from '../../helpers/balloons';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const LABEL_Y_FRACS = [0.89, 0.885, 0.88, 0.895, 0.875, 0.9];
let labelYFrac = 0;

async function canvasRect(page: Page): Promise<{x: number; y: number; w: number; h: number}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: c.x, y: c.y, w: c.width, h: c.height};
  });
}

async function bpProp(page: Page, prop: string): Promise<any> {
  return page.evaluate((p) => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    return bp?.props?.[p];
  }, prop);
}

async function selectionCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.t.selection.trueCount);
}

async function clearSelection(page: Page): Promise<void> {
  await page.evaluate(() => grok.shell.t.selection.setAll(false));

  await page.waitForTimeout(300);
}

async function raceCatCount(page: Page, cat: string): Promise<number> {
  return page.evaluate((c) => {
    const col = grok.shell.t.col('RACE');
    let n = 0;
    for (let i = 0; i < grok.shell.t.rowCount; i++) if (col.get(i) === c) n++;
    return n;
  }, cat);
}

async function raceCategoriesInAxisOrder(page: Page): Promise<string[]> {
  return page.evaluate(() => grok.shell.t.col('RACE').categories.slice());
}

/** Category slots do NOT span the full canvas — the value axis takes a left margin, so
 *  `(i + 0.5) / cats.length` lands outside the plot area. Measured once by clicking each
 *  candidate centre in the label band and keeping the offsets that select a whole category. */
let catCenters: number[] = [];

async function measureGeometry(page: Page): Promise<void> {
  if (labelYFrac > 0 && catCenters.length > 0) return;
  const cats = await raceCategoriesInAxisOrder(page);
  const counts = await Promise.all(cats.map((c) => raceCatCount(page, c)));
  const r = await canvasRect(page);

  for (const frac of LABEL_Y_FRACS) {
    // the plot area starts after the value-axis margin; scan plausible left edges
    for (const left of [0.15, 0.12, 0.18, 0.10, 0.20, 0.0]) {
      const width = 1 - left - 0.05;
      const centers = cats.map((_, i) => left + width * ((i + 0.5) / cats.length));
      await clearSelection(page);
      await page.mouse.click(r.x + r.w * centers[0], r.y + r.h * frac);
      const selected = await v.pollValue(() => selectionCount(page), (n) => n === counts[0], 600, 50);
      if (selected === counts[0]) {
        labelYFrac = frac;
        catCenters = centers;
        break;
      }
    }
    if (catCenters.length > 0) break;
  }
  await clearSelection(page);
  if (labelYFrac === 0 || catCenters.length === 0)
    throw new Error('category-label band not found: no candidate Y fraction selected the full first category');
}

async function findLabelYFrac(page: Page): Promise<number> {
  await measureGeometry(page);
  return labelYFrac;
}

/** X fraction of category [i]'s centre, measured rather than assumed. */
async function catCenter(page: Page, i: number): Promise<number> {
  await measureGeometry(page);
  return catCenters[i];
}

async function shiftDragBand(
  page: Page, r: {x: number; y: number; w: number; h: number},
  x0: number, y0: number, x1: number, y1: number,
): Promise<void> {
  await page.keyboard.down('Shift');
  await page.mouse.move(r.x + r.w * x0, r.y + r.h * y0);
  await page.mouse.down();
  await page.mouse.move(r.x + r.w * x1, r.y + r.h * y1, {steps: 14});
  await page.mouse.up();
  await page.keyboard.up('Shift');
  await v.waitForViewerRendered(page, 'Box plot', 800);
}

test('Box plot pointer selection and highlight', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const tv = grok.shell.tv;
    const bp = tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
    bp.props.category1ColumnName = 'RACE';
    bp.props.markerSize = 10;
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await v.waitForViewerRendered(page, 'Box plot', 1500);
  await v.waitForCanvasQuiet(page, 'Box plot');

  await softStep('Scenario 1 / Step 3: marker click sets currentRowIdx and fires d4-boxplot-point-click', async () => {
    await page.evaluate(() => {
      const w = window as any;
      w.__pointClickRowId = null;
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      w.__pointClickSub = bp.onEvent('d4-boxplot-point-click').subscribe((a: any) => {
        w.__pointClickRowId = a?.args?.rowId ?? a?.rowId ?? null;
      });
    });
    const r = await canvasRect(page);
    let currentRow = -1;
    for (const [fx, fy] of [[0.62, 0.55], [0.62, 0.45], [0.62, 0.65], [0.55, 0.5], [0.7, 0.5]]) {
      await page.mouse.click(r.x + r.w * fx, r.y + r.h * fy);

      currentRow = await v.pollValue(
        () => page.evaluate(() => grok.shell.t.currentRowIdx), (n) => n >= 0, 500, 50);
      if (currentRow >= 0) break;
    }
    const eventRowId = await page.evaluate(() => (window as any).__pointClickRowId);
    await page.evaluate(() => (window as any).__pointClickSub?.unsubscribe());
    expect(currentRow).toBeGreaterThanOrEqual(0);
    expect(eventRowId).toBe(currentRow);
  });

  await softStep('Scenario 2 / Step 4: shift-drag selects with visible hue (baseline coloring)', async () => {
    await clearSelection(page);

    const hueNoSel = await v.countSelectionHuePixels(page, 'Box plot');
    const r = await canvasRect(page);

    await shiftDragBand(page, r, 0.55, 0.30, 0.72, 0.80);
    const sel = await selectionCount(page);
    const huePx = await v.pollValue(() => v.countSelectionHuePixels(page, 'Box plot'),
      (px) => px > hueNoSel + 2000, 3000, 150);
    expect(sel).toBeGreaterThan(0);
    expect(huePx).toBeGreaterThan(hueNoSel + 2000);
  });

  await softStep('Scenario 2 / Step 5: selection hue survives non-default categorical coloring (github-3066)', async () => {

    await clearSelection(page);

    await v.setViewerProps(page, 'Box plot', [{set: {markerColorColumnName: 'SEX'}, wait: 800}]);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const hueNoSel = await v.countSelectionHuePixels(page, 'Box plot');
    const r = await canvasRect(page);
    await shiftDragBand(page, r, 0.55, 0.30, 0.72, 0.80);
    const sel = await selectionCount(page);
    const huePx = await v.pollValue(() => v.countSelectionHuePixels(page, 'Box plot'),
      (px) => px > hueNoSel + 2000, 3000, 150);
    expect(sel).toBeGreaterThan(0);
    expect(huePx).toBeGreaterThan(hueNoSel + 2000);
    await v.setViewerProps(page, 'Box plot', [{set: {markerColorColumnName: ''}, wait: 500}]);
  });

  await softStep('Scenario 3 / Step 6: ctrl-click a DISJOINT category adds its rows (additive union)', async () => {
    const labelY = await findLabelYFrac(page);
    await clearSelection(page);
    const cats = await raceCategoriesInAxisOrder(page);
    expect(cats.length).toBeGreaterThanOrEqual(2);

    const r = await canvasRect(page);
    const firstCenter = await catCenter(page, 0);
    let baseline = 0;
    for (const [dx, y0, y1] of [[0.05, 0.42, 0.60], [0.08, 0.30, 0.80], [0.10, 0.20, 0.85]]) {
      await shiftDragBand(page, r, firstCenter - dx, y0, firstCenter + dx, y1);
      baseline = await selectionCount(page);
      if (baseline > 0) break;
    }
    expect(baseline).toBeGreaterThan(0);
    const secondCat = cats[1];
    const secondCount = await raceCatCount(page, secondCat);
    expect(secondCount).toBeGreaterThan(0);
    const cx = r.x + r.w * (await catCenter(page, 1));
    await page.keyboard.down('Control');
    await page.mouse.click(cx, r.y + r.h * labelY);
    await page.keyboard.up('Control');
    await v.waitForViewerRendered(page, 'Box plot', 700);
    expect(await selectionCount(page)).toBe(baseline + secondCount);
  });

  await softStep('Scenario 3 / Step 7: plain category-label click replaces selection with the full category', async () => {
    const labelY = await findLabelYFrac(page);
    const cats = await raceCategoriesInAxisOrder(page);
    const firstCount = await raceCatCount(page, cats[0]);
    const r = await canvasRect(page);
    const cx = r.x + r.w * (await catCenter(page, 0));
    await page.mouse.click(cx, r.y + r.h * labelY);
    await v.waitForViewerRendered(page, 'Box plot', 700);
    expect(await selectionCount(page)).toBe(firstCount);
  });

  await softStep('Scenario 3 / Step 8 and Step 9: single click on empty space clears selection, no reset-view', async () => {
    await page.evaluate(() => {
      const w = window as any;
      w.__resetEvents = [];
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      w.__resetSub = bp.onEvent('d4-boxplot-reset-view').subscribe(() => { w.__resetEvents.push(Date.now()); });
    });
    const r = await canvasRect(page);

    await page.mouse.click(r.x + r.w * 0.5, r.y + r.h * 0.05);
    await v.waitForViewerRendered(page, 'Box plot', 700);
    expect(await selectionCount(page)).toBe(0);

    const resetEventsAfterSingle = await page.evaluate(() => (window as any).__resetEvents.length);
    expect(resetEventsAfterSingle).toBe(0);
  });

  await softStep('Scenario 3 / Step 10: double-click on empty space resets the view (fires d4-boxplot-reset-view)', async () => {
    const r = await canvasRect(page);
    await page.mouse.dblclick(r.x + r.w * 0.5, r.y + r.h * 0.05);

    const resetEventsAfterDouble = await v.pollValue(
      () => page.evaluate(() => (window as any).__resetEvents.length), (n) => n > 0, 800, 50);
    expect(resetEventsAfterDouble).toBeGreaterThan(0);
    await page.evaluate(() => (window as any).__resetSub?.unsubscribe());
  });

  await softStep('Scenario 4 / Step 12: no selection leak into filtered-out (left-most) category (github-2764)', async () => {
    await clearSelection(page);
    const cats = await raceCategoriesInAxisOrder(page);
    const leftMost = cats[0];
    await page.evaluate((cat) => {
      const df = grok.shell.t;
      const col = df.col('RACE');
      df.filter.init((i: number) => col.get(i) !== cat);
    }, leftMost);
    await v.waitForViewerRendered(page, 'Box plot', 800);
    await v.waitForCanvasQuiet(page, 'Box plot');

    const r = await canvasRect(page);
    await page.keyboard.down('Shift');
    await page.mouse.move(r.x + r.w * 0.05, r.y + r.h * 0.30);
    await page.mouse.down();
    await page.mouse.move(r.x + r.w * 0.95, r.y + r.h * 0.80, {steps: 22});
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await v.waitForViewerRendered(page, 'Box plot', 800);
    const leaked = await page.evaluate((cat) => {
      const df = grok.shell.t;
      const col = df.col('RACE');
      let leak = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (col.get(i) === cat && df.selection.get(i)) leak++;
      return leak;
    }, leftMost);
    expect(leaked).toBe(0);
    await page.evaluate(() => grok.shell.t.filter.setAll(true));
    await v.waitForViewerRendered(page, 'Box plot', 500);
  });

  await softStep('Scenario 5: baseline (categorical coloring, category select adds selection hue)', async () => {
    const labelY = await findLabelYFrac(page);
    await clearSelection(page);
    await v.setViewerProps(page, 'Box plot', [{set: {markerColorColumnName: 'RACE'}, wait: 800}]);
    await v.waitForCanvasQuiet(page, 'Box plot');

    const hueBaseline = await v.countSelectionHuePixels(page, 'Box plot');
    await page.evaluate(() => (window as any).__hueBaseline = null);
    const cats = await raceCategoriesInAxisOrder(page);
    const r = await canvasRect(page);
    const cauIdx = cats.indexOf('Caucasian') >= 0 ? cats.indexOf('Caucasian') : cats.length - 1;
    const cauCount = await raceCatCount(page, cats[cauIdx]);
    await page.mouse.click(r.x + r.w * (await catCenter(page, cauIdx)), r.y + r.h * labelY);
    await v.waitForViewerRendered(page, 'Box plot', 800);
    await v.waitForCanvasQuiet(page, 'Box plot');
    expect(await selectionCount(page)).toBe(cauCount);
    const hueSelected = await v.countSelectionHuePixels(page, 'Box plot');
    expect(hueSelected).toBeGreaterThan(hueBaseline + 2000);
    await page.evaluate((b) => (window as any).__hueBaseline = b, hueBaseline);
  });

  await softStep('Scenario 5 / Step 14: showSelectedRows=false suppresses the selection hue (GROK-19950)', async () => {
    const hueBaseline = await page.evaluate(() => (window as any).__hueBaseline);
    const hueWithSel = await v.countSelectionHuePixels(page, 'Box plot');
    await v.setViewerProps(page, 'Box plot', [{set: {showSelectedRows: false}, wait: 800}]);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const hueSuppressed = await v.countSelectionHuePixels(page, 'Box plot');

    expect(hueSuppressed).toBeLessThan(hueWithSel);
    expect(hueSuppressed).toBeLessThanOrEqual(hueBaseline + Math.round(hueBaseline * 0.1) + 50);
    expect(hueSuppressed).toBeGreaterThan(0);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('RACE');
    console.log('Step 14 hue baseline/withSel/suppressed:', hueBaseline, hueWithSel, hueSuppressed);
    await v.setViewerProps(page, 'Box plot', [{set: {showSelectedRows: true}, wait: 500}]);
  });

  await softStep('Scenario 5 / Step 15: rowSource=Selected suppresses hue via showSelectedRows dependsOn gate (GROK-19950)', async () => {
    await v.setViewerProps(page, 'Box plot', [{set: {rowSource: 'Selected'}, wait: 900}]);
    await v.waitForCanvasQuiet(page, 'Box plot');

    expect(await v.countSelectionHuePixels(page, 'Box plot')).toBe(0);
    await v.setViewerProps(page, 'Box plot', [{set: {rowSource: 'All', markerColorColumnName: ''}, wait: 600}]);
  });

  await softStep('Scenario 6 / Step 17: hovering a marker shows a .d4-tooltip', async () => {
    await clearSelection(page);
    const r = await canvasRect(page);
    let tipPresent = false;
    for (const [fx, fy] of [[0.62, 0.55], [0.62, 0.45], [0.6, 0.6], [0.55, 0.5]]) {
      await page.mouse.move(r.x + r.w * fx, r.y + r.h * fy);

      tipPresent = await v.pollValue(
        () => page.evaluate(() => !!document.querySelector('.d4-tooltip')), (p) => p, 700, 50);
      if (tipPresent) break;
    }
    expect(tipPresent).toBe(true);
    await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
    await v.waitForViewerRendered(page, 'Box plot', 400);
  });

  await softStep('Scenario 6 / Step 18: showMouseOverPoint=false draws no hover overlay (no-error floor)', async () => {
    await v.setViewerProps(page, 'Box plot', [{set: {showMouseOverPoint: false}, wait: 600}]);
    const r = await canvasRect(page);
    await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
    await v.waitForCanvasQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await page.mouse.move(r.x + r.w * 0.62, r.y + r.h * 0.55);

    await page.waitForTimeout(700);
    const {deltaPx} = await v.diffCanvasColors(page, 'Box plot');

    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeLessThan(300);
    console.log('Step 18 showMouseOverPoint=false hover deltaPx:', deltaPx);
    await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
    await v.setViewerProps(page, 'Box plot', [{set: {showMouseOverPoint: true}, wait: 400}]);
  });

  await softStep('Scenario 7 / Step 20: cross-viewer bar hover does not repaint the box plot (github-3065)', async () => {
    await clearSelection(page);
    await page.evaluate(() => {
      const tv = grok.shell.tv;
      const bar = tv.addViewer('Bar chart');
      bar.props.splitColumnName = 'RACE';
      const bp = tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.markerColorColumnName = 'RACE';
    });
    await page.locator('[name="viewer-Bar-chart"]').waitFor({timeout: 10000});
    await v.waitForViewerRendered(page, 'Box plot', 1200);
    await v.waitForCanvasQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    const barRect = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Bar-chart"]')!;
      const c = root.querySelector('canvas')!.getBoundingClientRect();
      return {x: c.x, y: c.y, w: c.width, h: c.height};
    });
    await page.mouse.move(barRect.x + barRect.w * 0.3, barRect.y + barRect.h * 0.6);

    await page.waitForTimeout(900);
    const {deltaPx} = await v.diffCanvasColors(page, 'Box plot');
    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeLessThan(300);
    console.log('Step 20 cross-viewer hover box plot deltaPx:', deltaPx);
    await page.evaluate(() => {
      const tv = grok.shell.tv;
      const bar = tv.viewers.find((x: any) => x.type === 'Bar chart');
      if (bar) bar.close();
      const bp = tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.markerColorColumnName = '';
    });
    await v.waitForViewerRendered(page, 'Box plot', 600);
  });

  await softStep('Scenario 8 / Step 22: box plot survives categorical-coloring row deletion without error (GROK-20502)', async () => {
    const labelY = await findLabelYFrac(page);
    await clearSelection(page);
    await v.setViewerProps(page, 'Box plot', [{set: {markerColorColumnName: 'RACE'}, wait: 800}]);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const cats = await raceCategoriesInAxisOrder(page);
    const r = await canvasRect(page);
    const cauIdx = cats.indexOf('Caucasian') >= 0 ? cats.indexOf('Caucasian') : cats.length - 1;
    await page.mouse.click(r.x + r.w * (await catCenter(page, cauIdx)), r.y + r.h * labelY);
    await v.waitForViewerRendered(page, 'Box plot', 700);
    const selForDelete = await selectionCount(page);
    expect(selForDelete).toBeGreaterThan(0);
    // The probe is raised INSIDE the window this arm opens, so an empty reading below means nothing
    // was raised rather than nothing was watching. grok.shell.warnings does not exist (js-api/src/shell.ts:176 exposes only
    // warning()), so a before/after count on it is not an alternative — it compares undefined with
    // undefined. The probe's own balloon carries a marker the absence read filters out.
    await armBalloonRecorderProved(page, 'boxplot row-deletion');

    const rowsBefore = await page.evaluate(() => grok.shell.t.rowCount);

    await page.evaluate(() => {
      const df = grok.shell.t;
      const set = new Set(Array.from(df.selection.getSelectedIndexes()));
      df.rows.removeWhereIdx((i: number) => set.has(i));
    });
    await v.waitForViewerRendered(page, 'Box plot', 1000);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const rowsAfter = await page.evaluate(() => grok.shell.t.rowCount);
    expect(rowsAfter).toBe(rowsBefore - selForDelete);

    const px = await v.countCanvasPixels(page, 'Box plot', {canvasSelector: 'canvas[name="canvas"]'});
    expect(px.total).toBeGreaterThan(0);
    await expectNoBalloonSinceArmed(page, 'deleting the selected rows', /warning/);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('RACE');
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

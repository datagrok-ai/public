/* ---
realizes: [boxplot.cp.pointer-select-highlight]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// Category-axis label band (xAxisBox, resolved by _catTableHitTest,
// box_plot_core.dart:1975): a narrow strip just below the plot area; the ink
// rows further down are the click-inert statistics table. Its exact edge
// shifts with layout, so candidate Y fractions are probed empirically.
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

// Drives the category-label expectations from the dataframe, not the UI.
async function raceCatCount(page: Page, cat: string): Promise<number> {
  return page.evaluate((c) => {
    const col = grok.shell.t.col('RACE');
    let n = 0;
    for (let i = 0; i < grok.shell.t.rowCount; i++) if (col.get(i) === c) n++;
    return n;
  }, cat);
}

// The viewer exposes no public axis-order accessor; the dataframe category
// list is the axis (left→right) order for a single categorical column.
async function raceCategoriesInAxisOrder(page: Page): Promise<string[]> {
  return page.evaluate(() => grok.shell.t.col('RACE').categories.slice());
}

// Locate the category-label click band empirically (a candidate Y is confirmed
// when the click selects exactly the first category's row count). Cached
// across scenarios; leaves the selection cleared.
async function findLabelYFrac(page: Page): Promise<number> {
  if (labelYFrac > 0) return labelYFrac;
  const cats = await raceCategoriesInAxisOrder(page);
  const cat0Count = await raceCatCount(page, cats[0]);
  const r = await canvasRect(page);
  for (const frac of LABEL_Y_FRACS) {
    await clearSelection(page);
    await page.mouse.click(r.x + r.w * (0.5 / cats.length), r.y + r.h * frac);
    await page.waitForTimeout(500);
    if (await selectionCount(page) === cat0Count) { labelYFrac = frac; break; }
  }
  await clearSelection(page);
  if (labelYFrac === 0)
    throw new Error('category-label band not found: no candidate Y fraction selected the full first category');
  return labelYFrac;
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
  await page.waitForTimeout(800);
}

test('Box plot pointer selection and highlight', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // Setting Category1 auto-sets markerColorColumnName to the same column, so
  // the baseline coloring of this spec is the RACE categorical palette.
  await page.evaluate(() => {
    const tv = grok.shell.tv;
    const bp = tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
    bp.props.category1ColumnName = 'RACE';
    bp.props.markerSize = 10;
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await page.waitForTimeout(1500);
  await v.waitForCanvasQuiet(page, 'Box plot');

  // ==================================================================
  // Scenario 1: Marker click sets currentRowIdx and fires point-click
  // ==================================================================
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
      await page.waitForTimeout(500);
      currentRow = await page.evaluate(() => grok.shell.t.currentRowIdx);
      if (currentRow >= 0) break;
    }
    const eventRowId = await page.evaluate(() => (window as any).__pointClickRowId);
    await page.evaluate(() => (window as any).__pointClickSub?.unsubscribe());
    expect(currentRow).toBeGreaterThanOrEqual(0);
    expect(eventRowId).toBe(currentRow);
  });

  // ==================================================================
  // Scenario 2: Shift-drag creates a selection with visible hue
  // ==================================================================
  await softStep('Scenario 2 / Step 4: shift-drag selects with visible hue (baseline coloring)', async () => {
    await clearSelection(page);
    // The categorical palette itself paints pixels inside the orange
    // selection-hue range — the honest signal is the DELTA above the
    // no-selection baseline, not an absolute count.
    const hueNoSel = await v.countSelectionHuePixels(page, 'Box plot');
    const r = await canvasRect(page);
    // Band over the Caucasian column (rightmost, densest cloud).
    await shiftDragBand(page, r, 0.55, 0.30, 0.72, 0.80);
    const sel = await selectionCount(page);
    const huePx = await v.countSelectionHuePixels(page, 'Box plot');
    expect(sel).toBeGreaterThan(0);
    expect(huePx).toBeGreaterThan(hueNoSel + 2000);
  });

  await softStep('Scenario 2 / Step 5: selection hue survives non-default categorical coloring (github-3066)', async () => {
    // Clear the previous selection BEFORE measuring the no-selection baseline —
    // a stale selection band in the baseline would swallow the post-drag delta.
    await clearSelection(page);
    // SEX = a non-default categorical coloring (the auto-default is the
    // Category1 column); its palette paints its own orange-range pixels, so
    // the signal is again the delta above its no-selection baseline.
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.markerColorColumnName = 'SEX';
    });
    await page.waitForTimeout(800);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const hueNoSel = await v.countSelectionHuePixels(page, 'Box plot');
    const r = await canvasRect(page);
    await shiftDragBand(page, r, 0.55, 0.30, 0.72, 0.80);
    const sel = await selectionCount(page);
    const huePx = await v.countSelectionHuePixels(page, 'Box plot');
    expect(sel).toBeGreaterThan(0);
    expect(huePx).toBeGreaterThan(hueNoSel + 2000);
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.markerColorColumnName = '';
    });
    await page.waitForTimeout(500);
  });

  // ==================================================================
  // Scenario 3: Category label click — ctrl-click (additive), plain, clear
  // ==================================================================
  await softStep('Scenario 3 / Step 6: ctrl-click a DISJOINT category adds its rows (additive union)', async () => {
    const labelY = await findLabelYFrac(page);
    await clearSelection(page);
    const cats = await raceCategoriesInAxisOrder(page);
    expect(cats.length).toBeGreaterThanOrEqual(2);
    // Confine the shift-drag to the FIRST category's X band — the selected rows
    // are then disjoint from every other category by construction. Marker cloud
    // density varies, so widen the band until the drag catches rows.
    const r = await canvasRect(page);
    const firstCenter = 0.5 / cats.length;
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
    const cx = r.x + r.w * ((1 + 0.5) / cats.length);
    await page.keyboard.down('Control');
    await page.mouse.click(cx, r.y + r.h * labelY);
    await page.keyboard.up('Control');
    await page.waitForTimeout(700);
    expect(await selectionCount(page)).toBe(baseline + secondCount);
  });

  await softStep('Scenario 3 / Step 7: plain category-label click replaces selection with the full category', async () => {
    const labelY = await findLabelYFrac(page);
    const cats = await raceCategoriesInAxisOrder(page);
    const firstCount = await raceCatCount(page, cats[0]);
    const r = await canvasRect(page);
    const cx = r.x + r.w * ((0 + 0.5) / cats.length);
    await page.mouse.click(cx, r.y + r.h * labelY);
    await page.waitForTimeout(700);
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
    // Empty space = top strip above the boxes (no markers there).
    await page.mouse.click(r.x + r.w * 0.5, r.y + r.h * 0.05);
    await page.waitForTimeout(700);
    expect(await selectionCount(page)).toBe(0);
    // props.viewport is not populated by pointer interactions on this viewer —
    // the reset-view EVENT is the observable click-vs-dblclick signal.
    const resetEventsAfterSingle = await page.evaluate(() => (window as any).__resetEvents.length);
    expect(resetEventsAfterSingle).toBe(0);
  });

  await softStep('Scenario 3 / Step 10: double-click on empty space resets the view (fires d4-boxplot-reset-view)', async () => {
    const r = await canvasRect(page);
    await page.mouse.dblclick(r.x + r.w * 0.5, r.y + r.h * 0.05);
    await page.waitForTimeout(800);
    const resetEventsAfterDouble = await page.evaluate(() => (window as any).__resetEvents.length);
    expect(resetEventsAfterDouble).toBeGreaterThan(0);
    await page.evaluate(() => (window as any).__resetSub?.unsubscribe());
  });

  // ==================================================================
  // Scenario 4: No selection leak into filtered-out categories (github-2764)
  // ==================================================================
  await softStep('Scenario 4 / Step 12: no selection leak into filtered-out (left-most) category (github-2764)', async () => {
    await clearSelection(page);
    const cats = await raceCategoriesInAxisOrder(page);
    const leftMost = cats[0];
    await page.evaluate((cat) => {
      const df = grok.shell.t;
      const col = df.col('RACE');
      df.filter.init((i: number) => col.get(i) !== cat);
    }, leftMost);
    await page.waitForTimeout(800);
    await v.waitForCanvasQuiet(page, 'Box plot');
    // The drag covers the left edge where the filtered-out category previously rendered.
    const r = await canvasRect(page);
    await page.keyboard.down('Shift');
    await page.mouse.move(r.x + r.w * 0.05, r.y + r.h * 0.30);
    await page.mouse.down();
    await page.mouse.move(r.x + r.w * 0.95, r.y + r.h * 0.80, {steps: 22});
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await page.waitForTimeout(800);
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
    await page.waitForTimeout(500);
  });

  // ==================================================================
  // Scenario 5: showSelectedRows and rowSource gate (GROK-19950)
  // ==================================================================
  await softStep('Scenario 5: baseline (categorical coloring, category select adds selection hue)', async () => {
    const labelY = await findLabelYFrac(page);
    await clearSelection(page);
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.markerColorColumnName = 'RACE';
    });
    await page.waitForTimeout(800);
    await v.waitForCanvasQuiet(page, 'Box plot');
    // The categorical palette itself paints orange-band pixels — the signal is
    // the DELTA above this no-selection baseline, not an absolute > 0.
    const hueBaseline = await v.countSelectionHuePixels(page, 'Box plot');
    await page.evaluate(() => (window as any).__hueBaseline = null);
    const cats = await raceCategoriesInAxisOrder(page);
    const r = await canvasRect(page);
    const cauIdx = cats.indexOf('Caucasian') >= 0 ? cats.indexOf('Caucasian') : cats.length - 1;
    const cauCount = await raceCatCount(page, cats[cauIdx]);
    await page.mouse.click(r.x + r.w * ((cauIdx + 0.5) / cats.length), r.y + r.h * labelY);
    await page.waitForTimeout(800);
    await v.waitForCanvasQuiet(page, 'Box plot');
    expect(await selectionCount(page)).toBe(cauCount);
    const hueSelected = await v.countSelectionHuePixels(page, 'Box plot');
    expect(hueSelected).toBeGreaterThan(hueBaseline + 2000);
    await page.evaluate((b) => (window as any).__hueBaseline = b, hueBaseline);
  });

  await softStep('Scenario 5 / Step 14: showSelectedRows=false suppresses the selection hue (GROK-19950)', async () => {
    const hueBaseline = await page.evaluate(() => (window as any).__hueBaseline);
    const hueWithSel = await v.countSelectionHuePixels(page, 'Box plot');
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.showSelectedRows = false;
    });
    await page.waitForTimeout(800);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const hueSuppressed = await v.countSelectionHuePixels(page, 'Box plot');
    // The hue returns to the coloring-only baseline, NOT to 0 — the categorical
    // palette keeps painting orange-range pixels; that non-zero remainder is the
    // evidence that regular color coding is preserved (the coloring legend is
    // canvas-drawn, no legend DOM exists to inspect).
    expect(hueSuppressed).toBeLessThan(hueWithSel);
    expect(hueSuppressed).toBeLessThanOrEqual(hueBaseline + Math.round(hueBaseline * 0.1) + 50);
    expect(hueSuppressed).toBeGreaterThan(0);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('RACE');
    console.log('Step 14 hue baseline/withSel/suppressed:', hueBaseline, hueWithSel, hueSuppressed);
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.showSelectedRows = true;
    });
    await page.waitForTimeout(500);
  });

  await softStep('Scenario 5 / Step 15: rowSource=Selected suppresses hue via showSelectedRows dependsOn gate (GROK-19950)', async () => {
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.rowSource = 'Selected';
    });
    await page.waitForTimeout(900);
    await v.waitForCanvasQuiet(page, 'Box plot');
    // With rowSource=Selected the coloring-only remainder is not rendered, so
    // exactly zero is honest here — distinct from the baseline-return of the
    // showSelectedRows=false step.
    expect(await v.countSelectionHuePixels(page, 'Box plot')).toBe(0);
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.rowSource = 'All';
      bp.props.markerColorColumnName = '';
    });
    await page.waitForTimeout(600);
  });

  // ==================================================================
  // Scenario 6: Hover tooltip and showMouseOverPoint
  // ==================================================================
  await softStep('Scenario 6 / Step 17: hovering a marker shows a .d4-tooltip', async () => {
    await clearSelection(page);
    const r = await canvasRect(page);
    let tipPresent = false;
    for (const [fx, fy] of [[0.62, 0.55], [0.62, 0.45], [0.6, 0.6], [0.55, 0.5]]) {
      await page.mouse.move(r.x + r.w * fx, r.y + r.h * fy);
      await page.waitForTimeout(700);
      tipPresent = await page.evaluate(() => !!document.querySelector('.d4-tooltip'));
      if (tipPresent) break;
    }
    expect(tipPresent).toBe(true);
    await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
    await page.waitForTimeout(400);
  });

  await softStep('Scenario 6 / Step 18: showMouseOverPoint=false draws no hover overlay (no-error floor)', async () => {
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.showMouseOverPoint = false;
    });
    await page.waitForTimeout(600);
    const r = await canvasRect(page);
    await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
    await v.waitForCanvasQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await page.mouse.move(r.x + r.w * 0.62, r.y + r.h * 0.55);
    await page.waitForTimeout(700);
    const {deltaPx} = await v.diffCanvasColors(page, 'Box plot');
    // Fault guard BEFORE the ceiling (-1 passes any toBeLessThan).
    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeLessThan(300);
    console.log('Step 18 showMouseOverPoint=false hover deltaPx:', deltaPx);
    await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.showMouseOverPoint = true;
    });
    await page.waitForTimeout(400);
  });

  // ==================================================================
  // Scenario 7: Cross-viewer highlight must not recolor box plot (github-3065)
  // ==================================================================
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
    await page.waitForTimeout(1200);
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
    await page.waitForTimeout(600);
  });

  // ==================================================================
  // Scenario 8: Coloring survives row deletion (GROK-20502) — canvas legend
  // ==================================================================
  await softStep('Scenario 8 / Step 22: box plot survives categorical-coloring row deletion without error (GROK-20502)', async () => {
    const labelY = await findLabelYFrac(page);
    await clearSelection(page);
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      bp.props.markerColorColumnName = 'RACE';
    });
    await page.waitForTimeout(800);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const cats = await raceCategoriesInAxisOrder(page);
    const r = await canvasRect(page);
    const cauIdx = cats.indexOf('Caucasian') >= 0 ? cats.indexOf('Caucasian') : cats.length - 1;
    await page.mouse.click(r.x + r.w * ((cauIdx + 0.5) / cats.length), r.y + r.h * labelY);
    await page.waitForTimeout(700);
    const selForDelete = await selectionCount(page);
    expect(selForDelete).toBeGreaterThan(0);
    const warnBefore = await page.evaluate(() => grok.shell.warnings?.length ?? 0);
    const rowsBefore = await page.evaluate(() => grok.shell.t.rowCount);
    // The selected indices are non-contiguous — remove via removeWhereIdx over
    // the captured selection, not a contiguous removeAt block.
    await page.evaluate(() => {
      const df = grok.shell.t;
      const set = new Set(Array.from(df.selection.getSelectedIndexes()));
      df.rows.removeWhereIdx((i: number) => set.has(i));
    });
    await page.waitForTimeout(1000);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const rowsAfter = await page.evaluate(() => grok.shell.t.rowCount);
    expect(rowsAfter).toBe(rowsBefore - selForDelete);
    // The coloring legend is CANVAS-drawn (no legend DOM to inspect) — the
    // observable is that the coloring box plot still renders: canvas keeps
    // painting, no new shell warnings, coloring column still applied.
    const warnAfter = await page.evaluate(() => grok.shell.warnings?.length ?? 0);
    const px = await v.countCanvasPixels(page, 'Box plot', {canvasSelector: 'canvas[name="canvas"]'});
    expect(px.total).toBeGreaterThan(0);
    expect(warnAfter).toBe(warnBefore);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('RACE');
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

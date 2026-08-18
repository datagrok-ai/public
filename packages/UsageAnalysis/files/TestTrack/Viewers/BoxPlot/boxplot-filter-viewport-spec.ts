/* ---
realizes: [boxplot.cp.filter-viewport]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

const OVERLAY = 'canvas[name="overlay"]';

async function viewportRect(page: Page): Promise<{y: number; height: number}> {
  return page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const vp = bp.viewport;
    return {y: vp.y, height: vp.height};
  });
}

async function bpProp(page: Page, prop: string): Promise<any> {
  return page.evaluate((p) => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    return bp?.props?.[p];
  }, prop);
}

async function setBpProp(page: Page, prop: string, value: any, settleMs = 1200): Promise<void> {
  await v.setViewerProps(page, 'Box plot', [{set: {[prop]: value}, wait: settleMs}]);
}

async function narrowAverageMass(page: Page, lo: number, hi: number): Promise<number> {
  return page.evaluate(({lo, hi}) => {
    const df = grok.shell.t;
    const col = df.col('Average Mass');
    df.filter.init((i: number) => { const x = col.get(i); return x != null && x >= lo && x <= hi; });
    return df.filter.trueCount;
  }, {lo, hi});
}

async function resetFilter(page: Page): Promise<void> {
  await page.evaluate(() => grok.shell.t.filter.setAll(true));
  await v.pollValue(
    () => page.evaluate(() => ({n: grok.shell.t.filter.trueCount, total: grok.shell.t.rowCount})),
    (c) => c.n === c.total, 800, 100);
}

async function averageMassFilterCanvas(page: Page): Promise<{x: number; y: number; w: number; h: number}> {
  await page.evaluate(() => {
    const fg = grok.shell.tv.getFiltersGroup();
    fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'Average Mass'});
  });
  await v.pollValue(() => page.evaluate(() => {
    const filtersRoot = document.querySelector('[name="viewer-Filters"]');
    return Array.from(filtersRoot?.querySelectorAll('.d4-filter') ?? []).some((c) =>
      c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'Average Mass');
  }), (present) => present, 1500, 100);
  const readRect = () => page.evaluate(() => {
    const filtersRoot = document.querySelector('[name="viewer-Filters"]')!;
    const cards = Array.from(filtersRoot.querySelectorAll('.d4-filter'));
    const card = cards.find((c) =>
      c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'Average Mass');
    if (!card) return null;
    card.scrollIntoView();
    const cv = card.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
    if (!cv) return null;
    const r = cv.getBoundingClientRect();
    return {x: r.x, y: r.y, w: r.width, h: r.height};
  });
  let rect: {x: number; y: number; w: number; h: number} | null = null;
  let prev = '';
  for (let i = 0; i < 20; i++) {
    rect = await readRect();
    const sig = JSON.stringify(rect);
    if (rect && rect.w > 0 && rect.h > 0 && sig === prev) break;
    prev = sig;
    await page.waitForTimeout(400);
  }
  if (!rect) throw new Error('Average Mass histogram filter canvas not found');
  return rect;
}

async function openFilterPanelWithBudget(page: Page, budgetMs = 90_000): Promise<void> {
  const deadline = Date.now() + budgetMs;
  for (let attempt = 1; ; attempt++) {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    const hasCards = await page.waitForFunction(
      () => document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length > 0,
      null, {timeout: 10_000, polling: 300},
    ).then(() => true).catch(() => false);
    const count = await page.evaluate(() =>
      document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length);
    console.log(`Filter panel open attempt ${attempt}: ${count} filter cards`);
    if (hasCards) {

      await page.locator('[name="viewer-Filters"] .d4-filter').first().hover().catch(() => {});
      return;
    }
    if (Date.now() >= deadline)
      throw new Error(`Filter panel produced no filter cards within ${budgetMs}ms (${attempt} attempts)`);

    await page.evaluate(() => {
      const f = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Filters') as any;
      if (f) f.close();
    });
    await v.pollValue(
      () => page.evaluate(() => Array.from(grok.shell.tv.viewers).some((x: any) => x.type === 'Filters')),
      (present) => !present, 1000, 100);
  }
}

async function dragFilterHandle(
  page: Page, rect: {x: number; y: number; w: number; h: number},
  side: 'min' | 'max', targetFrac: number,
): Promise<boolean> {
  await page.mouse.move(rect.x + rect.w * 0.5, rect.y + rect.h * 0.5);

  await page.waitForTimeout(400);
  const cursorNow = () => page.evaluate(() => {
    const filtersRoot = document.querySelector('[name="viewer-Filters"]')!;
    const cards = Array.from(filtersRoot.querySelectorAll('.d4-filter'));
    const card = cards.find((c) =>
      c.querySelector('.d4-filter-column-name')?.textContent?.trim() === 'Average Mass');
    const cv = card?.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
    return cv?.style.cursor ?? '';
  });
  for (let fy = 0.5; fy <= 0.98; fy += 0.06) {
    for (let px = 2; px <= 24; px += 3) {
      const x = side === 'min' ? rect.x + px : rect.x + rect.w - px;
      const y = rect.y + rect.h * fy;
      await page.mouse.move(x, y);
      await page.waitForTimeout(60); 
      if (await cursorNow() === 'ew-resize') {
        const before = await filterTrueCount(page);
        await page.mouse.down();
        await page.mouse.move(rect.x + rect.w * targetFrac, y, {steps: 12});
        await page.mouse.up();

        await v.pollValue(() => filterTrueCount(page), (n) => n !== before, 1200, 100);
        return true;
      }
    }
  }
  return false;
}

async function filterTrueCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.t.filter.trueCount);
}

async function categoryBandColors(page: Page): Promise<Record<string, number>> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const cv = root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
    const ctx = cv.getContext('2d')!;
    const y0 = Math.floor(cv.height * 0.80);
    const img = ctx.getImageData(0, y0, cv.width, cv.height - y0).data;
    const colors: Record<string, number> = {};
    for (let i = 0; i < img.length; i += 4) {
      const key = String((img[i] << 16) | (img[i + 1] << 8) | img[i + 2]);
      colors[key] = (colors[key] ?? 0) + 1;
    }
    return colors;
  });
}

function bandColorsDelta(a: Record<string, number>, b: Record<string, number>): number {
  let d = 0;
  for (const k of new Set([...Object.keys(a), ...Object.keys(b)]))
    d += Math.abs((a[k] ?? 0) - (b[k] ?? 0));
  return d;
}

async function hasEmptyValuedCategory(page: Page): Promise<boolean> {
  return page.evaluate(() => {
    const df = grok.shell.t;
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const series = df.col('Series');
    const val = df.col(bp.props.valueColumnName);
    const per: Record<string, {n: number; nulls: number}> = {};
    for (let i = 0; i < df.rowCount; i++) {
      const c = String(series.get(i) ?? '');
      per[c] = per[c] ?? {n: 0, nulls: 0};
      per[c].n++;
      if (val.isNone(i)) per[c].nulls++;
    }
    return Object.values(per).some((s) => s.n > 0 && s.nulls === s.n);
  });
}

test('Box Plot filter semantics and viewport response', async ({page}) => {
  test.setTimeout(600_000);

  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => { pageErrors.push(String(e)); });

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 5000});

  const rowCount = await page.evaluate(() => grok.shell.t.rowCount);
  expect(rowCount).toBe(100);
  await page.evaluate(() => {
    const tv = grok.shell.tv;
    const bp = tv.addViewer('Box plot');
    bp.props.valueColumnName = 'Average Mass';
    bp.props.category1ColumnName = 'Series';
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await v.waitForViewerRendered(page, 'Box plot', 1500);
  await v.waitForCanvasQuiet(page, 'Box plot');

  await openFilterPanelWithBudget(page);

  await softStep('Scenario 1 Step 4: narrowing the filter narrows the viewport; trueCount < 100; zoomValuesByFilter reads true', async () => {
    await resetFilter(page);
    await v.waitForCanvasQuiet(page, 'Box plot');
    expect(await filterTrueCount(page)).toBe(100);
    expect(await bpProp(page, 'zoomValuesByFilter')).toBe(true);
    const baseline = await viewportRect(page);

    let tc = rowCount;
    for (let attempt = 0; attempt < 2 && tc >= rowCount; attempt++) {
      const fc = await averageMassFilterCanvas(page);
      const minDragged = await dragFilterHandle(page, fc, 'min', 0.30);
      const maxDragged = await dragFilterHandle(page, fc, 'max', 0.70);
      console.log('S1 range-handle drags (min/max found):', minDragged, maxDragged);
      tc = await filterTrueCount(page);
    }
    await v.waitForCanvasQuiet(page, 'Box plot');
    const narrowed = await viewportRect(page);
    console.log('S1 viewport height baseline/narrowed:', baseline.height, narrowed.height, 'trueCount:', tc);
    expect(narrowed.height).toBeLessThan(baseline.height);
    expect(tc).toBeLessThan(100);
  });

  await softStep('Scenario 2 Step 6: with zoomValuesByFilter false and an active filter the viewport is unchanged', async () => {
    await setBpProp(page, 'zoomValuesByFilter', false);
    await resetFilter(page);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const baseline = await viewportRect(page);
    const tc = await narrowAverageMass(page, 300, 400);
    await v.waitForViewerRendered(page, 'Box plot', 1400);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const after = await viewportRect(page);
    console.log('S2 zoom-off viewport height baseline/after:', baseline.height, after.height, 'trueCount:', tc);
    expect(tc).toBeLessThan(100);
    expect(after.height).toBe(baseline.height);
    expect(after.y).toBe(baseline.y);
  });

  await softStep('Scenario 2 Step 8: switching zoomValuesByFilter back to true narrows the viewport to the filtered range', async () => {

    const before = await viewportRect(page);
    await setBpProp(page, 'zoomValuesByFilter', true);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const after = await viewportRect(page);
    console.log('S2 switch-back viewport height before/after:', before.height, after.height);
    expect(after.height).toBeLessThan(before.height);
  });

  await softStep('Scenario 3 Step 5: the viewer-local formula filter leaves df.filter.trueCount unchanged; the canvas repaints', async () => {
    await resetFilter(page);
    const tcBefore = await narrowAverageMass(page, 300, 400);
    await v.waitForViewerRendered(page, 'Box plot', 1200);
    await v.waitForCanvasQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');

    await setBpProp(page, 'filter', '${Average Mass} > 350', 1500);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const tcAfter = await filterTrueCount(page);
    const {deltaPx} = await v.diffCanvasColors(page, 'Box plot');
    console.log('S3 trueCount before/after:', tcBefore, tcAfter, 'canvas deltaPx:', deltaPx);
    expect(tcAfter).toBe(tcBefore);

    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeGreaterThan(1000);
    await setBpProp(page, 'filter', '', 800);
    await resetFilter(page);
  });

  let fixtureAdded = false;
  try {
    await softStep('Scenario 4 Step 5: disabling showEmptyCategories drops the empty-valued category (axis re-layout)', async () => {

      const fx = await page.evaluate(async () => {
        const df = grok.shell.t;
        const series = df.col('Series');
        const am = df.col('Average Mass');
        const per: Record<string, {n: number; nulls: number}> = {};
        for (let i = 0; i < df.rowCount; i++) {
          const c = String(series.get(i) ?? '');
          per[c] = per[c] ?? {n: 0, nulls: 0};
          per[c].n++;
          if (am.isNone(i)) per[c].nulls++;
        }
        if (Object.entries(per).some(([c, s]) => c.length > 0 && s.n > 0 && s.nulls === s.n))
          return {added: false};
        const cat = series.categories.find((c: string) => c && c.length > 0);
        await df.columns.addNewCalculated('AverageMassFixture',
          `if(\${Series} == "${cat}", null, \${Average Mass})`);
        return {added: true, category: cat};
      });
      fixtureAdded = fx.added;
      console.log('S4 empty-valued-category fixture:', JSON.stringify(fx));
      if (fx.added)
        await setBpProp(page, 'valueColumnName', 'AverageMassFixture', 1200);
      await setBpProp(page, 'showEmptyCategories', true, 1200);
      await v.waitForCanvasQuiet(page, 'Box plot');

      expect(await hasEmptyValuedCategory(page)).toBe(true);
      await v.snapshotCanvasColors(page, 'Box plot');
      const bandBefore = await categoryBandColors(page);
      await setBpProp(page, 'showEmptyCategories', false, 1400);
      await v.waitForCanvasQuiet(page, 'Box plot');
      const {deltaPx} = await v.diffCanvasColors(page, 'Box plot');
      const bandDelta = bandColorsDelta(bandBefore, await categoryBandColors(page));
      console.log('S4 re-layout deltaPx (off):', deltaPx, 'label-band delta:', bandDelta);
      expect(await bpProp(page, 'showEmptyCategories')).toBe(false);

      expect(deltaPx).toBeGreaterThan(4000);

      expect(bandDelta).toBeGreaterThan(100);
    });

    await softStep('Scenario 4 Step 7: re-enabling showEmptyCategories restores the empty-valued category (axis re-layout)', async () => {
      await v.waitForCanvasQuiet(page, 'Box plot');
      await v.snapshotCanvasColors(page, 'Box plot');
      const bandBefore = await categoryBandColors(page);
      await setBpProp(page, 'showEmptyCategories', true, 1400);
      await v.waitForCanvasQuiet(page, 'Box plot');
      const {deltaPx} = await v.diffCanvasColors(page, 'Box plot');
      const bandDelta = bandColorsDelta(bandBefore, await categoryBandColors(page));
      console.log('S4 re-layout deltaPx (back on):', deltaPx, 'label-band delta:', bandDelta);
      expect(await bpProp(page, 'showEmptyCategories')).toBe(true);

      expect(deltaPx).toBeGreaterThan(4000);

      expect(bandDelta).toBeGreaterThan(100);
    });
  } finally {
    if (fixtureAdded) {
      await setBpProp(page, 'valueColumnName', 'Average Mass', 1000).catch(() => {});
      await page.evaluate(() => { try { grok.shell.t.columns.remove('AverageMassFixture'); } catch (_) {} }).catch(() => {});
      await v.pollValue(() => page.evaluate(() => {
        try { return grok.shell.t.columns.names().includes('AverageMassFixture'); }
        catch (_) { return false; }
      }), (present) => !present, 500, 100);
    }
  }

  await softStep('Scenario 5 Step 4: setting Marker Color to TPSA with an active filter repaints without console errors (GROK-20171 floor)', async () => {
    await resetFilter(page);
    await narrowAverageMass(page, 300, 400);
    await v.waitForViewerRendered(page, 'Box plot', 1200);
    await v.waitForCanvasQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    const errBefore = consoleErrors.length + pageErrors.length;
    await setBpProp(page, 'markerColorColumnName', 'TPSA', 1500);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const {deltaPx} = await v.diffCanvasColors(page, 'Box plot');
    const errAfter = consoleErrors.length + pageErrors.length;

    const scale = await v.countCanvasPixels(page, 'Box plot', {canvasSelector: OVERLAY});
    console.log('S5 color repaint deltaPx:', deltaPx, 'overlay scale px:', JSON.stringify(scale),
      'errors before/after:', errBefore, errAfter);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('TPSA');
    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(scale.total).toBeGreaterThan(500);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 5 Step 5: color scale recomputes to the filtered TPSA range (overlay repaint)', async () => {

    const ranges = await page.evaluate(() => {
      const df = grok.shell.t;
      const t = df.col('TPSA');
      let mn = Infinity, mx = -Infinity;
      for (let i = 0; i < df.rowCount; i++) {
        if (!df.filter.get(i) || t.isNone(i)) continue;
        const x = t.get(i);
        if (x < mn) mn = x;
        if (x > mx) mx = x;
      }
      return {visMin: mn, visMax: mx, fullMin: t.stats.min, fullMax: t.stats.max};
    });
    console.log('S5 TPSA visible range vs full:', JSON.stringify(ranges));
    expect(ranges.visMin !== ranges.fullMin || ranges.visMax !== ranges.fullMax).toBe(true);
    await v.waitForCanvasQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot', OVERLAY);
    await resetFilter(page);

    await page.waitForTimeout(1200);
    await v.waitForCanvasQuiet(page, 'Box plot');
    const {deltaPx} = await v.diffCanvasColors(page, 'Box plot', OVERLAY);
    console.log('S5 overlay scale deltaPx on filter widen:', deltaPx);
    expect(deltaPx).toBeGreaterThan(1000);
    await setBpProp(page, 'markerColorColumnName', '', 600);
    await resetFilter(page);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

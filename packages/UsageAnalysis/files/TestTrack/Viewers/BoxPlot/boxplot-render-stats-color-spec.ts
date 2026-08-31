/* ---
realizes: [boxplot.cp.render-stats-color-sync, boxplot.int.violin-needs-bins-and-style, boxplot.int.pvalue-toggle-key-and-menu]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

async function bpProp(page: Page, prop: string): Promise<any> {
  return page.evaluate((p) => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    return bp?.props?.[p];
  }, prop);
}

async function setBpProp(page: Page, prop: string, value: any, settleMs = 900): Promise<void> {
  await v.setViewerProps(page, 'Box plot', [{set: {[prop]: value}, wait: settleMs}]);
}

async function canvasRect(page: Page): Promise<{x: number; y: number; w: number; h: number}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: c.x, y: c.y, w: c.width, h: c.height};
  });
}

async function revealGroupStatsIcon(page: Page): Promise<string> {
  const origin = await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: c.x, y: c.y};
  });
  for (const [dx, dy] of [[35, 15], [40, 17], [30, 14], [45, 16]]) {
    const shown1 = await v.armEvent(page, 'grok.events.onTooltipShown', 150);

    await page.mouse.move(origin.x + dx, origin.y + dy);
    await shown1();
  }

  return v.pollValue(() => page.evaluate(() => {
    const el = document.querySelector('[name="show-group-stats"]') as HTMLElement | null;
    return el ? getComputedStyle(el).visibility : 'absent';
  }), (vis) => vis === 'visible', 300, 100);
}

async function waitIconHidden(page: Page, capMs: number): Promise<string> {
  return v.pollValue(() => page.evaluate(() => {
    const el = document.querySelector('[name="show-group-stats"]') as HTMLElement | null;
    return el ? getComputedStyle(el).visibility : 'absent';
  }), (vis) => vis !== 'visible', capMs, 100);
}

async function colMax(page: Page, col: string): Promise<number> {
  return page.evaluate((c) => grok.shell.t.col(c).stats.max, col);
}

test('Box Plot rendering, statistics, and grid color synchronization', async ({page}) => {
  test.setTimeout(600_000);

  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => { pageErrors.push(String(e)); });

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const tv = grok.shell.tv;
    const bp = tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
    bp.props.category1ColumnName = 'SEX';
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await v.waitForViewerRendered(page, 'Box plot', 1500);
  await v.waitForViewerQuiet(page, 'Box plot');

  await softStep('[anchor: Scenario 1 Step 2] Box coloring baseline: whiskerColor=null gives per-category sequential hues', async () => {

    expect(await bpProp(page, 'whiskerColor')).toBeNull();
    const distinctHues = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      const cv = bp.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const colors = new Map<number, number>();
      for (let i = 0; i < img.length; i += 4) {
        const r = img[i], g = img[i + 1], b = img[i + 2], a = img[i + 3];
        if (a === 0 || (r >= 250 && g >= 250 && b >= 250)) continue;
        colors.set((r << 16) | (g << 8) | b, (colors.get((r << 16) | (g << 8) | b) ?? 0) + 1);
      }

      const top = [...colors.entries()].filter(([, n]) => n >= 500).map(([k]) => k);
      const bluish = top.some((k) => ((k >> 16) & 255) < ((k) & 255));      
      const orangish = top.some((k) => ((k >> 16) & 255) > ((k) & 255) + 8); 
      return {dominantCount: top.length, bluish, orangish};
    });
    console.log('Scenario 1 Step 2 dominant hue families:', JSON.stringify(distinctHues));
    expect(distinctHues.bluish).toBe(true);
    expect(distinctHues.orangish).toBe(true);
  });

  await softStep('[anchor: Scenario 1 Step 3] Resize wide so the statistics strip is visible (_statsValuesFit gate) with default stats present', async () => {

    expect(await bpProp(page, 'showStatistics')).toBe(true);
    const geom = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      const rect = bp.root.getBoundingClientRect();
      const closeStats = bp.root.querySelector('[name="close-stats"]');
      return {width: rect.width, hasCloseStats: !!closeStats};
    });
    console.log('Scenario 1 Step 3 viewer width / close-stats:', JSON.stringify(geom));
    expect(geom.width).toBeGreaterThan(300);
    expect(geom.hasCloseStats).toBe(true);
    expect(await bpProp(page, 'showAvg')).toBe(true);
    expect(await bpProp(page, 'showMed')).toBe(true);
  });

  await softStep('[anchor: Scenario 1 Step 6] Enable the statistics ladder: total/inliers/outliers/stdev/Q1/Q3 each read true, no console error', async () => {
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await setBpProp(page, 'showTotalCount', true, 300);
    await setBpProp(page, 'showInliersCount', true, 300);
    await setBpProp(page, 'showOutliersCount', true, 300);
    await setBpProp(page, 'showStdev', true, 300);
    await setBpProp(page, 'showQ1', true, 300);
    await setBpProp(page, 'showQ3', true, 500);
    await v.waitForViewerQuiet(page, 'Box plot');
    expect(await bpProp(page, 'showTotalCount')).toBe(true);
    expect(await bpProp(page, 'showInliersCount')).toBe(true);
    expect(await bpProp(page, 'showOutliersCount')).toBe(true);
    expect(await bpProp(page, 'showStdev')).toBe(true);
    expect(await bpProp(page, 'showQ1')).toBe(true);
    expect(await bpProp(page, 'showQ3')).toBe(true);
    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('Scenario 1 Step 6 error deltas:', JSON.stringify(errDelta), JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
  });

  await softStep('[anchor: Scenario 1 Step 8] Statistics Format sets and does not affect the p-value (statisticsFormat dependsOn showStatistics only)', async () => {
    await setBpProp(page, 'showPValue', true, 400);
    const pvBefore = await bpProp(page, 'showPValue');
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await setBpProp(page, 'statisticsFormat', '#,##0.00', 700);
    expect(await bpProp(page, 'statisticsFormat')).toBe('#,##0.00');
    expect(await bpProp(page, 'showPValue')).toBe(pvBefore);
    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('Scenario 1 Step 8 statisticsFormat error deltas:', JSON.stringify(errDelta), JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
    await setBpProp(page, 'statisticsFormat', 'auto', 400);
  });

  await softStep('[anchor: Scenario 1 Step 9] Explicit whiskerColor repaints the boxes to one uniform color (settle-gated canvas diff)', async () => {
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'whiskerColor', 0xFF1F77B4, 300);
    const deltaPx = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 2000, timeoutMs: 15000});
    console.log('Scenario 1 Step 9 whiskerColor canvas deltaPx:', deltaPx);

    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeGreaterThan(2000);
    expect(await bpProp(page, 'whiskerColor')).toBe(0xFF1F77B4);
    await setBpProp(page, 'whiskerColor', null, 600);
  });

  await softStep('[anchor: Scenario 2 Step 3] Bare p-value hover reveals the show-group-stats icon, no test-name tooltip', async () => {

    await setBpProp(page, 'showGroupComparison', false, 500);
    await setBpProp(page, 'showPValue', true, 600);
    expect(await bpProp(page, 'showGroupComparison')).toBe(false);
    const iconVis = await revealGroupStatsIcon(page);
    console.log('Scenario 2 Step 3 show-group-stats visibility on bare-p hover:', iconVis);
    expect(iconVis).toBe('visible');

    const tooltipText = await page.evaluate(() => {
      const tip = document.querySelector('.d4-tooltip');
      return (tip?.textContent ?? '').trim();
    });
    console.log('Scenario 2 Step 3 bare-p hover tooltip text:', JSON.stringify(tooltipText));
    expect(tooltipText).toBe('');
    const r = await canvasRect(page);
    await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
    await waitIconHidden(page, 300);
  });

  await softStep('[anchor: Scenario 2 Step 5] T key toggles Show P Value (transition signal; pvalue-toggle-key-and-menu)', async () => {

    const r = await canvasRect(page);
    await setBpProp(page, 'showPValue', false, 500);
    expect(await bpProp(page, 'showPValue')).toBe(false);
    await page.mouse.click(r.x + r.w * 0.5, r.y + r.h * 0.05);

    await v.waitForViewerRendered(page, 'Box plot', 300);
    await page.keyboard.press('t');
    await v.pollValue(() => bpProp(page, 'showPValue'), (val) => val === true, 700, 100);
    expect(await bpProp(page, 'showPValue')).toBe(true);
    await page.keyboard.press('t');
    await v.pollValue(() => bpProp(page, 'showPValue'), (val) => val === false, 700, 100);
    expect(await bpProp(page, 'showPValue')).toBe(false);
    await setBpProp(page, 'showPValue', true, 500);
  });

  await softStep('[anchor: Scenario 2 Step 7] Category1=RACE activates the 3+ category (Alexander-Govern) branch; bare p present', async () => {
    await setBpProp(page, 'category1ColumnName', 'RACE', 1200);

    await setBpProp(page, 'markerColorColumnName', 'SEX', 800);
    await setBpProp(page, 'showPValue', true, 600);
    const raceCats = await page.evaluate(() => grok.shell.t.col('RACE').categories.length);
    console.log('Scenario 2 Step 7 RACE category count:', raceCats);
    expect(raceCats).toBeGreaterThanOrEqual(3);
    expect(await bpProp(page, 'category1ColumnName')).toBe('RACE');
    const iconVis = await revealGroupStatsIcon(page);
    console.log('Scenario 2 Step 7 icon visibility (3+ cats):', iconVis);
    expect(iconVis).toBe('visible');
    const r = await canvasRect(page);
    await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
    await waitIconHidden(page, 300);
    await setBpProp(page, 'category1ColumnName', 'SEX', 1000);
  });

  await softStep('[anchor: Scenario 2 Step 9] Statistics-region right-click menu toggles the same prop as the panel (third toggle path)', async () => {
    await setBpProp(page, 'showStatistics', true, 500);
    const before = await bpProp(page, 'showTotalCount');

    const r = await canvasRect(page);
    let toggled = false;
    for (const fy of [0.93, 0.95, 0.91, 0.96]) {
      await page.mouse.click(r.x + r.w * 0.5, r.y + r.h * fy, {button: 'right'});

      await v.pollValue(() => page.evaluate(() =>
        (Array.from(document.querySelectorAll('.d4-menu-item')) as HTMLElement[]).some((i) => {
          const b = i.getBoundingClientRect();
          return b.width > 0 && b.height > 0
            && /^Show Total Count$/i.test((i.getAttribute('d4-name') ?? i.textContent ?? '').trim());
        })), (ready) => ready, 600, 100);
      const clicked = await page.evaluate(() => {
        const items = (Array.from(document.querySelectorAll('.d4-menu-item')) as HTMLElement[])
          .filter((i) => { const b = i.getBoundingClientRect(); return b.width > 0 && b.height > 0; });
        const label = (i: HTMLElement) => (i.getAttribute('d4-name') ?? i.textContent ?? '').trim();
        const item = items.find((i) => /^Show Total Count$/i.test(label(i)));
        if (!item) return false;
        item.click();
        return true;
      });
      await page.keyboard.press('Escape');

      await v.pollStable(() => page.evaluate(() => document.querySelectorAll('.d4-menu-popup').length),
        (a, b) => a === b, 400, 100);
      if (clicked) { toggled = true; break; }
    }
    console.log('Scenario 2 Step 9 menu toggle reached:', toggled);
    expect(toggled).toBe(true);
    expect(await bpProp(page, 'showTotalCount')).toBe(!before);
    await setBpProp(page, 'showTotalCount', before, 400);
  });

  await softStep('[anchor: Scenario 3 Step 3] Plot Style violin: large canvas diff and both SEX distributions present (github-2966)', async () => {
    await setBpProp(page, 'category1ColumnName', 'SEX', 800);
    await setBpProp(page, 'plotStyle', 'box', 600);
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'plotStyle', 'violin', 300);
    const deltaPx = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 5000, timeoutMs: 15000});
    console.log('Scenario 3 Step 3 box→violin deltaPx:', deltaPx);
    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeGreaterThan(5000);
    const strips = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      const cv = bp.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const w = cv.width, h = cv.height;
      const strip = (x0f: number, x1f: number) => {
        let n = 0;
        for (let y = 0; y < h; y++)
          for (let x = Math.floor(w * x0f); x < Math.floor(w * x1f); x++) {
            const i = (y * w + x) * 4;
            const r = img[i], g = img[i + 1], b = img[i + 2], a = img[i + 3];
            if (a !== 0 && !(r >= 250 && g >= 250 && b >= 250)) n++;
          }
        return n;
      };
      return {left: strip(0.25, 0.45), right: strip(0.55, 0.75)};
    });
    console.log('Scenario 3 Step 3 violin strip ink (left/right):', JSON.stringify(strips));

    expect(strips.left).toBeGreaterThan(0);
    expect(strips.right).toBeGreaterThan(0);
  });

  await softStep('[anchor: Scenario 3 Step 5] Bins 50→500 and Interquartile Line Width each yield a canvas diff (GROK-18245)', async () => {
    await setBpProp(page, 'bins', 50, 800);
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'bins', 500, 300);
    const binsDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Scenario 3 Step 5 bins 50→500 deltaPx:', binsDelta);
    expect(binsDelta).toBeGreaterThanOrEqual(0);
    expect(binsDelta).toBeGreaterThan(0);
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'interquartileLineWidth', 10, 300);
    const iqrDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Scenario 3 Step 5 interquartileLineWidth deltaPx:', iqrDelta);
    expect(iqrDelta).toBeGreaterThanOrEqual(0);
    expect(iqrDelta).toBeGreaterThan(0);
  });

  await softStep('[anchor: Scenario 3 Step 7] Violin Line Width and Violin Whisker Color each yield a distinct canvas diff', async () => {
    await setBpProp(page, 'plotStyle', 'box', 700);
    expect(await bpProp(page, 'plotStyle')).toBe('box');
    await setBpProp(page, 'plotStyle', 'violin', 900);
    expect(await bpProp(page, 'plotStyle')).toBe('violin');
    await setBpProp(page, 'violinLineWidth', 1, 700);
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'violinLineWidth', 4, 300);
    const vlwDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Scenario 3 Step 7 violinLineWidth deltaPx:', vlwDelta);
    expect(vlwDelta).toBeGreaterThanOrEqual(0);
    expect(vlwDelta).toBeGreaterThan(0);
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'violinWhiskerColor', 0xFF00AA00, 300);
    const vwcDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Scenario 3 Step 7 violinWhiskerColor deltaPx:', vwcDelta);
    expect(vwcDelta).toBeGreaterThanOrEqual(0);
    expect(vwcDelta).toBeGreaterThan(0);
  });

  await softStep('[anchor: Scenario 4 Step 4] Marker Color=WEIGHT with a linear grid scheme; min/max change repaints the box plot (GROK-17506)', async () => {
    await setBpProp(page, 'plotStyle', 'box', 600);
    await setBpProp(page, 'valueColumnName', 'AGE', 700);
    await setBpProp(page, 'category1ColumnName', 'SEX', 800);
    await page.evaluate(() => {
      const weight = grok.shell.t.col('WEIGHT');
      weight.meta.colors.setLinear([DG.Color.blue, DG.Color.red]);
    });
    await setBpProp(page, 'markerColorColumnName', 'WEIGHT', 1200);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('WEIGHT');
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await page.evaluate(() => {
      const weight = grok.shell.t.col('WEIGHT');
      weight.meta.colors.setLinear([DG.Color.blue, DG.Color.red], {min: 60, max: 120});
    });
    const deltaPx = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 2000, timeoutMs: 15000});
    console.log('Scenario 4 Step 4 grid min/max change deltaPx:', deltaPx);
    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeGreaterThan(2000);
  });

  await softStep('[anchor: Scenario 4 Step 6] Linear→conditional coloring repaints the box plot; selection pixels survive (github-3066)', async () => {
    const r = await canvasRect(page);
    await page.evaluate(() => grok.shell.t.selection.setAll(false));
    await v.waitForViewerRendered(page, 'Box plot', 300);
    await page.keyboard.down('Shift');
    await page.mouse.move(r.x + r.w * 0.55, r.y + r.h * 0.30);
    await page.mouse.down();
    await page.mouse.move(r.x + r.w * 0.72, r.y + r.h * 0.80, {steps: 14});
    await page.mouse.up();
    await page.keyboard.up('Shift');

    const selCount = await v.pollValue(() => page.evaluate(() => grok.shell.t.selection.trueCount),
      (n) => n > 0, 800, 100);
    console.log('Scenario 4 Step 6 selection count:', selCount);
    expect(selCount).toBeGreaterThan(0);
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await page.evaluate(() => {
      const weight = grok.shell.t.col('WEIGHT');
      weight.meta.colors.setConditional({'50-90': '#00ff00', '90-150': '#ffa500'});
    });
    const deltaPx = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 2000, timeoutMs: 15000});
    console.log('Scenario 4 Step 6 linear→conditional deltaPx:', deltaPx);
    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeGreaterThan(2000);
    const huePx = await v.countSelectionHuePixels(page, 'Box plot');
    console.log('Scenario 4 Step 6 selection-hue pixels under conditional coloring:', huePx);
    expect(huePx).toBeGreaterThan(0);
    await page.evaluate(() => {
      grok.shell.t.selection.setAll(false);
      grok.shell.t.col('WEIGHT').meta.colors.setDisabled();
    });
    await v.waitForViewerRendered(page, 'Box plot', 500);
  });

  await softStep('[anchor: Scenario 4 Step 8] Marker Color=SEX (categorical); changing a category color repaints the box plot (palette sync)', async () => {
    await setBpProp(page, 'markerColorColumnName', 'SEX', 1000);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('SEX');
    await v.waitForViewerQuiet(page, 'Box plot');
    await v.snapshotCanvasColors(page, 'Box plot');
    await page.evaluate(() => {
      const sex = grok.shell.t.col('SEX');
      const cats = sex.categories;
      const map: Record<string, string> = {};
      map[cats[0]] = '#e41a1c';
      sex.meta.colors.setCategorical(map, {fallbackColor: '#808080'});
    });
    const deltaPx = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 500, timeoutMs: 15000});
    console.log('Scenario 4 Step 8 categorical color change deltaPx:', deltaPx);
    expect(deltaPx).toBeGreaterThanOrEqual(0);
    expect(deltaPx).toBeGreaterThan(500);
    await page.evaluate(() => grok.shell.t.col('SEX').meta.colors.setDisabled());
    await setBpProp(page, 'markerColorColumnName', 'SEX', 500);
  });

  await softStep('[anchor: Scenario 5 Step 3] Value=STARTED (datetime): statistics render with no new console error (datetime statisticsFormat fix)', async () => {
    await setBpProp(page, 'showStatistics', true, 400);
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await setBpProp(page, 'valueColumnName', 'STARTED', 1500);
    expect(await bpProp(page, 'valueColumnName')).toBe('STARTED');
    await v.waitForViewerQuiet(page, 'Box plot');

    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('Scenario 5 Step 3 datetime error deltas:', JSON.stringify(errDelta), JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
    const errBefore2 = consoleErrors.length;
    const pageErrBefore2 = pageErrors.length;
    await setBpProp(page, 'statisticsFormat', 'yyyy-MM-dd HH:mm', 800);
    expect(await bpProp(page, 'statisticsFormat')).toBe('yyyy-MM-dd HH:mm');
    expect(consoleErrors.slice(errBefore2)).toEqual([]);
    expect(pageErrors.slice(pageErrBefore2)).toEqual([]);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});

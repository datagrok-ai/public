/* ---
realizes: [barchart.cp.setup-and-interact]
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitCol = 'Primary Series Name';
const dominantCategory = 'Triazoles';

function barBands(page: import('@playwright/test').Page) {
  return page.evaluate(() => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
    const data = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
    const rows: {count: number, minX: number, maxX: number}[] = [];
    for (let y = 0; y < cv.height; y++) {
      let count = 0, minX = cv.width, maxX = 0;
      for (let x = 0; x < cv.width; x++) {
        const i = (y * cv.width + x) * 4;
        const r = data[i], g = data[i + 1], b = data[i + 2];
        if (r >= 120 && r <= 180 && g >= 190 && g <= 240 && b >= 120 && b <= 180) {
          count++;
          if (x < minX) minX = x;
          if (x > maxX) maxX = x;
        }
      }
      rows.push({count, minX, maxX});
    }
    const minRun = 10; 
    const bands: {yTopFrac: number, yBottomFrac: number, yMidFrac: number, xMinFrac: number, xMidFrac: number}[] = [];
    let start = -1;
    for (let y = 0; y <= cv.height; y++) {
      const inBar = y < cv.height && rows[y].count >= minRun;
      if (inBar && start < 0) start = y;
      else if (!inBar && start >= 0) {
        const mid = Math.floor((start + y - 1) / 2);
        bands.push({
          yTopFrac: start / cv.height,
          yBottomFrac: (y - 1) / cv.height,
          yMidFrac: mid / cv.height,
          xMinFrac: rows[mid].minX / cv.width,
          xMidFrac: (rows[mid].minX + rows[mid].maxX) / 2 / cv.width,
        });
        start = -1;
      }
    }
    const rect = cv.getBoundingClientRect();
    return {bands, rect: {x: rect.x, y: rect.y, w: rect.width, h: rect.height}};
  });
}

test('Bar Chart — Setup and core interaction (On Click Filter vs Select)', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 4000});
  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await v.pollValue(() => page.locator('.property-grid').count(), (n) => n > 0, 500, 100);

  await v.setViewerProps(page, 'Bar chart', [{set: {
    splitColumnName: splitCol, valueColumnName: 'CAST Idea ID', valueAggrType: 'count',
  }, wait: 900}]);

  await softStep('Scenario 1 Step 1: per-category bars render for at least two distinct categories', async () => {
    const info = await page.evaluate(({split}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      const counts: Record<string, number> = {};
      for (const c of col.categories) counts[c] = 0;
      for (let i = 0; i < df.rowCount; i++) counts[col.get(i)]++;
      const nonEmpty = Object.values(counts).filter((n) => (n as number) > 0).length;
      return {
        split: bc.props.splitColumnName,
        value: bc.props.valueColumnName,
        aggr: bc.props.valueAggrType,
        nonEmptyCategories: nonEmpty,
        hasCanvas: !!bc.root.querySelector('canvas'),
        counts,
      };
    }, {split: splitCol});
    expect(info.split).toBe(splitCol);
    expect(info.aggr).toBe('count');
    expect(info.hasCanvas).toBe(true);
    expect(info.nonEmptyCategories).toBeGreaterThanOrEqual(2);
    expect(info.counts[dominantCategory]).toBeGreaterThan(0);
  });

  await softStep('Scenario 1 Step 4: On Click=Filter — clicking a bar filters the grid to that category', async () => {
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.onClick = 'Filter';
      grok.shell.tv.dataFrame.filter.setAll(true);
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} } 
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);
    const geo = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv.getBoundingClientRect();
      return {
        onClick: bc.props.onClick,
        rect: {x: rect.x, y: rect.y, w: rect.width, h: rect.height},
        before: grok.shell.tv.dataFrame.filter.trueCount,
        rowCount: grok.shell.tv.dataFrame.rowCount,
      };
    });
    expect(geo.onClick).toBe('Filter');

    const px = geo.rect.x + geo.rect.w * 0.55;
    const py = geo.rect.y + geo.rect.h * 0.392;
    await page.mouse.click(px, py);
    const probe = await v.pollValue(() => page.evaluate(({split, cat}) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      let survivorsForCat = 0;
      const catsSeen = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) {
        if (df.filter.get(i)) {
          catsSeen.add(col.get(i));
          if (col.get(i) === cat) survivorsForCat++;
        }
      }
      return {
        trueCount: df.filter.trueCount,
        rowCount: df.rowCount,
        distinctCats: catsSeen.size,
        survivorsForCat,
      };
    }, {split: splitCol, cat: dominantCategory}),
    (p) => p.trueCount > 0 && p.trueCount < p.rowCount && p.distinctCats === 1, 3000, 150);
    expect(probe.trueCount).toBeLessThan(probe.rowCount);
    expect(probe.trueCount).toBeGreaterThan(0);
    expect(probe.distinctCats).toBe(1);
    expect(probe.survivorsForCat).toBe(probe.trueCount);
  });

  await softStep('Scenario 1 Step 5: clicking blank canvas space clears the click-filter (grid returns to full range)', async () => {

    await page.evaluate(() => {
      grok.shell.tv.dataFrame.filter.setAll(true);
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} } 
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const geo = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv.getBoundingClientRect();
      return {onClick: bc.props.onClick, rect: {x: rect.x, y: rect.y, w: rect.width, h: rect.height}};
    });
    expect(geo.onClick).toBe('Filter');
    const filterState = () => page.evaluate(({split}) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      const cats = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) cats.add(col.get(i));
      return {trueCount: df.filter.trueCount, distinctCats: cats.size, rowCount: df.rowCount};
    }, {split: splitCol});
    await page.mouse.click(geo.rect.x + geo.rect.w * 0.55, geo.rect.y + geo.rect.h * 0.392);
    const filtered = await v.pollValue(filterState,
      (f) => f.distinctCats === 1 && f.trueCount < f.rowCount, 3000, 150);
    expect(filtered.distinctCats).toBe(1);
    expect(filtered.trueCount).toBeLessThan(filtered.rowCount);
    await page.mouse.click(geo.rect.x + geo.rect.w * 0.55, geo.rect.y + geo.rect.h * 0.6);
    const afterBlank = await v.pollValue(filterState,
      (a) => a.trueCount === a.rowCount && a.distinctCats > filtered.distinctCats, 3000, 150);

    expect(afterBlank.trueCount).toBe(afterBlank.rowCount);
    expect(afterBlank.distinctCats).toBeGreaterThan(filtered.distinctCats);
  });

  await softStep('Scenario 1 Step 9: Double-click canvas (Reset View) leaves the chart error-free', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    await page.evaluate(() => {
      grok.shell.tv.dataFrame.filter.setAll(true);
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} } 
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const rect = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const r = cv.getBoundingClientRect();
      return {x: r.x, y: r.y, w: r.width, h: r.height};
    });
    await page.mouse.dblclick(rect.x + rect.w * 0.55, rect.y + rect.h * 0.5);

    await page.waitForTimeout(700);
    const state = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      return {hasCanvas: !!bc.root.querySelector('canvas'), split: bc.props.splitColumnName};
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(state.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 3: On Click=Select — clicking a bar selects rows without filtering', async () => {
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.onClick = 'Select';
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      df.selection.setAll(false);
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const geo = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv.getBoundingClientRect();
      return {
        onClick: bc.props.onClick,
        rect: {x: rect.x, y: rect.y, w: rect.width, h: rect.height},
        rowCount: grok.shell.tv.dataFrame.rowCount,
      };
    });
    expect(geo.onClick).toBe('Select');

    await page.mouse.click(geo.rect.x + geo.rect.w * 0.55, geo.rect.y + geo.rect.h * 0.392);
    const probe = await v.pollValue(() => page.evaluate(({split, cat}) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      let selForCat = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.selection.get(i) && col.get(i) === cat) selForCat++;
      return {sel: df.selection.trueCount, selForCat, filt: df.filter.trueCount, rowCount: df.rowCount};
    }, {split: splitCol, cat: dominantCategory}),
    (p) => p.sel > 0 && p.selForCat === p.sel, 3000, 150);
    expect(probe.sel).toBeGreaterThan(0);
    expect(probe.selForCat).toBe(probe.sel);
    expect(probe.filt).toBe(probe.rowCount);
  });

  await softStep('Scenario 2 Step 5: pressing Esc clears the selection', async () => {

    const selBefore = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selBefore).toBeGreaterThan(0);

    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLElement;
      cv?.focus?.();
    });
    await page.keyboard.press('Escape');
    const selAfter = await v.pollValue(
      () => page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount as number),
      (n) => n === 0, 3000, 150);
    expect(selAfter).toBe(0);
  });

  const selState = () => page.evaluate(({split}) => {
    const df = grok.shell.tv.dataFrame;
    const col = df.col(split);
    const cats = new Set<string>();
    for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) cats.add(col.get(i));
    let catTotal = 0;
    for (let i = 0; i < df.rowCount; i++) if (cats.has(col.get(i))) catTotal++;
    return {sel: df.selection.trueCount, distinctCats: cats.size, catTotal, filt: df.filter.trueCount, rowCount: df.rowCount};
  }, {split: splitCol});

  await softStep('Scenario 2 Step 6: Ctrl+click adds a second category to the selection (additive, not replacement)', async () => {

    await page.mouse.move(5, 5);
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.onClick = 'Select';
      const df = grok.shell.tv.dataFrame;
      df.rows.highlight(null);
      df.filter.setAll(true);
      df.selection.setAll(false);
    });

    const {bands, rect} = await v.pollValue(() => barBands(page), (b) => b.bands.length >= 2, 3000, 150);
    expect(bands.length).toBeGreaterThanOrEqual(2);

    const ctrlClickBand = async (band: {yMidFrac: number, xMidFrac: number}) => {
      await page.keyboard.down('Control');
      await page.mouse.click(rect.x + band.xMidFrac * rect.w, rect.y + band.yMidFrac * rect.h);
      await page.keyboard.up('Control');
    };

    await ctrlClickBand(bands[0]);
    const s1 = await v.pollValue(selState, (s) => s.sel > 0 && s.distinctCats === 1, 3000, 150);
    expect(s1.sel).toBeGreaterThan(0);
    expect(s1.distinctCats).toBe(1);
    expect(s1.sel).toBe(s1.catTotal);
    expect(s1.filt).toBe(s1.rowCount);

    await ctrlClickBand(bands[1]);
    const s2 = await v.pollValue(selState, (s) => s.distinctCats === 2 && s.sel > s1.sel, 3000, 150);
    expect(s2.sel).toBeGreaterThan(s1.sel);
    expect(s2.distinctCats).toBe(2);
    expect(s2.sel).toBe(s2.catTotal);
    expect(s2.filt).toBe(s2.rowCount);
  });

  await softStep('Scenario 2 Step 8: Shift+drag rubber-band over two bars selects the covered categories', async () => {

    await page.mouse.move(5, 5);
    await page.evaluate(() => {
      grok.shell.tv.dataFrame.rows.highlight(null);
      grok.shell.tv.dataFrame.selection.setAll(false);
    });
    const {bands, rect} = await v.pollValue(() => barBands(page), (b) => b.bands.length >= 2, 3000, 150);
    expect(bands.length).toBeGreaterThanOrEqual(2);

    const xStart = rect.x + Math.max(bands[0].xMinFrac, bands[1].xMinFrac) * rect.w + 3;
    const xEnd = xStart + 24;
    const yStart = rect.y + bands[0].yTopFrac * rect.h + 1;
    const yEnd = rect.y + bands[1].yBottomFrac * rect.h - 1;
    await page.mouse.move(xStart, yStart);
    await page.keyboard.down('Shift');
    await page.mouse.down();
    await page.mouse.move(xEnd, yEnd, {steps: 8});
    await page.mouse.up();
    await page.keyboard.up('Shift');

    const s = await v.pollValue(selState, (x) => x.sel > 0 && x.distinctCats >= 2, 3000, 150);
    expect(s.sel).toBeGreaterThan(0);
    expect(s.distinctCats).toBeGreaterThanOrEqual(2);
    expect(s.sel).toBe(s.catTotal);
    expect(s.filt).toBe(s.rowCount);
  });

  await softStep('Scenario 2 Step 9: Esc clears the gesture-built selection (round-trip to zero)', async () => {

    const selBefore = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selBefore).toBeGreaterThan(0);

    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLElement;
      cv?.focus?.();
    });
    await page.keyboard.press('Escape');
    const selAfter = await v.pollValue(
      () => page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount as number),
      (n) => n === 0, 3000, 150);
    expect(selAfter).toBe(0);
  });

  await softStep('Scenario 3 Step 1: grid color coding on the Split column drives bar colors and survives a layout round-trip', async () => {

    const CEIL = 250;
    const FLOOR = 8000;

    await page.evaluate(({split}) => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      df.selection.setAll(false);
      const col = df.col(split);
      delete col.tags['.color-coding-categorical'];
      delete col.tags['.color-coding-type'];
      col.meta.colors.setCategorical({});
    }, {split: splitCol});

    await page.waitForTimeout(500);

    await v.snapshotCanvasColors(page, 'Bar chart');

    await page.waitForTimeout(500);
    const precheck = (await v.diffCanvasColors(page, 'Bar chart')).deltaPx;

    const scheme = await page.evaluate(({split}) => {
      const col = grok.shell.tv.dataFrame.col(split);
      const palette = ['#ff0000', '#00ff00', '#0000ff', '#ff00ff', '#00ffff', '#ffff00', '#ff8000', '#8000ff'];
      const s: Record<string, string> = {};
      col.categories.forEach((c: string, i: number) => { s[c] = palette[i % palette.length]; });
      col.meta.colors.setCategorical(s);
      for (const vw of grok.shell.tv.viewers) if (vw.type !== 'Grid') try { vw.invalidate?.(); } catch (_) {}
      return s;
    }, {split: splitCol});
    await v.waitForViewerRendered(page, 'Bar chart', 800);
    const colorDelta = (await v.diffCanvasColors(page, 'Bar chart')).deltaPx;

    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id;
    });

    await page.waitForTimeout(800);
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.close();
    });
    await v.pollValue(
      () => page.evaluate(() => !Array.from(grok.shell.tv.viewers).some((x: any) => x.type === 'Bar chart')),
      (closed) => closed, 500, 100);
    await page.evaluate(({split}) => {
      const col = grok.shell.tv.dataFrame.col(split);
      delete col.tags['.color-coding-categorical'];
      delete col.tags['.color-coding-type'];
      col.meta.colors.setCategorical({});
    }, {split: splitCol});
    const clearedKeys = await v.pollValue(() => page.evaluate(({split}) => {
      const col = grok.shell.tv.dataFrame.col(split);
      return Object.keys(JSON.parse(col.tags['.color-coding-categorical'] ?? '{}')).length;
    }, {split: splitCol}), (n) => n === 0, 400, 100);
    await page.evaluate(async ({id}) => {
      const saved = await grok.dapi.layouts.find(id);
      (window as any).__savedLayout = saved;
      grok.shell.tv.loadLayout(saved);
    }, {id: layoutId});
    const rt = await v.pollValue(() => page.evaluate(({split}) => {
      const col2 = grok.shell.tv.dataFrame.col(split);
      const bc2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      return {restored: JSON.parse(col2.tags['.color-coding-categorical'] ?? '{}'), reopened: !!bc2};
    }, {split: splitCol}),
    (x) => x.reopened && Object.keys(x.restored).length > 0, 3000, 150);
    await page.evaluate(async () => {
      const w = window as any;
      await grok.dapi.layouts.delete(w.__savedLayout);
      delete w.__savedLayout;
    });

    expect(precheck).toBeGreaterThanOrEqual(0); 
    expect(precheck).toBeLessThan(CEIL);
    expect(colorDelta).toBeGreaterThan(FLOOR);
    expect(clearedKeys).toBe(0);
    expect(rt.restored).toEqual(scheme);
    expect(rt.reopened).toBe(true);
  });

  v.finishSpec();
});

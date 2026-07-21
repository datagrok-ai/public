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
  await page.waitForTimeout(500);

  await page.evaluate(async ({split}) => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    bc.props.splitColumnName = split;
    bc.props.valueColumnName = 'CAST Idea ID';
    bc.props.valueAggrType = 'count';
    await new Promise((r) => setTimeout(r, 900));
  }, {split: splitCol});

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
    const geo = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.onClick = 'Filter';
      grok.shell.tv.dataFrame.filter.setAll(true);
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} } // best-effort cleanup
      await new Promise((r) => setTimeout(r, 600));
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
    await page.waitForTimeout(800);
    const probe = await page.evaluate(({split, cat}) => {
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
    }, {split: splitCol, cat: dominantCategory});
    expect(probe.trueCount).toBeLessThan(probe.rowCount);
    expect(probe.trueCount).toBeGreaterThan(0);
    expect(probe.distinctCats).toBe(1);
    expect(probe.survivorsForCat).toBe(probe.trueCount);
  });

  await softStep('Scenario 1 Step 5: clicking blank canvas space clears the click-filter (grid returns to full range)', async () => {
    // With On Click=Filter still active, re-apply the dominant-bar filter, then
    // click an empty canvas zone. Only the dominant bar (Triazoles) is
    // hit-testable headless (recon 2026-07-21: every other canvas position
    // registers no bar hit), so the "switch to a different bar" leg is a
    // documented reduction; this asserts the reachable half — a blank-zone
    // click releases the click-filter and the grid returns to the full range
    // (observed dev behavior 2026-07-21: tc 64 → 100, distinct 1 → 5).
    const geo = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      grok.shell.tv.dataFrame.filter.setAll(true);
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} } // best-effort cleanup
      await new Promise((r) => setTimeout(r, 500));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv.getBoundingClientRect();
      return {onClick: bc.props.onClick, rect: {x: rect.x, y: rect.y, w: rect.width, h: rect.height}};
    });
    expect(geo.onClick).toBe('Filter');
    await page.mouse.click(geo.rect.x + geo.rect.w * 0.55, geo.rect.y + geo.rect.h * 0.392);
    await page.waitForTimeout(700);
    const filtered = await page.evaluate(({split}) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      const cats = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) cats.add(col.get(i));
      return {trueCount: df.filter.trueCount, distinctCats: cats.size, rowCount: df.rowCount};
    }, {split: splitCol});
    expect(filtered.distinctCats).toBe(1);
    expect(filtered.trueCount).toBeLessThan(filtered.rowCount);
    await page.mouse.click(geo.rect.x + geo.rect.w * 0.55, geo.rect.y + geo.rect.h * 0.6);
    await page.waitForTimeout(700);
    const afterBlank = await page.evaluate(({split}) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      const cats = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) cats.add(col.get(i));
      return {trueCount: df.filter.trueCount, distinctCats: cats.size, rowCount: df.rowCount};
    }, {split: splitCol});
    // Blank-canvas click under On Click=Filter clears the active bar filter and
    // the grid returns to the full range (probe 2026-07-21: tc 64 → 100,
    // distinct 1 → 5). This is the reachable half of md Scenario 1 Step 5.
    expect(afterBlank.trueCount).toBe(afterBlank.rowCount);
    expect(afterBlank.distinctCats).toBeGreaterThan(filtered.distinctCats);
  });

  await softStep('Scenario 1 Step 9: Double-click canvas (Reset View) leaves the chart error-free', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    // md Scenario 1 Step 8 (Alt+drag to zoom into a range) is a documented
    // reduction: a headless Alt+drag over the canvas produces no zoom (recon
    // 2026-07-21: zero canvas delta before/after), so Reset View cannot be
    // preceded by a real zoom here. This exercises Reset View directly and
    // asserts it stays error-free with the chart intact.
    await page.evaluate(async () => {
      grok.shell.tv.dataFrame.filter.setAll(true);
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} } // best-effort cleanup
      await new Promise((r) => setTimeout(r, 500));
    });
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
    const geo = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.onClick = 'Select';
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      df.selection.setAll(false);
      await new Promise((r) => setTimeout(r, 500));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv.getBoundingClientRect();
      return {
        onClick: bc.props.onClick,
        rect: {x: rect.x, y: rect.y, w: rect.width, h: rect.height},
        rowCount: df.rowCount,
      };
    });
    expect(geo.onClick).toBe('Select');

    await page.mouse.click(geo.rect.x + geo.rect.w * 0.55, geo.rect.y + geo.rect.h * 0.392);
    await page.waitForTimeout(700);
    const probe = await page.evaluate(({split, cat}) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      let selForCat = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.selection.get(i) && col.get(i) === cat) selForCat++;
      return {sel: df.selection.trueCount, selForCat, filt: df.filter.trueCount, rowCount: df.rowCount};
    }, {split: splitCol, cat: dominantCategory});
    expect(probe.sel).toBeGreaterThan(0);
    expect(probe.selForCat).toBe(probe.sel);
    expect(probe.filt).toBe(probe.rowCount);
  });

  await softStep('Scenario 2 Step 5: pressing Esc clears the selection', async () => {
    const selBefore = await page.evaluate(({split, cat}) => {
      const df = grok.shell.tv.dataFrame;
      if (df.selection.trueCount === 0) {
        const col = df.col(split);
        for (let i = 0; i < df.rowCount; i++) if (col.get(i) === cat) df.selection.set(i, true);
      }
      return df.selection.trueCount;
    }, {split: splitCol, cat: dominantCategory});
    expect(selBefore).toBeGreaterThan(0);

    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLElement;
      cv?.focus?.();
    });
    await page.keyboard.press('Escape');
    await page.waitForTimeout(600);
    const selAfter = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selAfter).toBe(0);
  });

  v.finishSpec();
});

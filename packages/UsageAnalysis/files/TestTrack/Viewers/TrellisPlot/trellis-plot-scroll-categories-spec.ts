/* ---
realizes: [trellisplot.cp.scroll-categories, trellisplot.int.pack-categories-vs-viewport-scroll]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text);

async function xAxisState(page: Page): Promise<{
  cells: number; plusEnabled: boolean; minusEnabled: boolean;
}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const icons = root.querySelector('[name="x-axis-icons"]');
    const plus = icons?.querySelector('[name="icon-plus"]') as HTMLElement | null;
    const minus = icons?.querySelector('[name="icon-minus"]') as HTMLElement | null;
    const enabled = (el: HTMLElement | null) =>
      !!el && getComputedStyle(el).pointerEvents !== 'none';
    return {
      cells: root.querySelectorAll('.d4-trellis-plot-cell').length,
      plusEnabled: enabled(plus),
      minusEnabled: enabled(minus),
    };
  });
}

async function clickXIcon(page: Page, which: 'plus' | 'minus'): Promise<void> {
  await page.evaluate((w) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    (root.querySelector(`[name="x-axis-icons"] [name="icon-${w}"]`) as HTMLElement)?.click();
  }, which);
  await v.waitForViewerRendered(page, 'Trellis plot', 800);
}

async function realClickXIcon(page: Page, which: 'plus' | 'minus', timeout: number): Promise<boolean> {
  try {
    await page.locator(`[name="viewer-Trellis-plot"] [name="x-axis-icons"] [name="icon-${which}"]`)
      .click({timeout});
  }
  catch {
    return false;
  }
  await v.waitForViewerRendered(page, 'Trellis plot', 800);
  return true;
}

async function resetXColumnsViaMenu(page: Page): Promise<boolean> {
  const clicked = await page.evaluate(async () => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const combos = Array.from(root.querySelectorAll('[name="div-column-combobox-"]'));

    const xCombo = combos.find((c) => {
      const r = c.getBoundingClientRect();
      return r.width > r.height;
    }) as HTMLElement | undefined;
    if (!xCombo) return false;
    const r = xCombo.getBoundingClientRect();
    xCombo.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, button: 2,
    }));

    await new Promise((res) => setTimeout(res, 600));
    const label = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
      .find((el) => {
        if ((el.textContent || '').trim() !== 'Reset X columns') return false;
        const box = el.getBoundingClientRect();
        return box.width > 0 && box.height > 0;
      });
    if (!label) { document.body.click(); return false; }
    (label.closest('.d4-menu-item') as HTMLElement).click();
    return true;
  });
  await v.waitForViewerRendered(page, 'Trellis plot', 900);
  return clicked;
}

test('Trellis plot: category viewport paging (+/- icons and pack-categories coupling)', async ({page}) => {
  test.setTimeout(600_000);
  page.setDefaultTimeout(120_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });

  await loginToDatagrok(page);
  await page.waitForTimeout(5000); 

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = true; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  const cardinalities = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {disPop: df.col('DIS_POP').categories.length, race: df.col('RACE').categories.length};
  });
  expect(cardinalities).toEqual({disPop: 6, race: 4});

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');

  await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.xColumnNames = ['DIS_POP'];
    tp.props.yColumnNames = ['RACE'];
    tp.props.packCategories = true;
  });
  await v.waitForViewerRendered(page, 'Trellis plot', 900);

  const yCatCount = await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    return tp.yCategoriesCount as number;
  });
  expect(yCatCount).toBe(4);

  const entry = await xAxisState(page);
  expect(entry.plusEnabled).toBe(false);
  expect(entry.minusEnabled).toBe(true);

  const pagingWindow = {pageErrors: pageErrors.length, consoleErrors: consoleErrors.length};

  await softStep('Scenario 1 Step 4', async () => {

    await clickXIcon(page, 'minus');
    const paged = await xAxisState(page);
    expect(paged.plusEnabled).toBe(true);
    const cellsBeforePlus = paged.cells;

    await clickXIcon(page, 'plus');
    const afterPlus = await xAxisState(page);
    expect(afterPlus.cells).toBe(cellsBeforePlus + yCatCount);
  });

  await softStep('Scenario 1 Step 5', async () => {

    await clickXIcon(page, 'minus');
    const preExpansion = (await xAxisState(page)).cells;
    await clickXIcon(page, 'plus');
    const expanded = (await xAxisState(page)).cells;
    expect(expanded).toBe(preExpansion + yCatCount);
    await clickXIcon(page, 'minus');
    const afterMinus = (await xAxisState(page)).cells;
    expect(afterMinus).toBe(preExpansion);
  });

  await softStep('Scenario 1 Step 6', async () => {

    for (let i = 0; i < 8; i++) {
      const s = await xAxisState(page);
      if (!s.minusEnabled) break;
      await clickXIcon(page, 'minus');
    }
    const atMin = await xAxisState(page);
    expect(atMin.cells).toBe(yCatCount); 
    expect(atMin.minusEnabled).toBe(false);
    expect(atMin.plusEnabled).toBe(true);

    expect(await realClickXIcon(page, 'minus', 5000)).toBe(false);
    expect((await xAxisState(page)).cells).toBe(atMin.cells);
    expect(await realClickXIcon(page, 'plus', 20000)).toBe(true);
    expect((await xAxisState(page)).cells).toBe(atMin.cells + yCatCount);

    for (let i = 0; i < 8; i++) {
      const s = await xAxisState(page);
      if (!s.plusEnabled) break;
      await clickXIcon(page, 'plus');
    }
    const atMax = await xAxisState(page);
    expect(atMax.cells).toBe(cardinalities.disPop * yCatCount); 
    expect(atMax.plusEnabled).toBe(false);
    expect(atMax.minusEnabled).toBe(true);

    expect(await realClickXIcon(page, 'plus', 5000)).toBe(false);
    expect((await xAxisState(page)).cells).toBe(atMax.cells);
    expect(await realClickXIcon(page, 'minus', 20000)).toBe(true);
    expect((await xAxisState(page)).cells).toBe(atMax.cells - yCatCount);
  });

  await softStep('Scenario 1 Step 7', async () => {

    expect(pageErrors.slice(pagingWindow.pageErrors)).toEqual([]);
    expect(consoleErrors.slice(pagingWindow.consoleErrors)).toEqual([]);
  });

  await softStep('Scenario 1 Step 8', async () => {

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX', 'DIS_POP'];
      tp.props.yColumnNames = ['RACE'];
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const pagedCells = (await xAxisState(page)).cells;
    const pagedXCat = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return tp.xCategoriesCount as number;
    });
    expect(pagedCells).toBeLessThan(pagedXCat * yCatCount);

    const errorsBefore = pageErrors.length + consoleErrors.length;
    const reset = await resetXColumnsViaMenu(page);
    expect(reset).toBe(true);

    const after = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      return {
        xNames: tp.props.xColumnNames as string[],
        xCat: tp.xCategoriesCount as number,
        cells: root.querySelectorAll('.d4-trellis-plot-cell').length,
      };
    });

    expect(after.xNames).toEqual([]);
    expect(after.cells).toBe(yCatCount);
    const errorsAfter = pageErrors.length + consoleErrors.length;

    expect(errorsAfter).toBe(errorsBefore);

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['DIS_POP'];
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
  });

  let packedMaxCells = -1;
  let unpackedMaxCells = -1;

  await softStep('Scenario 2 Step 4', async () => {

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['RACE', 'SEVERITY'];
      tp.props.yColumnNames = ['SEX'];
      tp.props.packCategories = true;
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const truth = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE'), severity = df.col('SEVERITY');
      const seen = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) seen.add(`${race.get(i)}|${severity.get(i)}`);
      return {
        full: race.categories.length * severity.categories.length,
        populated: seen.size,
        sexCats: df.col('SEX').categories.length as number,
        xCat: tp.xCategoriesCount as number,
        yCat: tp.yCategoriesCount as number,
      };
    });

    expect(truth.xCat).toBe(truth.full);
    expect(truth.yCat).toBe(truth.sexCats);
    expect(truth.populated).toBeLessThan(truth.full);

    const paged = await xAxisState(page);
    expect(paged.plusEnabled).toBe(true);
    expect(paged.cells % truth.yCat).toBe(0); 
    await clickXIcon(page, 'plus');
    const afterPlus = await xAxisState(page);

    expect(afterPlus.cells).toBe(paged.cells + truth.yCat);

    for (let i = 0; i < 30; i++) {
      const s = await xAxisState(page);
      if (!s.plusEnabled) break;
      await clickXIcon(page, 'plus');
    }
    const atMax = await xAxisState(page);
    expect(atMax.plusEnabled).toBe(false);
    expect(atMax.cells % truth.yCat).toBe(0);
    expect(atMax.cells).toBeGreaterThan(0);

    expect(atMax.cells).toBeLessThanOrEqual(truth.populated * truth.yCat);
    packedMaxCells = atMax.cells;
    unpackedMaxCells = truth.full * truth.yCat;
  });

  await softStep('Scenario 2 Step 6', async () => {

    expect(packedMaxCells).toBeGreaterThan(0);
    const errorsBefore = pageErrors.length + consoleErrors.length;

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.packCategories = false;
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);

    for (let i = 0; i < 30; i++) {
      const s = await xAxisState(page);
      if (!s.plusEnabled) break;
      await clickXIcon(page, 'plus');
    }
    const offMax = await xAxisState(page);
    expect(offMax.plusEnabled).toBe(false);

    expect(offMax.cells).toBe(unpackedMaxCells);

    expect(offMax.cells).toBeGreaterThan(packedMaxCells);

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.packCategories = true;
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const errorsAfter = pageErrors.length + consoleErrors.length;
    expect(errorsAfter).toBe(errorsBefore);
  });

  v.finishSpec();
});

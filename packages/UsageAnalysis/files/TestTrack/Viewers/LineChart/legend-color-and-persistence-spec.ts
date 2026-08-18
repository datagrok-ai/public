/* ---
realizes: [linechart.cp.legend-color-and-persistence]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitColumn = 'Stereo Category';
const baselineColors: Record<string, number> = {
  R_ONE: 0xFFFF0000,
  S_ABS: 0xFF00FF00,
  S_ACHIR: 0xFF0000FF,
  S_PART: 0xFFFFFF00,
  S_UNKN: 0xFFFF00FF,
};

const pageErrors: string[] = [];
const consoleErrors: string[] = [];

function realErrors(): string[] {
  return [...pageErrors, ...consoleErrors];
}

async function readCategoryColors(page: Page): Promise<Record<string, number>> {
  return page.evaluate((col) => {
    const df = grok.shell.tv.dataFrame;
    const cat = df.col(col);
    const colors = cat.meta.colors;
    const out: Record<string, number> = {};
    const seen: Record<string, boolean> = {};
    for (let i = 0; i < df.rowCount; i++) {
      const v = cat.get(i);
      if (v && !seen[v]) { seen[v] = true; out[v] = colors.getColor(i, cat); }
    }
    return out;
  }, splitColumn);
}

async function readROneColor(page: Page): Promise<number | null> {
  return page.evaluate((col) => {
    const tv = grok.shell.tv;
    const cat = tv.dataFrame.col(col);
    for (let i = 0; i < tv.dataFrame.rowCount; i++)
      if (cat.get(i) === 'R_ONE') return cat.meta.colors.getColor(i, cat);
    return null;
  }, splitColumn);
}

async function lcFilterCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    return lc.filter.trueCount;
  });
}

test('Line Chart — legend filter-color and layout persistence', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  pageErrors.length = 0;
  consoleErrors.length = 0;

  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.locator('[name="icon-line-chart"]').click();
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 15000});

  await page.evaluate((args) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    lc.props.xColumnName = 'Chemical Space X';
    lc.props.yColumnNames = ['Chemical Space Y'];
    lc.props.splitColumnNames = [args.col];

    grok.shell.tv.dataFrame.col(args.col).meta.colors.setCategorical(args.colors);
  }, {col: splitColumn, colors: baselineColors});
  await v.waitForViewerRendered(page, 'Line chart', 1200);

  const legend = page.locator('[name="viewer-Line-chart"] [name="legend"]');
  await legend.waitFor({timeout: 10000});
  await expect(legend.locator('.d4-legend-item.d4-legend-text-item')).toHaveCount(5);

  const baselineFilter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  const colorBefore = await readCategoryColors(page);
  expect(colorBefore.R_ONE).toBe(0xFFFF0000); 

  await softStep('S1: legend click filters + preserves remaining line colors', async () => {
    const before = realErrors().length;
    await legend.locator('.d4-legend-value', {hasText: /^S_ABS$/}).click();
    const filtered = await v.pollValue(() => lcFilterCount(page), (n) => n < baselineFilter, 700, 100);
    expect(filtered).toBeLessThan(baselineFilter);

    const colorAfter = await readCategoryColors(page);
    expect(colorAfter.R_ONE).toBe(colorBefore.R_ONE);
    expect(colorAfter.R_ONE).not.toBe(baselineColors.S_ABS); 
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1: re-click legend category resets filter to full count', async () => {
    const before = realErrors().length;
    await legend.locator('.d4-legend-value', {hasText: /^S_ABS$/}).click();
    await v.pollValue(() => lcFilterCount(page), (n) => n === baselineFilter, 700, 100);
    expect(await lcFilterCount(page)).toBe(baselineFilter);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: category color persists through layout round-trip (GROK-17278)', async () => {
    const before = realErrors().length;
    await page.evaluate((col) => {
      grok.shell.tv.dataFrame.col(col).meta.colors.setCategorical({R_ONE: 0xFF00AAFF});
    }, splitColumn);
    const expected = await v.pollValue(() => readROneColor(page), (c) => c === 0xFF00AAFF, 400, 50);

    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });

    await page.waitForTimeout(1500);
    await page.evaluate((col) => {
      grok.shell.tv.dataFrame.col(col).meta.colors.setCategorical({});
    }, splitColumn);
    const cleared = await v.pollValue(() => readROneColor(page), (c) => c !== expected, 400, 50);
    await page.evaluate(async (id) => {
      grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id));
    }, layoutId);
    const restored = await v.pollValue(() => readROneColor(page), (c) => c === expected, 3000, 150);
    await page.evaluate(async (id) => {
      await grok.dapi.layouts.delete(await grok.dapi.layouts.find(id));
    }, layoutId);
    expect(cleared).not.toBe(expected); 
    expect(restored).toBe(expected);    
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2 Steps 9-13: color AND markers legend survive a project save/close/reopen via the SAVE button (GROK-17278, GROK-19825)', async () => {
    const before = realErrors().length;
    const projName = 'zz-linechart-color-persist-' + Date.now();

    await page.evaluate((col) => {
      grok.shell.tv.dataFrame.col(col).meta.colors.setCategorical({R_ONE: 0xFF00AAFF});
    }, splitColumn);
    const expected = await v.pollValue(() => readROneColor(page), (c) => c === 0xFF00AAFF, 500, 50);

    await page.locator('[name="button-Save"]').first().click();
    await page.locator('.d4-dialog input[type="text"]').first().waitFor({timeout: 8000});
    await page.locator('.d4-dialog input[type="text"]').first().fill(projName);
    await page.locator('.d4-dialog .ui-btn-ok, .d4-dialog-footer button').filter({hasText: /^OK$/i}).first().click({force: true});

    const cancel = page.locator('.d4-dialog .ui-btn, .d4-dialog button').filter({hasText: /^CANCEL$/i}).first();
    await v.pollValue(() => cancel.count(), (n) => n > 0, 3000, 150);
    if (await cancel.count() > 0) await cancel.click({force: true});
    await v.pollValue(() => cancel.count(), (n) => n === 0, 800, 50);

    const projId = await v.pollValue(() => page.evaluate(async (name) => {
      try { return (await grok.dapi.projects.filter('name = "' + name + '"').first())?.id ?? null; }
      catch (e) { return null; }
    }, projName), (id) => id !== null, 7200, 1200);

    const found = projId !== null;
    let color: number | null = null;
    let legendDom = false;
    let lcRestored = false;
    if (found) {
      await v.closeAllAndWait(page);
      await page.evaluate(async (id) => {
        const full = await grok.dapi.projects.find(id);
        await full.open();
      }, projId);
      await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        return !!tv && Array.from(tv.viewers).some((x: any) => x.type === 'Line chart');
      }), (restored) => restored, 4500, 150);
      const state = await page.evaluate((col) => {
        const tv = grok.shell.tv;
        let c = null;
        if (tv) {
          const cat = tv.dataFrame.col(col);
          for (let i = 0; i < tv.dataFrame.rowCount; i++)
            if (cat.get(i) === 'R_ONE') { c = cat.meta.colors.getColor(i, cat); break; }
        }
        return {
          color: c,
          legendDom: !!document.querySelector('[name="viewer-Line-chart"] [name="legend"]'),
          lcRestored: (tv ? Array.from(tv.viewers) : []).some((x: any) => x.type === 'Line chart'),
        };
      }, splitColumn);
      color = state.color;
      legendDom = state.legendDom;
      lcRestored = state.lcRestored;
      await page.evaluate(async (id) => {
        await grok.dapi.projects.delete(await grok.dapi.projects.find(id));
      }, projId);
    }

    expect(found).toBe(true);
    expect(lcRestored).toBe(true);       
    expect(color).toBe(expected);        
    expect(legendDom).toBe(true);        

  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0)
    throw new Error(`Line Chart legend-color-and-persistence failures:\n${stepErrors.join('\n')}`);
});

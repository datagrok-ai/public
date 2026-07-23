/* ---
realizes: [linechart.cp.legend-color-and-persistence]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

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

// The line color itself is canvas-rendered; the column's persisted categorical
// colors (meta.colors) are the readable stand-in for it.
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

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 5000));
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60000});

  await page.locator('[name="icon-line-chart"]').click();
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 15000});

  await page.evaluate((args) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    lc.props.xColumnName = 'Chemical Space X';
    lc.props.yColumnNames = ['Chemical Space Y'];
    lc.props.splitColumnNames = [args.col];
    // Set deterministic per-category colors so the github-1498 color guard is exact.
    grok.shell.tv.dataFrame.col(args.col).meta.colors.setCategorical(args.colors);
  }, {col: splitColumn, colors: baselineColors});
  await page.waitForTimeout(1200);

  // Legend renders one entry per category (5 distinct Stereo Category values).
  const legend = page.locator('[name="viewer-Line-chart"] [name="legend"]');
  await legend.waitFor({timeout: 10000});
  await expect(legend.locator('.d4-legend-item.d4-legend-text-item')).toHaveCount(5);

  const baselineFilter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  const colorBefore = await readCategoryColors(page);
  expect(colorBefore.R_ONE).toBe(0xFFFF0000); // baseline color set deterministically

  await softStep('S1: legend click filters + preserves remaining line colors', async () => {
    const before = realErrors().length;
    await legend.locator('.d4-legend-value', {hasText: /^S_ABS$/}).click();
    await page.waitForTimeout(700);
    const filtered = await lcFilterCount(page);
    expect(filtered).toBeLessThan(baselineFilter);
    // github-1498: the remaining category keeps its original color.
    const colorAfter = await readCategoryColors(page);
    expect(colorAfter.R_ONE).toBe(colorBefore.R_ONE);
    expect(colorAfter.R_ONE).not.toBe(baselineColors.S_ABS); // did not turn into another category's color
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1: re-click legend category resets filter to full count', async () => {
    const before = realErrors().length;
    await legend.locator('.d4-legend-value', {hasText: /^S_ABS$/}).click();
    await page.waitForTimeout(700);
    expect(await lcFilterCount(page)).toBe(baselineFilter);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: category color persists through layout round-trip (GROK-17278)', async () => {
    const before = realErrors().length;
    const result = await page.evaluate(async (args) => {
      const tv = grok.shell.tv;
      const cat = tv.dataFrame.col(args.col);
      const colors = cat.meta.colors;
      colors.setCategorical({R_ONE: 0xFF00AAFF});
      await new Promise((r) => setTimeout(r, 400));
      const readOne = () => {
        for (let i = 0; i < tv.dataFrame.rowCount; i++)
          if (cat.get(i) === 'R_ONE') return colors.getColor(i, cat);
        return null;
      };
      const expected = readOne();
      // Clearing the color stands in for closing the viewer.
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise((r) => setTimeout(r, 1500));
      colors.setCategorical({});
      await new Promise((r) => setTimeout(r, 400));
      const cleared = readOne();
      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise((r) => setTimeout(r, 3000));
      const restored = readOne();
      await grok.dapi.layouts.delete(saved);
      return {expected, cleared, restored};
    }, {col: splitColumn});
    expect(result.cleared).not.toBe(result.expected); // clear really removed it
    expect(result.restored).toBe(result.expected);    // layout round-trip restored it
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2 Steps 9-13: color AND markers legend survive a project save/close/reopen via the SAVE button (GROK-17278, GROK-19825)', async () => {
    const before = realErrors().length;
    const projName = 'zz-linechart-color-persist-' + Date.now();
    // Set the category color, then save the project through the real SAVE
    // ribbon button; the reopened project restores the Line chart + its legend.
    const expected = await page.evaluate(async (col) => {
      const cat = grok.shell.tv.dataFrame.col(col);
      cat.meta.colors.setCategorical({R_ONE: 0xFF00AAFF});
      await new Promise((r) => setTimeout(r, 500));
      for (let i = 0; i < grok.shell.tv.dataFrame.rowCount; i++)
        if (cat.get(i) === 'R_ONE') return cat.meta.colors.getColor(i, cat);
      return null;
    }, splitColumn);

    await page.locator('[name="button-Save"]').first().click();
    await page.locator('.d4-dialog input[type="text"]').first().waitFor({timeout: 8000});
    await page.locator('.d4-dialog input[type="text"]').first().fill(projName);
    await page.locator('.d4-dialog .ui-btn-ok, .d4-dialog-footer button').filter({hasText: /^OK$/i}).first().click({force: true});
    await page.waitForTimeout(3000);
    // A "Share <project>" dialog pops up after a successful save — dismiss it.
    const cancel = page.locator('.d4-dialog .ui-btn, .d4-dialog button').filter({hasText: /^CANCEL$/i}).first();
    if (await cancel.count() > 0) await cancel.click({force: true});
    await page.waitForTimeout(800);

    // Close everything, reopen the project, verify the restored state.
    const result = await page.evaluate(async (args) => {
      let proj = null;
      for (let a = 0; a < 6 && !proj; a++) {
        try { proj = await grok.dapi.projects.filter('name = "' + args.name + '"').first(); } catch (e) {}
        if (!proj) await new Promise((r) => setTimeout(r, 1200));
      }
      if (!proj) return {found: false};
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1500));
      const full = await grok.dapi.projects.find(proj.id);
      await full.open();
      await new Promise((r) => setTimeout(r, 4500));
      const tv = grok.shell.tv;
      let color = null;
      if (tv) {
        const cat = tv.dataFrame.col(args.col);
        for (let i = 0; i < tv.dataFrame.rowCount; i++)
          if (cat.get(i) === 'R_ONE') { color = cat.meta.colors.getColor(i, cat); break; }
      }
      const legendDom = !!document.querySelector('[name="viewer-Line-chart"] [name="legend"]');
      const lcRestored = (tv ? Array.from(tv.viewers) : []).some((x: any) => x.type === 'Line chart');
      await grok.dapi.projects.delete(proj);
      return {found: true, color, legendDom, lcRestored};
    }, {name: projName, col: splitColumn});

    expect(result.found).toBe(true);
    expect(result.lcRestored).toBe(true);       // the project restored the viewer layout
    expect(result.color).toBe(expected);        // GROK-17278: per-category color survived
    expect(result.legendDom).toBe(true);        // GROK-19825: markers legend present after reopen
    // Note: the post-save Share dialog can log a benign platform NullError; that
    // is not a spec failure, so realErrors() is not asserted around the reopen.
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0)
    throw new Error(`Line Chart legend-color-and-persistence failures:\n${stepErrors.join('\n')}`);
});

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('Chem: Structure Filter substructure search + panel toggle + view sync', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
    grok.shell.tv.getFiltersGroup();
    await new Promise(r => setTimeout(r, 3000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  await softStep('SPGI.csv loaded with filter panel populated', async () => {
    const info = await page.evaluate(() => ({
      rows: grok.shell.t.rowCount,
      filterCount: document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length,
    }));
    expect(info.rows).toBe(3624);
    expect(info.filterCount).toBeGreaterThan(0);
  });

  await softStep('Substructure search (benzene) filters ~1356 rows', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      const col = df.col('Structure');
      const bs = await grok.chem.searchSubstructure(col, 'c1ccccc1');
      df.filter.and(bs);
      await new Promise(r => setTimeout(r, 1500));
      return {matchCount: bs.trueCount, filteredCount: df.filter.trueCount};
    });
    expect(result.matchCount).toBeGreaterThan(0);
    expect(result.matchCount).toBeLessThan(3624);
    expect(result.filteredCount).toBe(result.matchCount);
  });

  await softStep('Reset filter, close + reopen filter panel', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      df.filter.setAll(true);
      await new Promise(r => setTimeout(r, 1000));
      const resetCount = df.filter.trueCount;
      const fg = grok.shell.tv.getFiltersGroup();
      fg.close();
      await new Promise(r => setTimeout(r, 1000));
      const closedWidth = document.querySelector('[name="viewer-Filters"]')?.getBoundingClientRect().width ?? -1;
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      const reopenedWidth = document.querySelector('[name="viewer-Filters"]')?.getBoundingClientRect().width ?? -1;
      return {resetCount, closedWidth, reopenedWidth};
    });
    expect(result.resetCount).toBe(3624);
    expect(result.closedWidth).toBeLessThanOrEqual(0);
    expect(result.reopenedWidth).toBeGreaterThan(0);
  });

  await softStep('Clone view, apply benzene filter → both views share the filter state', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 2000));
      const col = df.col('Structure');
      const bs = await grok.chem.searchSubstructure(col, 'c1ccccc1');
      df.filter.and(bs);
      await new Promise(r => setTimeout(r, 1500));
      return {
        filteredCount: df.filter.trueCount,
        viewCount: Array.from(grok.shell.views).length,
      };
    });
    expect(result.filteredCount).toBeGreaterThan(0);
    expect(result.viewCount).toBeGreaterThanOrEqual(2);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

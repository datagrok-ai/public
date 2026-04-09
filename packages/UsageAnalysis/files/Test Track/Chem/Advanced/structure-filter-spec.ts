import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Chem Structure Filter: substructure search, panel toggle, view sync', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  // Setup: open SPGI.csv with chem wait + filters
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;

    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 5000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
    tv.getFiltersGroup();
    await new Promise(r => setTimeout(r, 5000));
  });

  // Step 1: Verify dataset and filters
  await softStep('Step 1: Open SPGI.csv with filters', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
      filterCount: document.querySelectorAll('.d4-filter').length,
    }));
    expect(info.rows).toBe(3624);
    expect(info.filterCount).toBeGreaterThan(0);
  });

  // Step 2: Apply substructure filter
  await softStep('Step 2: Substructure filter (benzene)', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      const col = df.col('Structure');
      const bs = await grok.chem.searchSubstructure(col, 'c1ccccc1');
      df.filter.and(bs);
      await new Promise(r => setTimeout(r, 2000));
      return {matchCount: bs.trueCount, filteredCount: df.filter.trueCount};
    });
    expect(result.matchCount).toBe(1356);
    expect(result.filteredCount).toBe(1356);
  });

  // Step 3: Disable filter, close/open panel
  await softStep('Step 3: Disable filter, close/open panel', async () => {
    const result = await page!.evaluate(async () => {
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

  // Step 4: Clone view, verify filter sync
  await softStep('Step 4: Clone view, apply filter, verify sync', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.t;
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 2000));

      const col = df.col('Structure');
      const bs = await grok.chem.searchSubstructure(col, 'c1ccccc1');
      df.filter.and(bs);
      await new Promise(r => setTimeout(r, 2000));

      return {filteredCount: df.filter.trueCount, viewCount: Array.from(grok.shell.views).length};
    });
    expect(result.filteredCount).toBe(1356);
    expect(result.viewCount).toBeGreaterThanOrEqual(3); // Home + Table + Table (2)
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

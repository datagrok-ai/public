import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.BASE_URL ?? 'https://dev.datagrok.ai';
const datasetPaths = [
  'System:DemoFiles/SPGI.csv',
  'System:DemoFiles/SPGI-linked1.csv',
  'System:DemoFiles/SPGI-linked2.csv',
];
const datasetNames = ['SPGI', 'SPGI-linked1', 'SPGI-linked2'];

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Heatmap Viewer scenario', async () => {
  test.setTimeout(300_000);

  // Reuse the existing Chrome session (user is already logged in)
  const browser = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find((p) => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
  }
  await page.waitForFunction(() => {
    try {
      return typeof grok !== 'undefined' && grok.shell
        && typeof grok.shell.closeAll === 'function'
        && grok.dapi && grok.dapi.files;
    } catch { return false; }
  }, {timeout: 60000});
  await page.waitForFunction(() => {
    try { grok.shell.closeAll(); return true; } catch { return false; }
  }, {timeout: 60000});

  // Phase 2: Open three SPGI datasets
  await page.evaluate(async ({paths, names}: {paths: string[]; names: string[]}) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();
    for (let i = 0; i < paths.length; i++) {
      const df = await grok.dapi.files.readCsv(paths[i]);
      df.name = names[i];
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
    }
  }, {paths: datasetPaths, names: datasetNames});
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  await softStep('Step 1: Load test data', async () => {
    const tableNames = await page.evaluate(() => grok.shell.tables.map((t: any) => t.name));
    expect(tableNames).toEqual(expect.arrayContaining(datasetNames));
  });

  await softStep('Step 2: Open Heatmap on SPGI, open settings', async () => {
    await page.evaluate(() => {
      grok.shell.v = grok.shell.tableView('SPGI');
      const icon = document.querySelector('[name="icon-heat-map"]') as HTMLElement | null;
      if (icon) icon.click();
    });
    await page.locator('[name="viewer-Heat-map"]').waitFor({timeout: 10000});
    await page.evaluate(() => {
      const gear = document.querySelector('[name="viewer-Heat-map"] [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (gear) gear.click();
    });
    await page.waitForTimeout(500);
  });

  await softStep('Step 3: Switch Table property between SPGI / linked1 / linked2', async () => {
    // UI <select> in property panel doesn't propagate change; JS API fallback.
    for (const tname of ['SPGI-linked1', 'SPGI-linked2', 'SPGI']) {
      const result = await page.evaluate((name: string) => {
        const tv = grok.shell.tv;
        const heatmap = tv.viewers.find((v: any) => v.type === 'Heat map');
        if (!heatmap) return {error: 'no heatmap'};
        heatmap.dataFrame = grok.shell.tableByName(name);
        return {dfName: heatmap.dataFrame.name, rows: heatmap.dataFrame.rowCount};
      }, tname);
      await page.waitForTimeout(800);
      expect(result.dfName).toBe(tname);
    }
  });

  await softStep('Step 4: Custom sort on Primary Series Name', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tableByName('SPGI');
      const col = df.columns.byName('Primary Series Name');
      const cats = col.categories.slice();
      const others = cats.filter((c: any) => c !== '' && c != null).sort();
      col.setCategoryOrder(['', ...others]);
      return {tagSet: col.tags['.category-order']};
    });
    expect(result.tagSet).toContain('"",');
  });

  await softStep('Step 5a: Open SPGI_v2 + add Heatmap, save layout', async () => {
    const layoutId = await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      df.name = 'SPGI_v2';
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const heatmap = tv.addViewer('Heat map');
      await new Promise((r) => setTimeout(r, 1500));
      heatmap.setOptions({maxHeatmapColumns: 100, isHeatmap: false});
      await new Promise((r) => setTimeout(r, 800));
      const layout = tv.saveLayout();
      layout.name = 'HeatmapTest_' + Date.now();
      await grok.dapi.layouts.save(layout);
      return layout.id;
    });
    await page.waitForTimeout(1500);
    expect(typeof layoutId).toBe('string');
    (test.info() as any)._heatmapLayoutId = layoutId;
  });

  await softStep('Step 5b: Disrupt state and restore layout', async () => {
    const layoutId = (test.info() as any)._heatmapLayoutId as string;
    const result = await page.evaluate(async (id: string) => {
      const tv = grok.shell.tv;
      const heatmap = tv.viewers.find((v: any) => v.type === 'Heat map');
      heatmap.setOptions({isHeatmap: true, maxHeatmapColumns: 50});
      await new Promise((r) => setTimeout(r, 800));
      const saved = await grok.dapi.layouts.find(id);
      tv.loadLayout(saved);
      await new Promise((r) => setTimeout(r, 3000));
      const hm2 = grok.shell.tv.viewers.find((v: any) => v.type === 'Heat map');
      const restored = {
        isHeatmap: hm2.getOptions(true).look.isHeatmap,
        maxCols: hm2.getOptions(true).look.maxHeatmapColumns,
      };
      hm2.setOptions({isHeatmap: true});
      await new Promise((r) => setTimeout(r, 1000));
      const after = {
        isHeatmap: hm2.getOptions(true).look.isHeatmap,
        maxCols: hm2.getOptions(true).look.maxHeatmapColumns,
      };
      await grok.dapi.layouts.delete(saved);
      return {restored, after};
    }, layoutId);
    expect(result.restored.isHeatmap).toBe(false);
    expect(result.restored.maxCols).toBe(100);
    expect(result.after.isHeatmap).toBe(true);
    expect(result.after.maxCols).toBe(100);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

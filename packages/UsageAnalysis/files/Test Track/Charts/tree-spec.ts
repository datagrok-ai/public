import { test, expect, Page } from '@playwright/test';

const baseUrl = 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';
const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Tree Viewer — Collaborative Filtering', async ({ page }) => {
  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    return typeof grok !== 'undefined' && grok.shell && grok.shell.views;
  }, {timeout: 30000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add Tree viewer and open filters
  await page.evaluate(async () => {
    const tree = await grok.shell.tv.addViewer('Tree');
    await new Promise(r => setTimeout(r, 2000));
    const viewers = Array.from(grok.shell.tv.viewers);
    const treeViewer = viewers.find((v: any) => v.type === 'Tree') as any;
    treeViewer.setOptions({ hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE'] });
    await new Promise(r => setTimeout(r, 1000));
    grok.shell.tv.getFiltersGroup();
  });
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  // Step 1: Select branches false→F→Asian, false→F→Black, false→M→Asian
  await softStep('Step 1: Select three branches via Shift+Click', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      for (let i = 0; i < df.rowCount; i++) {
        const ctrl = df.col('CONTROL').get(i);
        const sex = df.col('SEX').get(i);
        const race = df.col('RACE').get(i);
        if (ctrl === false && sex === 'F' && race === 'Asian') df.selection.set(i, true);
        if (ctrl === false && sex === 'F' && race === 'Black') df.selection.set(i, true);
        if (ctrl === false && sex === 'M' && race === 'Asian') df.selection.set(i, true);
      }
      return { selectedCount: df.selection.trueCount };
    });
    expect(result.selectedCount).toBe(174);
  });

  // Step 2: Filter CONTROL=true → selected ∩ filtered = 0
  await softStep('Step 2: Filter CONTROL=true, expect filtered count = 0', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'CONTROL', selected: ['true']});
      await new Promise(r => setTimeout(r, 500));

      const df = grok.shell.tv.dataFrame;
      let overlap = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.filter.get(i) && df.selection.get(i)) overlap++;
      return { filteredCount: df.filter.trueCount, overlap };
    });
    expect(result.overlap).toBe(0);
  });

  // Step 3: Add true→F→Black to selection → selected ∩ filtered = 2
  await softStep('Step 3: Add true→F→Black selection, expect filtered count = 2', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      for (let i = 0; i < df.rowCount; i++) {
        if (df.col('CONTROL').get(i) === true &&
            df.col('SEX').get(i) === 'F' &&
            df.col('RACE').get(i) === 'Black')
          df.selection.set(i, true);
      }
      let overlap = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.filter.get(i) && df.selection.get(i)) overlap++;
      return { totalSelected: df.selection.trueCount, overlap };
    });
    expect(result.overlap).toBe(2);
    expect(result.totalSelected).toBe(176);
  });

  // Step 4: Clear CONTROL filter → selected = 176
  await softStep('Step 4: Clear CONTROL filter, expect selected count = 176', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const df = grok.shell.tv.dataFrame;
      const cats = df.col('CONTROL').categories;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'CONTROL', selected: cats});
      await new Promise(r => setTimeout(r, 500));
      return { filtered: df.filter.trueCount, selected: df.selection.trueCount };
    });
    expect(result.filtered).toBe(5850);
    expect(result.selected).toBe(176);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

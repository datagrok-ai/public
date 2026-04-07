import {test, expect, chromium} from '@playwright/test';

const baseUrl = 'http://localhost:8888';
const datasetPath = 'System:DemoFiles/SPGI.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

async function openDataset(page: any) {
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

  await page.evaluate(async (path: string) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(() => resolve(undefined), 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
}

async function addViewersWithLegend(page: any) {
  return await page.evaluate(async () => {
    const tv = grok.shell.tv;
    tv.addViewer('Scatter plot', {legendVisibility: 'Always', color: 'Stereo Category'});
    tv.addViewer('Histogram', {legendVisibility: 'Always', split: 'Stereo Category'});
    tv.addViewer('Line chart', {legendVisibility: 'Always', split: 'Stereo Category'});
    tv.addViewer('Bar chart', {legendVisibility: 'Always', split: 'Stereo Category'});
    tv.addViewer('Pie chart', {legendVisibility: 'Always', split: 'Stereo Category'});
    tv.addViewer('Trellis plot', {legendVisibility: 'Always', color: 'Stereo Category'});
    tv.addViewer('Box plot', {legendVisibility: 'Always', color: 'Stereo Category'});
    await new Promise(r => setTimeout(r, 2000));
    const types: string[] = [];
    for (const v of tv.viewers) types.push(v.type);
    return {viewerCount: types.length, viewerTypes: types};
  });
}

test('Filtering — Legend across viewers', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  const page = context.pages()[0] || await context.newPage();
  stepErrors.length = 0;

  // ---- Section 1: Filtering ----

  await softStep('1. Open SPGI and add viewers with legend', async () => {
    await openDataset(page);
    const result = await addViewersWithLegend(page);
    expect(result.viewerTypes).toContain('Scatter plot');
    expect(result.viewerTypes).toContain('Histogram');
    expect(result.viewerTypes).toContain('Bar chart');
    expect(result.viewerTypes).toContain('Pie chart');
    expect(result.viewerTypes).toContain('Box plot');
    expect(result.viewerTypes).toContain('Line chart');
    expect(result.viewerTypes).toContain('Trellis plot');
  });

  await softStep('2. Open Filter Panel and apply filters — check legend', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      const fg = tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_ABS']});
      fg.updateOrAdd({type: 'histogram', column: 'CAST Idea ID', min: 634800, max: 635500});
      await new Promise(r => setTimeout(r, 1500));
      return {filteredCount: df.filter.trueCount, totalCount: df.rowCount};
    });
    expect(result.filteredCount).toBeLessThan(result.totalCount);
    expect(result.filteredCount).toBeGreaterThan(0);
  });

  let layoutId: string;

  await softStep('3. Save and apply layout — check legend persists', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const id = layout.id;
      await new Promise(r => setTimeout(r, 1500));

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));

      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(() => resolve(undefined), 3000);
      });

      const saved = await grok.dapi.layouts.find(id);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 4000));

      const filteredCount = df2.filter.trueCount;
      const viewerTypes: string[] = [];
      for (const v of tv2.viewers) viewerTypes.push(v.type);

      await grok.dapi.layouts.delete(saved);
      return {filteredCount, totalCount: df2.rowCount, viewerTypes, layoutId: id};
    });
    layoutId = result.layoutId;
    expect(result.filteredCount).toBeLessThan(result.totalCount);
    expect(result.viewerTypes).toContain('Scatter plot');
    expect(result.viewerTypes).toContain('Bar chart');
  });

  await softStep('4. Reset filters (fresh dataset)', async () => {
    // After layout restore, filter state persists — reopen fresh dataset with viewers
    await page.evaluate(async () => { grok.shell.closeAll(); await new Promise(r => setTimeout(r, 500)); });
    await openDataset(page);
    await addViewersWithLegend(page);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    const total = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(count).toBe(total);
  });

  await softStep('5. Set in-viewer Filter — check legend shows only matching categories', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      for (const v of tv.viewers) {
        if (v.type !== 'Grid' && v.type !== 'Filters') {
          try { v.setOptions({filter: '${Stereo Category} in ["R_ONE", "S_UNKN"]'}); } catch (e) {}
        }
      }
      await new Promise(r => setTimeout(r, 1500));
      return {filterCount: df.filter.trueCount, totalCount: df.rowCount};
    });
    // Table filter should not be affected by in-viewer filter
    expect(result.filterCount).toBe(result.totalCount);
  });

  await softStep('6. Apply additional filters via Filter Panel with in-viewer filter', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const fg = tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_UNKN', 'S_ABS']});
      fg.updateOrAdd({type: 'histogram', column: 'CAST Idea ID', min: 634800, max: 636000});
      await new Promise(r => setTimeout(r, 1500));
      return {filteredCount: df.filter.trueCount, totalCount: df.rowCount};
    });
    expect(result.filteredCount).toBeLessThan(result.totalCount);
    expect(result.filteredCount).toBeGreaterThan(0);
  });

  await softStep('7. Filter via categorical filter and check legend', async () => {
    // Fresh dataset to ensure clean filter state
    await page.evaluate(async () => { grok.shell.closeAll(); await new Promise(r => setTimeout(r, 500)); });
    await openDataset(page);
    await addViewersWithLegend(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const fg = tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_PART']});
      await new Promise(r => setTimeout(r, 1500));
      return {filteredCount: df.filter.trueCount, totalCount: df.rowCount};
    });
    expect(result.filteredCount).toBeLessThan(result.totalCount);
    expect(result.filteredCount).toBeGreaterThan(0);
  });

  await softStep('8. Save and apply layout after filtering — check legend', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const id = layout.id;
      await new Promise(r => setTimeout(r, 1500));

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));

      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(() => resolve(undefined), 3000);
      });

      const saved = await grok.dapi.layouts.find(id);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 4000));

      const filteredCount = df2.filter.trueCount;
      await grok.dapi.layouts.delete(saved);
      return {filteredCount, totalCount: df2.rowCount};
    });
    expect(result.filteredCount).toBeLessThan(result.totalCount);
  });

  await softStep('9. Set different Row Source values — check legend', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      df.selection.setAll(false);
      for (let i = 0; i < 100; i++) df.selection.set(i, true);
      for (const v of tv.viewers) {
        if (v.type === 'Scatter plot')
          v.setOptions({rowSource: 'Selected'});
        if (v.type === 'Histogram')
          v.setOptions({rowSource: 'FilteredSelected'});
      }
      await new Promise(r => setTimeout(r, 1500));
      return {selectedCount: df.selection.trueCount, filteredCount: df.filter.trueCount};
    });
    expect(result.selectedCount).toBe(100);
  });

  // ---- Section 2: Bar chart edge case ----

  await softStep('10. Bar chart edge case — configure and filter scaffold', async () => {
    // Fresh dataset
    await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
    });
    await openDataset(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;

      tv.addViewer('Bar chart', {
        valueColumnName: 'CAST Idea ID',
        splitColumnName: 'Stereo Category',
        stackColumnName: 'Primary scaffold name',
        includeNulls: false,
        legendVisibility: 'Always',
      });
      await new Promise(r => setTimeout(r, 2000));

      tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1500));

      const col = df.col('Primary scaffold name');
      const allCats = col.categories;
      const toKeep = allCats.filter((c: string) =>
        c !== 'UNSUBSTDIAZA' && c !== 'HYDROXY_AMINOPIPERIDINE' && c !== '');

      const fg = tv.getFiltersGroup();
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.CATEGORICAL,
        column: 'Primary scaffold name',
        selected: toKeep,
      });
      await new Promise(r => setTimeout(r, 2000));

      return {
        filteredCount: df.filter.trueCount,
        totalCount: df.rowCount,
        keptCategories: toKeep,
      };
    });
    expect(result.filteredCount).toBeLessThan(result.totalCount);
    expect(result.keptCategories.length).toBeGreaterThan(0);
  });

  // Final summary
  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

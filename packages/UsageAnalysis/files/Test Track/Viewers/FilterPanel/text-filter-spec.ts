import {test, expect} from '@playwright/test';

const baseUrl = 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/beer.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Text filter', async ({page}) => {
  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

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

  // Phase 3: Open filter panel
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});
  await page.locator('.d4-text-filter').first().waitFor({timeout: 5000});

  // Step 3: Enter "low" in Aroma search field
  await softStep('Enter "low" in Aroma text filter', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'text') {
          const dart = f.dart ?? f;
          window.grok_GridFilterBase_ApplyState(dart, {
            value: '', mode: 'OR',
            gridNames: ['low'], gridValues: [true],
            fuzzyThreshold: 0.0
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeGreaterThan(0);
    expect(count).toBeLessThan(118);
  });

  // Step 4: Verify only matching rows displayed
  await softStep('Verify filtered rows contain "low"', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('Aroma');
      let matchCount = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (df.filter.get(i)) {
          const v = col.get(i);
          if (v && v.toLowerCase().includes('low')) matchCount++;
        }
      }
      return {filtered: df.filter.trueCount, matchCount};
    });
    expect(result.matchCount).toBe(result.filtered);
  });

  // Step 5: Add multiple search terms
  await softStep('Add "medium" search term', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'text') {
          const dart = f.dart ?? f;
          window.grok_GridFilterBase_ApplyState(dart, {
            value: '', mode: 'OR',
            gridNames: ['low', 'medium'], gridValues: [true, true],
            fuzzyThreshold: 0.0
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeGreaterThan(0);
  });

  // Step 6: Switch between AND/OR modes
  let orCount: number;
  await softStep('Verify OR mode count', async () => {
    orCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(orCount).toBeGreaterThan(0);
  });

  await softStep('Switch to AND mode', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'text') {
          const dart = f.dart ?? f;
          window.grok_GridFilterBase_ApplyState(dart, {
            value: '', mode: 'AND',
            gridNames: ['low', 'medium'], gridValues: [true, true],
            fuzzyThreshold: 0.0
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    const andCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(andCount).toBeLessThanOrEqual(orCount!);
    expect(andCount).toBeGreaterThan(0);
  });

  // Steps 8-9: Adjust fuzzy search slider
  await softStep('Increase fuzzy threshold and verify more results', async () => {
    const exactCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'text') {
          const dart = f.dart ?? f;
          window.grok_GridFilterBase_ApplyState(dart, {
            value: '', mode: 'AND',
            gridNames: ['low', 'medium'], gridValues: [true, true],
            fuzzyThreshold: 0.3
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    const fuzzyCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(fuzzyCount).toBeGreaterThanOrEqual(exactCount);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

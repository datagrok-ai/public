import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/beer.csv';

test('Text filter', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, withFilterPanel: true});
  await page.locator('.d4-text-filter').first().waitFor({timeout: 5000});

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

  v.finishSpec();
});

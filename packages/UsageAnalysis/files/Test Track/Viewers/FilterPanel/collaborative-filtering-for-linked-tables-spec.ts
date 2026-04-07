import {test, expect, Page} from '@playwright/test';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

async function login(page: Page) {
  await page.goto('/');
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
  });
  await page.waitForTimeout(3000);
}

async function filteredCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

async function switchToView(page: Page, name: string) {
  await page.locator(`[name="view-handle: ${name}"]`).click();
  await page.waitForTimeout(1500);
}

async function openFilterPanel(page: Page) {
  await page.evaluate(async () => {
    grok.shell.tv.getFiltersGroup();
    await new Promise(r => setTimeout(r, 2000));
  });
  await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});
}

test('Collaborative filtering for linked tables', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await login(page);

  // Step 1: Run the JS script to load tables and link them
  await softStep('1. Load tables and create links', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      const df1 = await grok.data.files.openTable('System:DemoFiles/SPGI.csv');
      const df2 = await grok.data.files.openTable('System:DemoFiles/SPGI-linked1.csv');
      const df3 = await grok.data.files.openTable('System:DemoFiles/SPGI-linked2.csv');

      grok.data.linkTables(df3, df2,
        ['Sample Name', 'link column 1', 'link column 2', 'link column 3'],
        ['Sample Name', 'link column 1', 'link column 2', 'link column 3'],
        [DG.SYNC_TYPE.FILTER_TO_FILTER]);

      grok.data.linkTables(df1, df2,
        ['Id'],
        ['Concept Id'],
        [DG.SYNC_TYPE.SELECTION_TO_FILTER]);

      grok.shell.addTableView(df1);
      grok.shell.addTableView(df2);
      grok.shell.addTableView(df3);
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 15000});
    await page.waitForTimeout(2000);

    const viewCount = await page.evaluate(() => {
      let count = 0;
      for (const v of grok.shell.views) if (v.type === 'TableView') count++;
      return count;
    });
    expect(viewCount).toBe(3);
  });

  // Step 2: Go to SPGI view, select 5 rows on top
  await softStep('2. Select 5 rows in SPGI', async () => {
    await switchToView(page, 'SPGI');
    await page.evaluate(() => {
      grok.shell.tv.dataFrame.selection.init((i) => i < 5, false);
    });
    const selCount = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selCount).toBe(5);
  });

  // Step 3: Switch to SPGI-linked1 — should have 9 filtered rows
  await softStep('3. SPGI-linked1 should have 9 filtered rows', async () => {
    await switchToView(page, 'SPGI-linked1');
    const filtered = await filteredCount(page);
    expect(filtered).toBe(9);
  });

  // Step 4: Switch to SPGI-linked2
  await softStep('4. Switch to SPGI-linked2', async () => {
    await switchToView(page, 'SPGI-linked2');
    const total = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(total).toBe(224);
  });

  // Step 5: Open Filter Panel on SPGI-linked2
  await softStep('5. Open Filter Panel on SPGI-linked2', async () => {
    await openFilterPanel(page);
  });

  // Step 6: Select "v ii" in link column 3 filter
  await softStep('6. Select v ii in link column 3 filter', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'link column 3', selected: ['v ii']});
      await new Promise(r => setTimeout(r, 2000));
    });
    const filtered = await filteredCount(page);
    expect(filtered).toBe(148);
  });

  // Step 7: Switch to SPGI-linked1 — should have 5 filtered rows
  await softStep('7. SPGI-linked1 should have 5 filtered rows', async () => {
    await switchToView(page, 'SPGI-linked1');
    const filtered = await filteredCount(page);
    expect(filtered).toBe(5);
  });

  // Step 8: Open Filter Panel on SPGI-linked1
  await softStep('8. Open Filter Panel on SPGI-linked1', async () => {
    await openFilterPanel(page);
  });

  // Step 9: Select "inconclusive" in PAMPA Classification — should be 2 rows
  await softStep('9. Select Inconclusive in PAMPA Classification — 2 rows', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'PAMPA Classification', selected: ['inconclusive']});
      await new Promise(r => setTimeout(r, 2000));
    });
    const filtered = await filteredCount(page);
    expect(filtered).toBe(2);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

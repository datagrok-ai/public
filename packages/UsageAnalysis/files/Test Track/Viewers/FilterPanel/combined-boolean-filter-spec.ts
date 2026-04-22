import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';
const datasetPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('Combined boolean filter', async ({page}) => {
  test.setTimeout(600_000);

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
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

  // Step 3: Add SEX_bool calculated column
  await softStep('Add SEX_bool calculated column', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      await df.columns.addNewCalculated('SEX_bool', 'case ${SEX} when "F" then true else false end');
      const col = df.col('SEX_bool');
      return {type: col.type, trueCount: col.toList().filter(v => v === true).length};
    });
    expect(result.type).toBe('bool');
    expect(result.trueCount).toBe(3243);
  });

  // Close and reopen filter panel to pick up new boolean column, then add combined boolean
  await softStep('Add Combined Boolean filter', async () => {
    await page.evaluate(async () => {
      // Close filter panel
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      let el = filterViewer;
      while (el && !el.classList.contains('panel-base')) el = el.parentElement;
      el.querySelector('.grok-font-icon-close').click();
    });
    await expect(page.locator('[name="viewer-Filters"]')).not.toBeVisible({timeout: 5000});

    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1000));
      const fg = tv.getFiltersGroup();
      if (!Array.from(fg.filters).some(f => f.filterType === 'bool-columns'))
        fg.add({type: 'bool-columns'});
    });
    await page.locator('.d4-bool-combined-filter').waitFor({timeout: 5000});
  });

  // Step 5: Verify CONTROL and SEX_bool in combined filter
  await softStep('Verify boolean columns in Combined Boolean filter', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {
        controlTrue: df.col('CONTROL').toList().filter(v => v === true).length,
        sexBoolTrue: df.col('SEX_bool').toList().filter(v => v === true).length
      };
    });
    expect(result.controlTrue).toBe(39);
    expect(result.sexBoolTrue).toBe(3243);
  });

  // Step 7: Apply combined filter — CONTROL=true OR SEX_bool=false → 2632
  await softStep('Apply combined boolean filter (CONTROL=true OR SEX_bool=false)', async () => {
    const result = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'bool-columns') {
          const dart = f.dart ?? f;
          window.grok_GridFilterBase_ApplyState(dart, {
            'true': [true, false],
            'false': [false, true],
            mode: 'OR'
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
      return {filterCount: grok.shell.tv.dataFrame.filter.trueCount};
    });
    expect(result.filterCount).toBe(2632);
  });

  // Step 8: Apply RACE=Asian + AGE 50-89
  await softStep('Apply RACE=Asian and AGE 50-89 filters', async () => {
    const result = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']});
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 50, max: 89});
      return {filterCount: grok.shell.tv.dataFrame.filter.trueCount};
    });
    expect(result.filterCount).toBe(8);
  });

  // Step 9: Save layout
  let layoutId: string;
  await softStep('Save layout', async () => {
    layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });
    expect(layoutId).toBeTruthy();
  });

  // Step 10: Close filter panel
  await softStep('Close filter panel', async () => {
    await page.evaluate(() => {
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      let el = filterViewer;
      while (el && !el.classList.contains('panel-base')) el = el.parentElement;
      el.querySelector('.grok-font-icon-close').click();
    });
    await expect(page.locator('[name="viewer-Filters"]')).not.toBeVisible({timeout: 5000});
  });

  // Step 11: Apply saved layout — verify
  await softStep('Apply saved layout and verify restored state', async () => {
    const result = await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 3000));
      return {filterCount: grok.shell.tv.dataFrame.filter.trueCount};
    }, layoutId!);
    expect(result.filterCount).toBe(8);
  });

  // Step 12: Remove all filters via hamburger menu
  await softStep('Remove all filters via hamburger menu', async () => {
    await page.evaluate(() => {
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      let el = filterViewer;
      while (el && !el.classList.contains('panel-base')) el = el.parentElement;
      const titleBar = el.querySelector('.panel-titlebar');
      titleBar.querySelector('[name="icon-font-icon-menu"]').click();
    });
    await page.locator('.d4-menu-popup').waitFor({timeout: 3000});
    await page.evaluate(() => {
      const popup = document.querySelector('.d4-menu-popup');
      const items = popup.querySelectorAll('.d4-menu-item-label');
      for (const item of items) {
        if (item.textContent.trim() === 'Remove All') {
          item.closest('.d4-menu-item').click();
          return;
        }
      }
    });
    await page.waitForTimeout(500);
    const result = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(result).toBe(5850);
  });

  // Step 13: Close filter panel
  await softStep('Close filter panel after remove all', async () => {
    await page.evaluate(() => {
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      let el = filterViewer;
      while (el && !el.classList.contains('panel-base')) el = el.parentElement;
      el.querySelector('.grok-font-icon-close').click();
    });
    await expect(page.locator('[name="viewer-Filters"]')).not.toBeVisible({timeout: 5000});
  });

  // Step 14: Reopen filter panel — Combined Boolean should auto-add
  await softStep('Reopen filter panel — Combined Boolean auto-added', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1000));
    });
    const hasBool = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      return Array.from(fg.filters).some(f => f.filterType === 'bool-columns');
    });
    expect(hasBool).toBe(true);
  });

  // Step 15: Apply saved layout again — verify
  await softStep('Apply saved layout after reopen and verify', async () => {
    const result = await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 3000));
      return {filterCount: grok.shell.tv.dataFrame.filter.trueCount};
    }, layoutId!);
    expect(result.filterCount).toBe(8);
  });

  // Cleanup
  await page.evaluate(async (id) => {
    const saved = await grok.dapi.layouts.find(id);
    await grok.dapi.layouts.delete(saved);
  }, layoutId!);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

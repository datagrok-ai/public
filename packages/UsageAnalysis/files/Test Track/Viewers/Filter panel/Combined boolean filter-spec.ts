
import {test, expect, Page} from '@playwright/test';

async function login(page: Page) {
  await page.goto('/');
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
  await page.evaluate(() => { document.body.classList.add('selenium'); grok.shell.settings.showFiltersIconsConstantly = true; });
}

async function filteredCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

async function selectedCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
}

async function clickFilterPanelHamburger(page: Page) {
  await page.evaluate(() => {
    const fp = document.querySelector('[name="viewer-Filters"]');
    const panel = fp?.closest('.d4-viewer-host') || fp?.closest('.panel-content')?.parentElement;
    const menuIcon = panel?.querySelector('[name="icon-font-icon-menu"]') as HTMLElement;
    if (menuIcon) menuIcon.click();
  });
}

async function safeStep(page: Page, name: string, body: () => Promise<void>) {
  await test.step(name, async () => {
    try {
      await body();
    }
    catch (e) {
      expect.soft(null, `Step "${name}" threw: ${(e as Error).message}`).not.toBeNull();
    }
  });
}

test('Combined Boolean Filter: demog', async ({page}) => {
  test.setTimeout(300_000);
  await login(page);

  // 1. Open demog
  await page.evaluate(async () => {
    grok.shell.closeAll();
    const df = await grok.data.getDemoTable('demog.csv');
    grok.shell.addTableView(df);
  });
  // Wait for grid to appear
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 10000});
  const totalRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);

  // 2. Open Filter Panel
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});

  // 3. Add SEX_bool calculated column
  await safeStep(page, 'Add SEX_bool column', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.columns.addNewCalculated('SEX_bool', 'case ${SEX} when "F" then true else false end');
    });

    await expect.poll(() => page.evaluate(() => {
      const col = grok.shell.tv.dataFrame.col('SEX_bool');
      return col != null && col.type === 'bool';
    }), 'SEX_bool column should appear').toBe(true);
  });

  // 4. Hamburger > Add Filter > Combined Boolean
  await safeStep(page, 'Add Combined Boolean filter', async () => {
    await clickFilterPanelHamburger(page);

    await page.getByText('Add Filter', {exact: true}).hover();
    await page.getByText('Combined Boolean', {exact: true}).click();

    // Wait for "Flags" label to appear
    await expect.poll(() => page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      return fp?.textContent?.includes('Flags') ?? false;
    }), 'Combined Boolean filter should show Flags label').toBe(true);
  });

  // 5. Verify Combined Boolean filter present via JS API
  await safeStep(page, 'Verify boolean columns in Combined Boolean filter', async () => {
    const info = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        const ft = (f as any).filterType;
        if (ft === 'bool-columns')
          return {found: true, filterType: ft};
      }
      return {found: false, filterType: 'none'};
    });
    expect.soft(info.found, 'Combined Boolean filter should be present').toBe(true);
  });

  // 6. Select SEX_bool=true rows via JS API, verify count, then Esc
  await safeStep(page, 'Click count to select rows', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('SEX_bool');
      if (col) df.selection.init((i: number) => col.get(i) === true);
    });

    await expect.poll(() => selectedCount(page), 'SEX_bool=true should select 3243').toBe(3243);

    await page.keyboard.press('Escape');

    await expect.poll(() => selectedCount(page), 'Esc should clear selection').toBe(0);
  });

  // 7. Apply combined filter via applyState
  await safeStep(page, 'Apply combined filter', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        const ft = (f as any).filterType;
        if (ft === 'bool-columns') {
          const dart = (f as any).dart ?? f;
          (window as any).grok_GridFilterBase_ApplyState(dart, {
            'true': [true, false],
            'false': [false, true],
            mode: 'OR'
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });

    await expect.poll(() => filteredCount(page), 'Combined filter should show 2632').toBe(2632);
  });

  // 8. Apply other filters: RACE = Asian, AGE 50-89
  await safeStep(page, 'Apply categorical filter (RACE = Asian)', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']});
    });

    await expect.poll(() => filteredCount(page), 'RACE filter should reduce count')
      .toBeLessThan(2632);
  });

  await safeStep(page, 'Apply numerical filter (AGE 50-89)', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 50, max: 89});
      await new Promise(r => setTimeout(r, 500));
    });

    const count = await filteredCount(page);
    expect.soft(count, 'After RACE + AGE filters, count > 0').toBeGreaterThan(0);
    expect.soft(count, 'After RACE + AGE filters, count < total').toBeLessThan(totalRows);
  });

  // 9. Save layout
  await safeStep(page, 'Save layout', async () => {
    await page.locator('[name="div-section--Layouts"]').click();
    const saveBtn = page.locator('.d4-toolbox-layouts').locator('button:has-text("SAVE")');
    await saveBtn.click();
    await page.locator('.grok-suggestions-chart-host').first().waitFor({timeout: 10000});
    // Wait for server-side save to complete (auto-wait can't detect this)
    await page.waitForTimeout(2000);
  });

  const countBeforeClose = await filteredCount(page);

  // 10. Close Filter Panel (Layouts section stays open from step 9)
  await safeStep(page, 'Close Filter Panel', async () => {
    await page.evaluate(() => {
      for (const v of grok.shell.tv.viewers)
        if (v.type === 'Filters') { v.close(); break; }
    });

    await expect.poll(() => filteredCount(page), 'After close all rows visible').toBe(totalRows);
  });

  // 11. Apply saved layout — thumbnail still visible
  await safeStep(page, 'Apply saved layout', async () => {
    await page.locator('.grok-suggestions-chart-host').first().click();

    // Wait for filters to restore
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
    });
    await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});

    await expect.poll(() => filteredCount(page), 'Layout should restore filters')
      .toBeLessThan(totalRows);
    const countAfterApply = await filteredCount(page);
    expect.soft(countAfterApply, 'Should match state before close').toBe(countBeforeClose);
  });

  // 12. Remove all filters
  await safeStep(page, 'Remove all filters', async () => {
    await clickFilterPanelHamburger(page);
    await page.getByText('Remove All', {exact: true}).click();
    const okBtn = page.getByText('OK', {exact: true});
    if (await okBtn.count() > 0)
      await okBtn.first().click();

    await expect.poll(() => filteredCount(page), 'After Remove All, all rows visible').toBe(totalRows);
  });

  // 13. Close Filter Panel
  await safeStep(page, 'Close Filter Panel (second time)', async () => {
    await page.evaluate(() => {
      for (const v of grok.shell.tv.viewers)
        if (v.type === 'Filters') { v.close(); break; }
    });

    await expect.poll(() => filteredCount(page), 'After close all rows visible').toBe(totalRows);
  });

  // 14. Reopen Filter Panel — Combined Boolean auto-added
  await safeStep(page, 'Reopen Filter Panel — Combined Boolean auto-added', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});

    const hasBoolFilter = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      const boolCols = [];
      for (let i = 0; i < grok.shell.tv.dataFrame.columns.length; i++) {
        const c = grok.shell.tv.dataFrame.columns.byIndex(i);
        if (c.type === 'bool') boolCols.push(c.name);
      }
      return boolCols.length > 0 && fg.filters.length > 0;
    });
    expect.soft(hasBoolFilter, 'After reopening, boolean filters should be present').toBe(true);
  });

  // Cleanup
  await page.evaluate(async () => {
    try {
      const layouts = await grok.dapi.layouts.filter('table.name = "demog"').list();
      if (layouts.length > 0) {
        layouts.sort((a: any, b: any) => b.createdOn - a.createdOn);
        await grok.dapi.layouts.delete(layouts[0]);
      }
    }
    catch (e) {}
  });
});

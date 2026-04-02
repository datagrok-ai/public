import {test, expect, Page} from '@playwright/test';

async function login(page: Page) {
  await page.goto('/');
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
  await page.evaluate(() => { document.body.classList.add('selenium'); grok.shell.settings.showFiltersIconsConstantly = true; });
  await page.waitForTimeout(3000);
}

async function filteredCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

async function clickFilterPanelHamburger(page: Page) {
  await page.evaluate(() => {
    const fp = document.querySelector('[name="viewer-Filters"]');
    const panel = fp?.closest('.d4-viewer-host') || fp?.closest('.panel-content')?.parentElement;
    const menuIcon = panel?.querySelector('[name="icon-font-icon-menu"]') as HTMLElement;
    if (menuIcon) menuIcon.click();
  });
  await page.waitForTimeout(500);
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

/** Click a hierarchical filter checkbox by value text */
async function clickHierCheckbox(page: Page, valueName: string) {
  await page.evaluate((name) => {
    const nodes = document.querySelectorAll('.d4-tree-view-node');
    for (const node of nodes) {
      const val = node.querySelector('.d4-hierarchical-filter-caption-value');
      if (val?.textContent?.trim() === name) {
        (node.querySelector('.d4-hierarchical-filter-checkbox') as HTMLElement)?.click();
        break;
      }
    }
  }, valueName);
  await page.waitForTimeout(500);
}

test('Hierarchical Filter: demog', async ({page}) => {
  test.setTimeout(300_000);
  await login(page);

  // 1. Open demog
  await page.evaluate(async () => {
    grok.shell.closeAll();
    const df = await grok.data.getDemoTable('demog.csv');
    grok.shell.addTableView(df);
  });
  await page.waitForTimeout(3000);
  const totalRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);

  // 2. Open Filter Panel
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.waitForTimeout(2000);

  // 3. Hamburger > Add Filter > Hierarchical
  await safeStep(page, 'Add Hierarchical filter', async () => {
    await clickFilterPanelHamburger(page);
    await page.waitForTimeout(500);

    await page.getByText('Add Filter', {exact: true}).hover();
    await page.waitForTimeout(500);

    await page.getByText('Hierarchical', {exact: true}).click();
    await page.waitForTimeout(2000);

    const hasFilter = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      return fp?.querySelector('.svg-tree') != null;
    });
    expect.soft(hasFilter, 'Hierarchical filter should appear with .svg-tree icon').toBe(true);
  });

  // 4. Select columns SEX, RACE, SEVERITY via applyState
  await safeStep(page, 'Select columns: SEX, RACE, SEVERITY', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if ((f as any).filterType === 'hierarchical') {
          const dart = (f as any).dart ?? f;
          (window as any).grok_GridFilterBase_ApplyState(dart, {
            type: 'hierarchical',
            active: true,
            colNames: ['SEX', 'RACE', 'SEVERITY'],
            allEnabled: true
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    await page.waitForTimeout(2000);

    // Verify header shows column path
    const headerText = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      const label = fp?.querySelector('.d4-filter-column-name');
      return label?.textContent?.trim() ?? '';
    });
    expect.soft(headerText, 'Header should contain SEX').toContain('SEX');
    expect.soft(headerText, 'Header should contain SEVERITY').toContain('SEVERITY');
  });

  // 5. Expand F → Caucasian → SEVERITY values visible
  await safeStep(page, 'Expand F then Caucasian', async () => {
    // Expand "F" using name= selector
    await page.locator('[name="tree-expander-F"]').click();
    await page.waitForTimeout(1000);

    // Expand "Caucasian" by finding the node by value text
    await page.evaluate(() => {
      const nodes = document.querySelectorAll('.d4-tree-view-node');
      for (const node of nodes) {
        const val = node.querySelector('.d4-hierarchical-filter-caption-value');
        if (val?.textContent?.trim() === 'Caucasian') {
          const tri = node.querySelector('.d4-tree-view-tri') as HTMLElement;
          if (tri) tri.click();
          break;
        }
      }
    });
    await page.waitForTimeout(1000);

    // Verify SEVERITY values visible
    const hasSeverity = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      const values = fp?.querySelectorAll('.d4-hierarchical-filter-caption-value');
      const texts = Array.from(values || []).map(v => v.textContent?.trim());
      return texts.includes('None') && texts.includes('High') && texts.includes('Low');
    });
    expect.soft(hasSeverity, 'SEVERITY values should be visible after expanding').toBe(true);
  });

  // 6. Select Caucasian → all children checked, verify row count
  await safeStep(page, 'Select Caucasian — children follow', async () => {
    await clickHierCheckbox(page, 'Caucasian');
    await page.waitForTimeout(1000);

    const count = await filteredCount(page);
    expect.soft(count, 'Only Caucasian F rows should pass').toBeLessThan(totalRows);
    expect.soft(count, 'Some Caucasian rows should match').toBeGreaterThan(0);
  });

  // 7. Deselect Low and Medium → parent indeterminate, count decreases
  await safeStep(page, 'Deselect Low and Medium', async () => {
    const countBefore = await filteredCount(page);

    await clickHierCheckbox(page, 'Low');
    await clickHierCheckbox(page, 'Medium');
    await page.waitForTimeout(1000);

    const countAfter = await filteredCount(page);
    expect.soft(countAfter, 'Count should decrease after unchecking Low/Medium').toBeLessThan(countBefore);
    expect.soft(countAfter, 'Some rows should remain').toBeGreaterThan(0);

    // Verify Caucasian is indeterminate (char code 0xf146 = minus-square)
    const isIndeterminate = await page.evaluate(() => {
      const nodes = document.querySelectorAll('.d4-tree-view-node');
      for (const node of nodes) {
        const val = node.querySelector('.d4-hierarchical-filter-caption-value');
        if (val?.textContent?.trim() === 'Caucasian') {
          const sub = node.querySelector('.d4-hierarchical-filter-checkbox-substitute');
          return sub?.textContent?.charCodeAt(0) === 0xf146;
        }
      }
      return false;
    });
    expect.soft(isIndeterminate, 'Caucasian should be in indeterminate state').toBe(true);
  });

  // 8. Rearrange columns: SEVERITY before RACE via applyState
  await safeStep(page, 'Rearrange columns: SEX / SEVERITY / RACE', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of fg.filters) {
        if ((f as any).filterType === 'hierarchical') {
          const dart = (f as any).dart ?? f;
          (window as any).grok_GridFilterBase_ApplyState(dart, {
            type: 'hierarchical', active: true,
            colNames: ['SEX', 'SEVERITY', 'RACE'], allEnabled: true
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    await page.waitForTimeout(2000);

    const headerText = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      const label = fp?.querySelector('.d4-filter-column-name');
      return label?.textContent?.trim() ?? '';
    });
    const sevIdx = headerText.indexOf('SEVERITY');
    const raceIdx = headerText.indexOf('RACE');
    expect.soft(sevIdx, 'SEVERITY should be in header').toBeGreaterThan(-1);
    if (sevIdx > -1 && raceIdx > -1)
      expect.soft(sevIdx, 'SEVERITY should appear before RACE').toBeLessThan(raceIdx);
  });

  // 9. Collaborative filtering: F → None, then DIS_POP = AS + Indigestion → 339
  await safeStep(page, 'Collaborative filtering', async () => {
    // Expand F
    await page.locator('[name="tree-expander-F"]').click();
    await page.waitForTimeout(500);

    // Click text "None" (not checkbox) to select ONLY that value
    // Clicking the label text selects only that node (deselects everything else)
    await page.evaluate(() => {
      const nodes = document.querySelectorAll('.d4-tree-view-node');
      for (const node of nodes) {
        const val = node.querySelector('.d4-hierarchical-filter-caption-value');
        if (val?.textContent?.trim() === 'None') {
          (val as HTMLElement).click();
          break;
        }
      }
    });
    await page.waitForTimeout(1000);

    const countFNone = await filteredCount(page);
    expect.soft(countFNone, 'F + None should give ~1815 rows').toBeGreaterThan(0);
    expect.soft(countFNone, 'F + None should be less than total').toBeLessThan(totalRows);

    // Apply DIS_POP categorical filter: AS and Indigestion
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'DIS_POP', selected: ['AS', 'Indigestion']});
    });
    await page.waitForTimeout(1500);

    const countWithDisPop = await filteredCount(page);
    expect.soft(countWithDisPop, 'With DIS_POP AS+Indigestion, count should be 339').toBe(339);
  });

  // 10. Save layout
  await safeStep(page, 'Save layout', async () => {
    await page.locator('[name="div-section--Layouts"]').click();
    await page.waitForTimeout(500);
    const saveBtn = page.locator('.d4-toolbox-layouts').locator('button:has-text("SAVE")');
    await saveBtn.click();
    await page.locator('.grok-suggestions-chart-host').first().waitFor({timeout: 10000});
    await page.waitForTimeout(1000);
  });

  const countBeforeClose = await filteredCount(page);

  // 11. Close Filter Panel via close (X) icon on viewer title bar
  await safeStep(page, 'Close Filter Panel', async () => {
    await page.evaluate(() => {
      // Close filter panel viewer — removes all filters
      for (const v of grok.shell.tv.viewers) {
        if (v.type === 'Filters') { v.close(); break; }
      }
    });
    await page.waitForTimeout(2000);

    const count = await filteredCount(page);
    expect.soft(count, 'After closing filter panel all rows visible').toBe(totalRows);
  });

  // 12. Apply saved layout — verify restored
  await safeStep(page, 'Apply saved layout', async () => {
    await page.locator('.grok-suggestions-chart-host').first().click();
    await page.waitForTimeout(5000);

    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
    });

    const countAfterApply = await filteredCount(page);
    expect.soft(countAfterApply, 'Applied layout should restore filters').toBeLessThan(totalRows);
    expect.soft(countAfterApply, 'Applied layout should match state before close').toBe(countBeforeClose);
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
  await page.waitForTimeout(1000);
});

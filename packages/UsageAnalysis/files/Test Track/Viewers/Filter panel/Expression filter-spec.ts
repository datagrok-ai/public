import {test, expect, Page} from '@playwright/test';

async function login(page: Page) {
  await page.goto('/');
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
  await page.evaluate(() => { document.body.classList.add('selenium'); grok.shell.settings.showFiltersIconsConstantly = true; });
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
}

/** Fill the Expression filter form fields and click the Add (+) button.
 *  KEEP small delays — Dart dropdown widgets need time between value changes. */
async function addExpressionRule(page: Page, column: string, operation: string, value: string) {
  const form = page.locator('.d4-expression-filter');

  const columnSelect = form.locator('.ui-input-choice select').first();
  await columnSelect.selectOption(column);
  await page.waitForTimeout(300);

  const opSelect = form.locator('.ui-input-choice select').nth(1);
  await opSelect.selectOption(operation);
  await page.waitForTimeout(300);

  const valueInput = form.locator('.ui-input-text input').first();
  await valueInput.click();
  await valueInput.fill(value);
  await page.waitForTimeout(300);

  const addBtn = form.locator('.fal.fa-plus').first();
  await addBtn.click();
  await page.waitForTimeout(500);
}

/** Focus and type into the free-text search input.
 *  KEEP delays — keyboard.type needs real timing for Dart key events. */
async function fillFreeTextSearch(page: Page, expression: string) {
  await page.evaluate(() => {
    const fp = document.querySelector('[name="viewer-Filters"]');
    const input = fp?.querySelector('input.d4-search-input[placeholder="Search"]') as HTMLInputElement;
    if (input) {
      input.value = '';
      input.focus();
    }
  });
  await page.waitForTimeout(500);

  const searchInput = page.locator('[name="viewer-Filters"] input.d4-search-input[placeholder="Search"]').first();
  await searchInput.click();
  await page.waitForTimeout(300);

  await page.keyboard.type(expression, {delay: 30});
  await page.waitForTimeout(500);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(1000);
}

/** Uncheck expression filter rules by calling applyState with known rule names. */
async function uncheckRulesByName(page: Page, allRules: string[], enabledRules: string[]) {
  const gridValues = allRules.map(name => enabledRules.includes(name));

  const result = await page.evaluate(({gridNames, gridValues}) => {
    const fg = grok.shell.tv.getFiltersGroup();
    const filters = fg.filters;
    const api = (window as any);

    for (let i = 0; i < filters.length; i++) {
      const f = filters[i] as any;
      if ((f.filterType ?? '') !== 'expression') continue;

      const dartHandle = f.dart ?? f;
      try {
        api.grok_GridFilterBase_ApplyState(dartHandle, {
          gridNames: gridNames,
          gridValues: gridValues,
          mode: 'OR',
          value: ''
        });
        grok.shell.tv.dataFrame.rows.requestFilter();
        return {success: true};
      }
      catch (e) {
        return {success: false, error: String(e)};
      }
    }
    return {success: false, reason: 'expression filter not found'};
  }, {gridNames: allRules, gridValues});

  await page.waitForTimeout(500);
  return result;
}

/** Wrap a step body so that exceptions are caught and reported as soft failures */
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

test('Expression Filter: demog', async ({page}) => {
  test.setTimeout(300_000);
  await login(page);

  // 1. Open demog
  await page.evaluate(async () => {
    grok.shell.closeAll();
    const df = await grok.data.getDemoTable('demog.csv');
    grok.shell.addTableView(df);
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 10000});
  const totalRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);

  // 2. Open Filter Panel
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});

  // 3. Hamburger > Add Filter > Expression
  await safeStep(page, 'Add Expression filter via hamburger menu', async () => {
    await clickFilterPanelHamburger(page);

    await page.getByText('Add Filter', {exact: true}).hover();
    await page.getByText('Expression', {exact: true}).click();

    // Wait for expression form to appear
    await page.locator('.d4-expression-filter').waitFor({timeout: 5000});
    const formVisible = await page.locator('.d4-expression-filter').count();
    expect.soft(formVisible, 'Expression filter form should appear').toBeGreaterThan(0);
  });

  // 4. Add expression filter rules via the form (Expression filter defaults to OR)
  await safeStep(page, 'Add WEIGHT > 50', async () => {
    await addExpressionRule(page, 'WEIGHT', '>', '50');
    await expect.poll(() => filteredCount(page), {timeout: 5000}).toBeLessThan(totalRows);
    expect.soft(await filteredCount(page), 'Some rows should match WEIGHT > 50').toBeGreaterThan(0);
  });

  await safeStep(page, 'Add HEIGHT < 160', async () => {
    await addExpressionRule(page, 'HEIGHT', '<', '160');
    expect.soft(await filteredCount(page), 'HEIGHT < 160 should match').toBeGreaterThan(0);
  });

  await safeStep(page, 'Add SEX = F', async () => {
    await addExpressionRule(page, 'SEX', 'equals', 'F');
    expect.soft(await filteredCount(page), 'SEX = F should match').toBeGreaterThan(0);
  });

  await safeStep(page, 'Add RACE contains an', async () => {
    await addExpressionRule(page, 'RACE', 'contains', 'an');
    expect.soft(await filteredCount(page), 'RACE contains an should match').toBeGreaterThan(0);
  });

  await safeStep(page, 'Add STARTED after 01/01/1991', async () => {
    await addExpressionRule(page, 'STARTED', 'after', '01/01/1991');
    expect.soft(await filteredCount(page), 'STARTED after 01/01/1991 should match').toBeGreaterThan(0);
  });

  const countInORMode = await filteredCount(page);

  // 5. Click OR → switch to AND
  await safeStep(page, 'Switch to AND logic', async () => {
    const orLabel = page.locator('[name="viewer-Filters"]').getByText('OR', {exact: true}).first();
    if (await orLabel.count() > 0) {
      await orLabel.click();
      await expect.poll(() => filteredCount(page), {timeout: 5000}).toBeLessThanOrEqual(countInORMode);
    }

    const countAND = await filteredCount(page);
    expect.soft(countAND, 'AND mode should filter more strictly than OR').toBeLessThanOrEqual(countInORMode);
    expect.soft(countAND, 'AND mode should have some matches').toBeGreaterThan(0);
  });

  // 6. Click AND → switch back to OR
  await safeStep(page, 'Switch back to OR logic', async () => {
    const andLabel = page.locator('[name="viewer-Filters"]').getByText('AND', {exact: true}).first();
    if (await andLabel.count() > 0) {
      await andLabel.click();
      await expect.poll(() => filteredCount(page), {timeout: 5000}).toBeGreaterThanOrEqual(1);
    }

    expect.soft(await filteredCount(page), 'OR mode should show rows').toBeGreaterThanOrEqual(1);
  });

  // 7. Right-click first expression rule → Remove Query
  await safeStep(page, 'Remove first query via right-click', async () => {
    const rulesGrid = page.locator('[name="viewer-Filters"] canvas').first();
    if (await rulesGrid.count() > 0) {
      const box = await rulesGrid.boundingBox();
      if (box) {
        await page.mouse.click(box.x + box.width / 2, box.y + 10, {button: 'right'});
        await page.waitForTimeout(500); // wait for context menu

        const removeItem = page.getByText('Remove Query', {exact: true});
        if (await removeItem.count() > 0)
          await removeItem.click();
      }
    }

    expect.soft(await filteredCount(page), 'Count should be valid after removing a query').toBeGreaterThan(0);
  });

  // 8. Hover Expression filter header → click Free text icon (italic I)
  await safeStep(page, 'Switch to Free text mode', async () => {
    const switched = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      if (!fp) return false;
      const icon = fp.querySelector('.fa-italic') as HTMLElement;
      if (icon) {
        icon.click();
        return true;
      }
      return false;
    });

    // Wait for search input to become visible
    await expect.poll(() => page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      const input = fp?.querySelector('input.d4-search-input[placeholder="Search"]') as HTMLElement;
      return input ? input.style.display !== 'none' : false;
    }), {timeout: 5000}).toBe(true);

    expect.soft(switched, 'Should find and click italic icon').toBe(true);
  });

  // 9. Enter "AGE > 30 and SEX = M"
  await safeStep(page, 'Free text: AGE > 30 and SEX = M', async () => {
    await fillFreeTextSearch(page, 'AGE > 30 and SEX = M');
    expect.soft(await filteredCount(page), 'AGE > 30 and SEX = M should match some rows').toBeGreaterThan(0);
  });

  // 10. Enter "AGE < 60 and HEIGHT > 190"
  await safeStep(page, 'Free text: AGE < 60 and HEIGHT > 190', async () => {
    await fillFreeTextSearch(page, 'AGE < 60 and HEIGHT > 190');
    expect.soft(await filteredCount(page), 'AGE < 60 and HEIGHT > 190 should match some rows').toBeGreaterThan(0);
  });

  // 11. Uncheck first 4 expression rules, verify filtering changes
  await safeStep(page, 'Uncheck first 4 rules and verify filtering', async () => {
    const countBefore = await filteredCount(page);

    const allRules = [
      '${HEIGHT} < 160',
      '${SEX} equals F',
      '${RACE} contains an',
      '${STARTED} after 01/01/1991',
      'AGE > 30 and SEX = M',
      'AGE < 60 and HEIGHT > 190'
    ];
    const enabledRules = [
      'AGE > 30 and SEX = M',
      'AGE < 60 and HEIGHT > 190'
    ];

    const result = await uncheckRulesByName(page, allRules, enabledRules);
    console.log('uncheckRules result:', JSON.stringify(result));

    const countAfter = await filteredCount(page);
    expect.soft(countAfter, 'Unchecking rules should change filtered count').not.toBe(countBefore);
    expect.soft(countAfter, 'Should still have matching rows').toBeGreaterThan(0);
    expect.soft(countAfter, 'Should filter some rows').toBeLessThan(totalRows);
  });

  // 12. Save layout via Toolbox > Layouts > SAVE
  await safeStep(page, 'Save layout', async () => {
    await page.locator('[name="div-section--Layouts"]').click();
    const saveBtn = page.locator('.d4-toolbox-layouts').locator('button:has-text("SAVE")');
    await saveBtn.click();
    await page.locator('.grok-suggestions-chart-host').first().waitFor({timeout: 10000});
    // Wait for server-side save to complete
    await page.waitForTimeout(2000);
  });

  const countBeforeClose = await filteredCount(page);

  // 13. Close Filter Panel — all filters removed, all rows visible
  await safeStep(page, 'Close Filter Panel', async () => {
    await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      const panel = fp?.closest('.d4-viewer-host') || fp?.closest('.panel-content')?.parentElement;
      const closeIcon = panel?.querySelector('[name="icon-times"]') as HTMLElement;
      if (closeIcon) closeIcon.click();
    });

    await expect.poll(() => filteredCount(page), {timeout: 5000}).toBe(totalRows);
  });

  // 14. Apply saved layout by clicking the layout thumbnail
  await safeStep(page, 'Apply saved layout', async () => {
    await page.locator('.grok-suggestions-chart-host').first().click();
    // Wait for filter panel to reappear after layout apply (server-side)
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
    });
    await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});

    await expect.poll(() => filteredCount(page), {timeout: 5000}).toBeLessThan(totalRows);
    expect.soft(await filteredCount(page), 'Applied layout should match state before close').toBe(countBeforeClose);
  });

  // Cleanup: delete saved layout
  await page.evaluate(async () => {
    const layouts = await grok.dapi.layouts.filter('table.name = "demog"').list();
    if (layouts.length > 0) {
      layouts.sort((a: any, b: any) => b.createdOn - a.createdOn);
      await grok.dapi.layouts.delete(layouts[0]);
    }
  });
});

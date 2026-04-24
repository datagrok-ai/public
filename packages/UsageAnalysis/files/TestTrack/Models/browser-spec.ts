import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Predictive Models Browser — navigate, search, filter templates, multi-select Compare', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    (window as any).grok.shell.windows.showBrowse = true;
  });

  await softStep('1. Go to Browse > Platform > Predictive models', async () => {
    const node = page.locator('[name="tree-Platform---Predictive-models"]');
    if (!(await node.isVisible({timeout: 2000}))) {
      // Platform branch is collapsed on a fresh tree — expand it first.
      await page.locator('[name="tree-Platform"] .d4-tree-view-tri').click();
      await node.waitFor({state: 'visible', timeout: 15000});
    }
    await node.click();
    await page.waitForTimeout(1500);
    const path = await page.evaluate(() => (window as any).grok.shell.v?.path);
    expect(path).toBe('/models?');
    await expect(page.locator('input[placeholder*="Search models"]')).toBeVisible({timeout: 15000});
  });

  await softStep('2. Type TestDemog in the search field', async () => {
    const search = page.locator('input[placeholder*="Search models"]').first();
    await search.fill('TestDemog');
    await page.waitForTimeout(1500);
    const model = page.locator('text=TestDemog').first();
    await expect(model).toBeVisible({timeout: 8000});
  });

  await softStep('3. On the Context Panel, check all tabs for the model', async () => {
    const model = page.locator('text=TestDemog').first();
    if (!(await model.isVisible({timeout: 2000})))
      throw new Error('PREREQUISITE: TestDemog model not present (run Models/train.md first)');
    await model.click();
    await page.waitForTimeout(800);
    const panel = page.locator('.grok-prop-panel, [class*="context-panel"]').first();
    await expect(panel).toBeVisible({timeout: 5000});
    for (const tabName of ['Details', 'Performance', 'Activity', 'Sharing']) {
      const tab = panel.locator(`text=${tabName}`).first();
      if (await tab.isVisible({timeout: 1500}))
        await tab.click();
    }
  });

  await softStep('4. Click the Filter templates icon near the search field', async () => {
    await page.locator('i[name="icon-filter"]').first().click();
    await page.waitForTimeout(800);
    await expect(page.locator('text=Quick Filters').first()).toBeVisible({timeout: 5000});
  });

  await softStep('5. Select multiple models (CTRL+click)', async () => {
    const items = page.locator('.grok-gallery-grid-item, .d4-item-grid .d4-item');
    const count = await items.count();
    if (count < 2)
      throw new Error(`PREREQUISITE: need >=2 models for multi-select, found ${count}`);
    await items.nth(0).click();
    await items.nth(1).click({modifiers: ['Control']});
    await page.waitForTimeout(500);
  });

  await softStep('6. Open the Commands tab on the Context Pane and click Compare', async () => {
    const count = await page.locator('.grok-gallery-grid-item, .d4-item-grid .d4-item').count();
    if (count < 2)
      throw new Error(`PREREQUISITE: need >=2 models for Compare, found ${count}`);
    const panel = page.locator('.grok-prop-panel, [class*="context-panel"]').first();
    const commandsTab = panel.locator('text=Commands').first();
    await commandsTab.click();
    await page.waitForTimeout(400);
    await panel.locator('text=Compare').first().click();
    await page.waitForTimeout(1000);
  });

  if (stepErrors.length > 0) {
    throw new Error(`[${stepErrors.length} step(s) failed]\n` +
      stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});

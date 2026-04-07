import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai/';

async function waitForGrok(page: Page) {
  await page.waitForFunction(() => typeof (window as any).grok !== 'undefined', { timeout: 30000 });
}

async function closeAll(page: Page) {
  await page.evaluate(() => (window as any).grok.shell.closeAll());
}

test.describe('Predictive Models Browser', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE_URL);
    await waitForGrok(page);
    await closeAll(page);
  });

  test('Browse predictive models, search, context panel tabs', async ({ page }) => {
    // Prerequisite: TestDemog model must exist

    // Step 1: Navigate to Browse > Platform > Predictive models
    await page.goto(`${BASE_URL}models`);
    await page.waitForTimeout(2000);

    const modelsView = page.locator('.d4-item-grid, [class*="browser"], [class*="grid"]').first();
    await expect(modelsView).toBeVisible({ timeout: 10000 });

    // Step 2: Search for TestDemog
    const searchBox = page.locator('input[placeholder*="Search"], input[placeholder*="search"], .d4-search input').first();
    if (await searchBox.isVisible()) {
      await searchBox.fill('TestDemog');
      await page.waitForTimeout(1000);
    }

    const modelItem = page.locator('text=TestDemog').first();
    await expect(modelItem).toBeVisible({ timeout: 5000 });
    await modelItem.click();
    await page.waitForTimeout(1000);

    // Step 3: Check Context Panel tabs
    const contextPanel = page.locator('.grok-prop-panel, .d4-info-panel, [class*="context-panel"]').first();
    await expect(contextPanel).toBeVisible({ timeout: 5000 });

    // Check for expected tabs
    for (const tabName of ['Details', 'Performance', 'Activity', 'Sharing']) {
      const tab = contextPanel.locator(`text=${tabName}`).first();
      if (await tab.isVisible()) await tab.click();
      await page.waitForTimeout(300);
    }

    // Step 4: Check Filter templates icon
    const filterIcon = page.locator('[class*="filter"], button[title*="Filter"], .d4-filter-icon').first();
    if (await filterIcon.isVisible()) {
      await filterIcon.click();
      await page.waitForTimeout(500);
    }

    await page.screenshot({ path: 'models-browser-single.png' });
  });

  test('Multi-select and Compare (requires ≥2 models)', async ({ page }) => {
    await page.goto(`${BASE_URL}models`);
    await page.waitForTimeout(2000);

    const modelItems = page.locator('.d4-item-grid .d4-item, [class*="card"], [class*="item"]');
    const count = await modelItems.count();

    if (count < 2) {
      test.skip(); // Not enough models for multi-select test
      return;
    }

    // CTRL+click to multi-select
    await modelItems.nth(0).click();
    await modelItems.nth(1).click({ modifiers: ['Control'] });
    await page.waitForTimeout(500);

    // Open Commands tab
    const contextPanel = page.locator('.grok-prop-panel, [class*="context-panel"]').first();
    const commandsTab = contextPanel.locator('text=Commands').first();
    if (await commandsTab.isVisible()) {
      await commandsTab.click();
      await page.waitForTimeout(500);
      // Click Compare
      const compareBtn = contextPanel.locator('text=Compare').first();
      if (await compareBtn.isVisible()) await compareBtn.click();
    }

    await page.screenshot({ path: 'models-browser-compare.png' });
  });
});

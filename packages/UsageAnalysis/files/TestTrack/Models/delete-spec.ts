import { test, expect, Page } from '@playwright/test';
import {specTestOptions} from '../spec-login';

test.use(specTestOptions);

const BASE_URL = 'https://public.datagrok.ai/';

async function waitForGrok(page: Page) {
  await page.waitForFunction(() => typeof (window as any).grok !== 'undefined', { timeout: 30000 });
}

test.describe('Delete predictive model', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE_URL);
    await waitForGrok(page);
  });

  test('Delete model from Browse > Platform > Models', async ({ page }) => {
    // Prerequisite: a model from previous steps must exist

    // Step 1: Navigate to Browse > Platform > Models
    await page.goto(`${BASE_URL}models`);
    await page.waitForTimeout(2000);

    const modelsView = page.locator('.d4-item-grid, [class*="browser"]').first();
    await expect(modelsView).toBeVisible({ timeout: 10000 });

    // Step 2: Find the model
    const modelItem = page.locator('.d4-item, [class*="card"]').first();
    const modelCount = await modelItem.count();
    if (modelCount === 0) {
      test.skip(); // No models to delete
      return;
    }

    const firstModel = page.locator('.d4-item, [class*="card"]').first();

    // Step 3: Right-click and select Delete
    await firstModel.click({ button: 'right' });
    await page.waitForTimeout(500);

    const deleteMenuItem = page.locator('[class*="menu"] text=Delete, .d4-menu-item:has-text("Delete")').first();
    await expect(deleteMenuItem).toBeVisible({ timeout: 5000 });
    await deleteMenuItem.click();
    await page.waitForTimeout(500);

    // Step 4: Confirm deletion
    const confirmDialog = page.locator('.d4-dialog, [class*="dialog"]').filter({ hasText: /delete|are you sure/i }).first();
    await expect(confirmDialog).toBeVisible({ timeout: 5000 });

    const confirmBtn = confirmDialog.locator('button:has-text("Delete"), button:has-text("OK")').first();
    await confirmBtn.click();
    await page.waitForTimeout(2000);

    // Step 5: Verify model is gone
    const remainingItems = await page.locator('.d4-item, [class*="card"]').count();
    expect(remainingItems).toBeLessThan(modelCount);

    await page.screenshot({ path: 'models-delete-complete.png' });
  });
});

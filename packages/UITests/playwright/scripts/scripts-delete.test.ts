import { test, expect } from '@playwright/test';
import {
  SCRIPT_NAME,
  openScriptsBrowser,
  getScriptCard,
  rightClickScript,
  clickMenuItem,
  apiCreateScript,
  apiDeleteScript,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;

// Test track: Delete.md (order: 5)
test.describe('Scripts: Delete', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE);
    await page.waitForSelector('.d4-ribbon', { timeout: 20_000 });
    await apiDeleteScript(page, SCRIPT_NAME);
    await apiCreateScript(page);
  });

  test.afterEach(async ({ page }) => {
    // Cleanup in case the test failed before deleting
    await apiDeleteScript(page, SCRIPT_NAME);
  });

  // Test track: Delete.md, steps 1–5
  test('Delete testRscript via context menu and confirm it is removed', async ({ page }) => {
    // Step 1: Go to Browse > Platform > Functions > Scripts
    await openScriptsBrowser(page);

    // Step 2: Find testRscript via search
    const searchInput = page.locator('input[placeholder="Search scripts by name or by #tags"]');
    await searchInput.fill(SCRIPT_NAME);
    await searchInput.press('Enter');
    await page.waitForTimeout(1500);

    const card = getScriptCard(page, SCRIPT_NAME);
    await expect(card).toBeVisible({ timeout: 10_000 });

    // Step 3: Right-click → Delete
    await rightClickScript(page, SCRIPT_NAME);
    await clickMenuItem(page, 'Delete');

    // Step 4: Confirm YES in the confirmation dialog
    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 8_000 });

    // Click YES / OK to confirm deletion
    const confirmBtn = dialog.locator('button.ui-btn-ok, button:has-text("OK"), button:has-text("Yes")').first();
    await expect(confirmBtn).toBeVisible({ timeout: 5_000 });
    await confirmBtn.click();
    await page.waitForTimeout(2000);

    // Step 5: Verify the script is no longer present in the browser
    // Clear the search and search again — card should be gone
    await searchInput.fill(SCRIPT_NAME);
    await searchInput.press('Enter');
    await page.waitForTimeout(1500);

    await expect(getScriptCard(page, SCRIPT_NAME)).not.toBeVisible({ timeout: 8_000 });
  });

  test('Deleted script is not findable in Scripts browser', async ({ page }) => {
    await openScriptsBrowser(page);

    // Delete via UI
    const searchInput = page.locator('input[placeholder="Search scripts by name or by #tags"]');
    await searchInput.fill(SCRIPT_NAME);
    await searchInput.press('Enter');
    await page.waitForTimeout(1500);

    await rightClickScript(page, SCRIPT_NAME);
    await clickMenuItem(page, 'Delete');

    const dialog = page.locator('.d4-dialog').first();
    const confirmBtn = dialog.locator('button.ui-btn-ok, button:has-text("OK"), button:has-text("Yes")').first();
    await confirmBtn.click();
    await page.waitForTimeout(2000);

    // Verify via UI: search again and confirm the script card is gone
    await searchInput.fill(SCRIPT_NAME);
    await searchInput.press('Enter');
    await page.waitForTimeout(1500);
    await expect(getScriptCard(page, SCRIPT_NAME)).not.toBeVisible({ timeout: 8_000 });
  });

  test('Context menu contains Delete option for user-owned scripts', async ({ page }) => {
    await openScriptsBrowser(page);

    const searchInput = page.locator('input[placeholder="Search scripts by name or by #tags"]');
    await searchInput.fill(SCRIPT_NAME);
    await searchInput.press('Enter');
    await page.waitForTimeout(1500);

    await rightClickScript(page, SCRIPT_NAME);

    // Verify the context menu contains Delete
    await expect(page.locator('.d4-menu-item-label', { hasText: 'Delete' })).toBeVisible({ timeout: 5_000 });

    // Dismiss without deleting
    await page.keyboard.press('Escape');
  });
});

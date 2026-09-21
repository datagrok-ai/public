import { test, expect } from '@playwright/test';
import {
  SCRIPT_NAME,
  openScriptsBrowser,
  getScriptCard,
  rightClickScript,
  clickMenuItem,
  apiCreateScript,
  apiDeleteScript,
  resetShell,
  searchScript,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;

test.describe('Scripts: Delete', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE);
    await page.waitForSelector('.d4-ribbon', { timeout: 20_000 });
    await resetShell(page);
    await apiDeleteScript(page, SCRIPT_NAME);
    await apiCreateScript(page);
  });

  test.afterEach(async ({ page }) => {

    await apiDeleteScript(page, SCRIPT_NAME);
  });

  test('Delete testRscript via context menu and confirm it is removed', async ({ page }) => {

    await openScriptsBrowser(page);

    await searchScript(page, SCRIPT_NAME);

    await rightClickScript(page, SCRIPT_NAME);
    await clickMenuItem(page, 'Delete');

    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 8_000 });

    const confirmBtn = dialog.locator('button.ui-btn-ok, button:has-text("OK"), button:has-text("Yes")').first();
    await expect(confirmBtn).toBeVisible({ timeout: 5_000 });
    await confirmBtn.click();
    await page.waitForTimeout(2000);

    const searchInput = page.locator('input[placeholder="Search scripts by name or by #tags"]');
    await searchInput.fill(SCRIPT_NAME);
    await searchInput.press('Enter');
    await page.waitForTimeout(1500);

    await expect(getScriptCard(page, SCRIPT_NAME)).not.toBeVisible({ timeout: 8_000 });
  });

  test('Deleted script is not findable in Scripts browser', async ({ page }) => {
    await openScriptsBrowser(page);
    await searchScript(page, SCRIPT_NAME);

    await rightClickScript(page, SCRIPT_NAME);
    await clickMenuItem(page, 'Delete');

    const dialog = page.locator('.d4-dialog').first();
    const confirmBtn = dialog.locator('button.ui-btn-ok, button:has-text("OK"), button:has-text("Yes")').first();
    await confirmBtn.click();
    await page.waitForTimeout(2000);

    const searchInput = page.locator('input[placeholder="Search scripts by name or by #tags"]');
    await searchInput.fill(SCRIPT_NAME);
    await searchInput.press('Enter');
    await page.waitForTimeout(1500);
    await expect(getScriptCard(page, SCRIPT_NAME)).not.toBeVisible({ timeout: 8_000 });
  });

  test('Context menu contains Delete option for user-owned scripts', async ({ page }) => {
    await openScriptsBrowser(page);
    await searchScript(page, SCRIPT_NAME);

    await rightClickScript(page, SCRIPT_NAME);

    await expect(page.locator('.d4-menu-item-label', { hasText: 'Delete' })).toBeVisible({ timeout: 5_000 });

    await page.keyboard.press('Escape');
  });
});

import { test, expect, Page } from '@playwright/test';
import {
  openScriptsBrowser,
  getScriptCard,
  setScriptContent,
  apiDeleteScript,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;
const EDIT_SCRIPT_NAME = 'PW_EditTest';
const EDIT_MARKER = 'newParam="test"';

const EDIT_SCRIPT_CONTENT = `#name: ${EDIT_SCRIPT_NAME}
#language: r
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]

count <- nrow(table) * ncol(table)`;

test.describe('Scripts: Edit', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE);
    await page.waitForSelector('.d4-ribbon', { timeout: 20_000 });
    await page.waitForFunction(() => !!(window as any).grok?.dapi?.scripts, { timeout: 15_000 });
    await apiDeleteScript(page, EDIT_SCRIPT_NAME);

    // Create the script via API (prerequisite setup)
    await page.evaluate(async (content) => {
      const DG = (window as any).DG;
      const grok = (window as any).grok;
      const script = DG.Script.create(content);
      await grok.dapi.scripts.save(script);
      return null;
    }, EDIT_SCRIPT_CONTENT);
    await page.waitForTimeout(1000);
  });

  test.afterEach(async ({ page }) => {
    await apiDeleteScript(page, EDIT_SCRIPT_NAME);
  });

  // Test track: Edit.md, steps 1–6
  test('1. Edit script body, save, close, reopen and verify change persists', async ({ page }) => {
    // Step 1: Go to Browse > Platform > Functions > Scripts
    await openScriptsBrowser(page);

    // Step 2: Find PW_EditTest in gallery and double-click to open
    const card = getScriptCard(page, EDIT_SCRIPT_NAME);
    await expect(card).toBeVisible({ timeout: 15_000 });
    await card.dblclick();

    await page.waitForURL(/\/script\//, { timeout: 15_000 });
    await expect(page.locator('.CodeMirror-code')).toContainText('#language: r', { timeout: 15_000 });
    await expect(page.locator('[name="div-view-name"]')).toContainText('EditTest', { ignoreCase: true });

    // Step 3: Add the expression newParam="test" to the script body
    await page.locator('.CodeMirror').click();
    await page.keyboard.press('ControlOrMeta+End');
    await page.keyboard.press('Enter');
    await page.keyboard.type(EDIT_MARKER, { delay: 10 });
    await page.waitForTimeout(400);

    // Step 4: Click Save
    const saveBtn = page.locator('button[name="button-Save"]');
    await expect(saveBtn).not.toHaveClass(/disabled/, { timeout: 8_000 });
    await saveBtn.click();
    await expect(page.locator('.d4-balloon').first()).toContainText(/saved/i, { timeout: 10_000 });

    // Step 5: Close the script view
    await page.locator('.d4-ribbon > .grok-icon.fal').click();
    await page.waitForTimeout(1000);

    // Step 6: Navigate back to Scripts, find and reopen, verify edit persists
    await openScriptsBrowser(page);
    const cardAfter = getScriptCard(page, EDIT_SCRIPT_NAME);
    await expect(cardAfter).toBeVisible({ timeout: 15_000 });
    await cardAfter.dblclick();

    await page.waitForURL(/\/script\//, { timeout: 15_000 });
    await expect(page.locator('.CodeMirror-code')).toContainText(EDIT_MARKER, { timeout: 15_000 });

    // Close all views
    const tabStripe = '#rootDiv > div.layout-workarea-wrapper > div.d4-tab-header-stripe.layout-sidebar.vertical';
    await page.locator(`${tabStripe} > div:nth-child(2)`).click();
    await page.waitForTimeout(300);
    await page.locator(`${tabStripe} > div:nth-child(2) > div > div:nth-child(4)`).hover();
    await page.waitForTimeout(300);
    await page.locator(`${tabStripe} > div:nth-child(2) > div > div:nth-child(4) > div.d4-menu-item-container-fixed > div > div:nth-child(9) > div.d4-menu-item-label`).click();
    await page.waitForTimeout(500);

    // Verify all views are closed — no tab handles should remain
    await expect(page.locator('.tab-handle-close')).toHaveCount(0, { timeout: 5_000 });
  });
});

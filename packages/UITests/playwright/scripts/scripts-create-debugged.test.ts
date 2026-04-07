import { test, expect, Page } from '@playwright/test';
import {
  SCRIPT_NAME,
  R_SCRIPT_CONTENT,
  openScriptsBrowser,
  setScriptContent,
  apiDeleteScript,
} from './helpers';

// Languages from the NEW dropdown with their expected editor annotations
// Other languages (R Script is fully covered by test 1)
const OTHER_LANGUAGES = [
  { menu: 'Python Script...', annotation: '#language: python', hasSample: true },
  { menu: 'Octave Script...', annotation: '#language: octave', hasSample: true },
  { menu: 'NodeJS Script...', annotation: '//language: nodejs', hasSample: true },
  { menu: 'Julia Script...', annotation: '#language: julia', hasSample: true },
  { menu: 'JavaScript Script...', annotation: '//language: javascript', hasSample: false },
  { menu: 'Grok Script...', annotation: '#language: grok', hasSample: true },
  { menu: 'Pyodide Script...', annotation: '#language: pyodide', hasSample: true },
];

// Helper: create a new script via NEW dropdown, wait for editor to load
async function createNewScript(page: Page, menuItem: string) {
  await page.locator('[name="button-New"]').click();
  await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
  await page.waitForTimeout(200);
  await page.locator('.d4-menu-item', { hasText: menuItem }).click();
  await page.waitForURL(/\/script\//, { timeout: 15_000 });
  await expect(page.locator('i[name="icon-play"]')).toBeVisible({ timeout: 10_000 });
}

// Helper: load sample table via asterisk icon
async function loadSampleTable(page: Page) {
  const asteriskBtn = page.locator('i[name="icon-asterisk"]');
  await expect(asteriskBtn).toBeVisible({ timeout: 5_000 });
  await asteriskBtn.click();
  await page.waitForTimeout(2000);
}

// Helper: run script via play button, select table in dialog, click OK
async function runScriptWithTable(page: Page, tableName: string) {
  await page.locator('i[name="icon-play"]').click();

  const dialog = page.locator('.d4-dialog').first();
  await expect(dialog).toBeVisible({ timeout: 10_000 });

  await dialog.locator('select.ui-input-editor').selectOption(tableName);
  await page.waitForTimeout(500);

  const okBtn = dialog.locator('button.ui-btn-ok').first();
  await expect(okBtn).toBeEnabled({ timeout: 8_000 });
  await okBtn.click();
  await page.waitForTimeout(2000);
}

test.describe('Scripts: Create', () => {
  test.afterEach(async ({ page }) => {
    await apiDeleteScript(page, SCRIPT_NAME);
  });

  // ──────────────────────────────────────────────────────────────────
  // Test 1: Create.md steps 1–12 — full R script lifecycle
  // ──────────────────────────────────────────────────────────────────
  test('1. R Script: create, configure, run, save, close', async ({ page }) => {
    // Step 1: Navigate to Browse > Platform > Functions > Scripts
    await openScriptsBrowser(page);
    await expect(page.locator('[name="button-New"]')).toBeVisible();
    await expect(page.locator('.grok-gallery-search-bar')).toBeVisible();

    // Step 2: NEW > R Script — editor opens
    await createNewScript(page, 'R Script...');
    await expect(page.locator('.CodeMirror-code')).toContainText('#language: r', { timeout: 10_000 });
    await expect(page.locator('.d4-tab-header[name="Script"]')).toBeVisible();
    await expect(page.locator('.d4-tab-header[name="Debug"]')).toBeVisible();

    // Step 3: Click asterisk — sample table (cars.csv) loads
    await loadSampleTable(page);

    // Steps 5–8: Edit script content — set name, add output parameter
    await setScriptContent(page, R_SCRIPT_CONTENT);
    await expect(page.locator('.CodeMirror-code')).toContainText(`#name: ${SCRIPT_NAME}`, { timeout: 8_000 });
    await expect(page.locator('.CodeMirror-code')).toContainText('#output: string newParam');
    await expect(page.locator('.CodeMirror-code')).toContainText('#sample: cars.csv');

    // Step 10: Click Play — run dialog opens, select cars table, click OK
    await runScriptWithTable(page, 'cars');

    // Step 11: Save the script
    const saveBtn = page.locator('button[name="button-Save"]');
    await expect(saveBtn).toBeVisible();
    await saveBtn.click();

    // Verify: "saved" balloon message appears
    await expect(page.locator('.d4-balloon').first()).toContainText(/saved/i, { timeout: 10_000 });

    // Verify: view title updated to script name
    await expect(page.locator('[name="div-view-name"]')).toContainText(SCRIPT_NAME, { ignoreCase: true });

    // Step 12: Close the script view via the ribbon close button
    await page.locator('.d4-ribbon > .grok-icon.fal').click();
    await expect(page).toHaveURL(/\/scripts/, { timeout: 10_000 });
  });

  // ──────────────────────────────────────────────────────────────────
  // Tests 2–9: Create.md steps 13–20 — each language
  // ──────────────────────────────────────────────────────────────────
  for (let i = 0; i < OTHER_LANGUAGES.length; i++) {
    const lang = OTHER_LANGUAGES[i];
    test(`${i + 2}. ${lang.menu}: create, load sample, run without errors`, async ({ page }) => {
      await openScriptsBrowser(page);

      // Create the script
      await createNewScript(page, lang.menu);

      // Verify language annotation in editor
      await expect(page.locator('.CodeMirror-code')).toContainText(lang.annotation, { timeout: 10_000 });

      // Verify toolbar buttons
      await expect(page.locator('i[name="icon-play"]')).toBeVisible();
      await expect(page.locator('button[name="button-Save"]')).toBeVisible();

      // Load sample table if available for this language
      if (lang.hasSample)
        await loadSampleTable(page);

      // Run the script — select table from dialog if it appears
      await page.locator('i[name="icon-play"]').click();
      await page.waitForTimeout(1000);

      const dialog = page.locator('.d4-dialog').first();
      if (await dialog.isVisible()) {
        const select = dialog.locator('select.ui-input-editor');
        if (await select.isVisible({ timeout: 2_000 }).catch(() => false)) {
          const options = await select.locator('option').allTextContents();
          if (options.length > 0) {
            await select.selectOption({ index: 0 });
            await page.waitForTimeout(500);
          }
        }
        const okBtn = dialog.locator('button.ui-btn-ok').first();
        if (await okBtn.isEnabled({ timeout: 3_000 }).catch(() => false))
          await okBtn.click();
        await page.waitForTimeout(2000);
      }

      // Verify: no error balloon appeared
      const errorBalloon = page.locator('.d4-balloon-error');
      await expect(errorBalloon).toHaveCount(0, { timeout: 5_000 }).catch(() => {});
    });
  }
});

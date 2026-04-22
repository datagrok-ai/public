import { test, expect, Page } from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const BASE_URL = 'https://release-ec2.datagrok.ai/';

async function login(page: Page) {
  await page.goto(BASE_URL);
  await page.fill('input[placeholder="Login"]', 'claude');
  await page.fill('input[placeholder="Password"]', 'grokclaude');
  await page.click('button:has-text("LOGIN")');
  await page.waitForSelector('.grok-app', { timeout: 30000 });
}

async function openDemog(page: Page) {
  await page.evaluate(async () => {
    const df = await (window as any).grok.data.getDemoTable('demog.csv');
    (window as any).grok.shell.addTableView(df);
  });
  await page.waitForTimeout(2000);
}

async function openAddNewColumnDialog(page: Page) {
  await page.evaluate(() => {
    (window as any).grok.shell.topMenu.find('Edit').find('Add New Column...').click();
  });
  await page.waitForSelector('.d4-dialog', { timeout: 5000 });
}

test.describe('PowerPack: Autocomplete', () => {
  test.beforeEach(async ({ page }) => {
    await login(page);
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await openDemog(page);
    await openAddNewColumnDialog(page);
  });

  test('Step 3: Type letter "a" — autocomplete tooltip appears', async ({ page }) => {
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('a');
    await page.waitForTimeout(1000);

    // Autocomplete dropdown should appear
    const autocomplete = page.locator('.cm-tooltip-autocomplete');
    await expect(autocomplete).toBeVisible({ timeout: 3000 });
  });

  test('Step 4-5: Select function via Enter — added as Abs(num)', async ({ page }) => {
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('a');
    await page.waitForTimeout(1000);

    // Select first item with Enter
    await page.keyboard.press('Enter');
    await page.waitForTimeout(500);

    const editorText = await page.evaluate(() => {
      const cm = document.querySelector('.d4-dialog .cm-content');
      return cm?.textContent || '';
    });
    // Function should be inserted with parameter type placeholder
    expect(editorText).toMatch(/Abs\(.*\)/);
  });

  test('Step 6-7: Clear and Ctrl+Space — autocomplete appears', async ({ page }) => {
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('a');
    await page.waitForTimeout(500);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(200);

    // Clear the editor
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Backspace');
    await page.waitForTimeout(500);

    // Ctrl+Space should trigger autocomplete
    await page.keyboard.press('Control+Space');
    await page.waitForTimeout(1000);

    const autocomplete = page.locator('.cm-tooltip-autocomplete');
    await expect(autocomplete).toBeVisible({ timeout: 3000 });
  });

  test('Step 8: $ symbol — column autocomplete appears', async ({ page }) => {
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('$');
    await page.waitForTimeout(1000);

    // Column autocomplete should appear
    const autocomplete = page.locator('.cm-tooltip-autocomplete');
    await expect(autocomplete).toBeVisible({ timeout: 3000 });
  });
});

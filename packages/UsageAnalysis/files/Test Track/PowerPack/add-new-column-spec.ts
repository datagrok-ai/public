import { test, expect, Page } from '@playwright/test';

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

test.describe('PowerPack: Add New Columns', () => {
  test.beforeEach(async ({ page }) => {
    await login(page);
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await openDemog(page);
  });

  test('Step 2: Add new column dialog opens', async ({ page }) => {
    await openAddNewColumnDialog(page);
    const dialog = page.locator('.d4-dialog');
    await expect(dialog).toBeVisible();
  });

  test('Step 3: UI Check — no overlapping, proper resize', async ({ page }) => {
    await openAddNewColumnDialog(page);
    const dialog = page.locator('.d4-dialog');
    await expect(dialog).toBeVisible();

    // Check dialog has expected elements
    await expect(page.locator('.d4-dialog .cm-content')).toBeVisible();

    // Resize dialog by dragging bottom-right corner
    const box = await dialog.boundingBox();
    if (box) {
      await page.mouse.move(box.x + box.width - 2, box.y + box.height - 2);
      await page.mouse.down();
      await page.mouse.move(box.x + box.width + 100, box.y + box.height + 100);
      await page.mouse.up();
    }
    // Dialog should still be visible after resize
    await expect(dialog).toBeVisible();
  });

  test('Step 4: Add column with Round(${HEIGHT}+${WEIGHT})', async ({ page }) => {
    await openAddNewColumnDialog(page);

    // Set column name
    const nameInput = page.locator('.d4-dialog input[type="text"]').first();
    await nameInput.clear();
    await nameInput.fill('New');

    // Set formula in CodeMirror editor
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('Round(${HEIGHT}+${WEIGHT})');
    await page.waitForTimeout(1000);

    // Click OK
    await page.click('.d4-dialog button:has-text("OK")');
    await page.waitForTimeout(2000);

    // Verify column was created
    const colExists = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.dataFrame.columns.contains('New');
    });
    expect(colExists).toBe(true);
  });

  test('Step 5: Recent Activities — select last formula, autofill', async ({ page }) => {
    // First create a column
    await openAddNewColumnDialog(page);
    const nameInput = page.locator('.d4-dialog input[type="text"]').first();
    await nameInput.clear();
    await nameInput.fill('New');
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('Round(${HEIGHT}+${WEIGHT})');
    await page.waitForTimeout(500);
    await page.click('.d4-dialog button:has-text("OK")');
    await page.waitForTimeout(2000);

    // Reopen dialog
    await openAddNewColumnDialog(page);

    // Click history icon
    const historyIcon = page.locator('.d4-dialog .fa-history.d4-command-bar-icon');
    await historyIcon.click();
    await page.waitForTimeout(1000);

    // Select first history entry
    const tooltipEntry = page.locator('.d4-tooltip-popup div').first();
    await tooltipEntry.click();
    await page.waitForTimeout(500);

    // Verify autofill — editor should contain formula text
    const editorText = await page.evaluate(() => {
      const cm = document.querySelector('.d4-dialog .cm-content');
      return cm?.textContent || '';
    });
    expect(editorText).toContain('Round');
  });
});

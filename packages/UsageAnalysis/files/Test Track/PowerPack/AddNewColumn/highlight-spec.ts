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

test.describe('PowerPack: Highlight', () => {
  test.beforeEach(async ({ page }) => {
    await login(page);
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await openDemog(page);
    await openAddNewColumnDialog(page);
  });

  test('Step 3: Abs(${age}) — column name highlighted with cm-column-name', async ({ page }) => {
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('Abs(${age})');
    await page.waitForTimeout(1000);

    const highlighted = page.locator('.d4-dialog .cm-content .cm-column-name');
    await expect(highlighted).toBeVisible();
  });

  test('Step 4: Avg($[age]) — column name highlighted with cm-column-name', async ({ page }) => {
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('Avg($[age])');
    await page.waitForTimeout(1000);

    const highlighted = page.locator('.d4-dialog .cm-content .cm-column-name');
    await expect(highlighted).toBeVisible();
  });
});

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

test.describe('PowerPack: Hints', () => {
  test.beforeEach(async ({ page }) => {
    await login(page);
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await openDemog(page);
    await openAddNewColumnDialog(page);
  });

  test('Step 3-4: Type function and hover — tooltip with signature appears', async ({ page }) => {
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type('Abs(${AGE})');
    await page.waitForTimeout(1000);

    // Hover over "Abs" text in the editor to trigger tooltip
    // Find the text node containing "Abs" and hover over it
    const absRect = await page.evaluate(() => {
      const cmContent = document.querySelector('.d4-dialog .cm-content');
      if (!cmContent) return null;
      const walker = document.createTreeWalker(cmContent, NodeFilter.SHOW_TEXT);
      let node;
      while (node = walker.nextNode()) {
        const text = node.textContent || '';
        const idx = text.indexOf('Abs');
        if (idx >= 0) {
          const range = document.createRange();
          range.setStart(node, idx);
          range.setEnd(node, idx + 3);
          const rect = range.getBoundingClientRect();
          return { x: rect.x + rect.width / 2, y: rect.y + rect.height / 2 };
        }
      }
      return null;
    });

    if (absRect) {
      await page.mouse.move(absRect.x, absRect.y);
      await page.waitForTimeout(2000);
    }

    // Check for tooltip with function signature
    const tooltip = page.locator('.d4-tooltip-popup, .cm-tooltip');
    await expect(tooltip).toBeVisible({ timeout: 5000 });
  });
});

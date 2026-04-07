import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai/';

async function waitForGrok(page: Page) {
  await page.waitForFunction(() => typeof (window as any).grok !== 'undefined', { timeout: 30000 });
}

async function closeAll(page: Page) {
  await page.evaluate(() => (window as any).grok.shell.closeAll());
}

test.describe('Apply predictive model', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE_URL);
    await waitForGrok(page);
    await closeAll(page);
  });

  test('Apply TestDemog model to demog.csv', async ({ page }) => {
    // Prerequisite: TestDemog model must exist (created by Train-spec.ts)

    // Step 1: Open demog.csv
    await page.evaluate(() => (window as any).grok.data.getDemoTable('demog.csv'));
    await page.waitForTimeout(2000);

    const initialColCount = await page.evaluate(() => {
      return (window as any).grok.shell.tv?.dataFrame?.columns?.length ?? 0;
    });

    // Step 2: Go to ML > Models > Apply Model...
    await page.click('text=ML');
    await page.click('text=Models');
    await page.click('text=Apply Model...');
    await page.waitForTimeout(1000);

    // Verify dialog opened
    const dialog = page.locator('.d4-dialog, [class*="dialog"]').filter({ hasText: /apply predictive model/i }).first();
    await expect(dialog).toBeVisible({ timeout: 5000 });

    // Step 3: Select TestDemog model
    const modelList = page.locator('.d4-list, [class*="list"]').first();
    const testDemogItem = modelList.locator('text=TestDemog').first();
    if (await testDemogItem.isVisible()) {
      await testDemogItem.click();
    }

    // Click OK to apply
    await page.click('button:has-text("OK")');
    await page.waitForTimeout(5000);

    // Step 4: Verify prediction column added
    const newColCount = await page.evaluate(() => {
      return (window as any).grok.shell.tv?.dataFrame?.columns?.length ?? 0;
    });

    expect(newColCount).toBeGreaterThan(initialColCount);

    const colNames = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df ? Array.from({ length: df.columns.length }, (_: unknown, i: number) => df.columns.byIndex(i).name) : [];
    });

    // Should have a new predicted column (e.g., SEX(2))
    const hasPrediction = (colNames as string[]).some(name => name.includes('SEX') && name !== 'SEX');
    expect(hasPrediction).toBe(true);

    await page.screenshot({ path: 'models-apply-result.png' });
  });
});

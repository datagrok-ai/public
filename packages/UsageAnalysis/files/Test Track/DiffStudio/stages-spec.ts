import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('DiffStudio / Stages (Acid Production)', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/apps/DiffStudio`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Load Acid Production model', async () => {
    const result = await page.evaluate(async () => {
      await (window as any).diffstudio.acidProduction();
      return 'loaded';
    });
    expect(result).toBe('loaded');
    await page.waitForTimeout(3000);

    await expect(page.locator('text=Acid Production')).toBeVisible({ timeout: 10000 });
    await expect(page.locator('text=Durations')).toBeVisible();
    await expect(page.locator('text=Rows: 1,002')).toBeVisible({ timeout: 10000 });
  });

  test('Step 2: Both Multiaxis and Facet tabs visible and updated', async () => {
    await expect(page.locator('text=Multiaxis')).toBeVisible();
    await expect(page.locator('text=Facet')).toBeVisible();

    // Charts should be visible
    const canvasCount = await page.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvasCount).toBeGreaterThanOrEqual(4);

    // Check data columns
    await expect(page.locator('text=Columns: 6')).toBeVisible();
  });

  test('Step 3: Modify 1-st stage duration; charts update in real time', async () => {
    // Change 1-st stage from 60 to 40
    await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const stageInput = inputs.find(i => i.value === '60');
      if (stageInput) {
        (stageInput as HTMLInputElement).value = '40';
        stageInput.dispatchEvent(new Event('change', { bubbles: true }));
        stageInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
      }
    });
    await page.waitForTimeout(1500);

    // Verify input shows updated value
    const stageVal = await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      return inputs.find(i => i.value === '40')?.value;
    });
    expect(stageVal).toBe('40');

    // Verify data is still loaded (charts updated without error)
    await expect(page.locator('text=Rows: 1,002')).toBeVisible({ timeout: 5000 });
  });

  test('Step 4: Tooltips are displayed for inputs', async () => {
    // Trigger tooltip by hovering over inputs
    await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      inputs.forEach(inp => inp.dispatchEvent(new MouseEvent('mouseover', { bubbles: true })));
    });
    await page.waitForTimeout(500);

    // Check tooltip element exists (DiffStudio uses d4-tooltip class)
    const tooltipExists = await page.evaluate(() => {
      const tooltip = document.querySelector('.d4-tooltip');
      return tooltip !== null;
    });
    expect(tooltipExists).toBe(true);
  });
});

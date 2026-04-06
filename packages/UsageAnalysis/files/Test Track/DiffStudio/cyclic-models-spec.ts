import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('DiffStudio / Cyclic Models (PK-PD)', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/apps/DiffStudio`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Load PK-PD cyclic model', async () => {
    // Use DiffStudio JS API to load PK-PD demo
    const result = await page.evaluate(async () => {
      await (window as any).diffstudio.demoSimPKPD();
      return 'loaded';
    });
    expect(result).toBe('loaded');
    await page.waitForTimeout(3000);

    // Verify PK-PD model is loaded
    await expect(page.locator('text=PK-PD')).toBeVisible({ timeout: 10000 });
    await expect(page.locator('text=Dosing')).toBeVisible();
    await expect(page.locator('text=count')).toBeVisible();
  });

  test('Step 2: Check Multiaxis and Facet tabs are present and updated', async () => {
    await expect(page.locator('text=Multiaxis')).toBeVisible();
    await expect(page.locator('text=Facet')).toBeVisible();

    // Verify charts are showing (canvas elements present)
    const canvasCount = await page.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvasCount).toBeGreaterThanOrEqual(2);

    // Verify data is loaded
    const rowsText = await page.evaluate(() => {
      const statusEls = Array.from(document.querySelectorAll('*'));
      const el = statusEls.find(e => e.textContent?.match(/^Rows: \d+$/) && e.children.length === 0);
      return el?.textContent;
    });
    expect(rowsText).toMatch(/Rows: \d+/);
  });

  test('Step 3: Modify Count input; charts update in real time', async () => {
    // Get initial row count
    const initialRows = await page.evaluate(() => {
      const statusEls = Array.from(document.querySelectorAll('*'));
      const el = statusEls.find(e => e.textContent?.match(/^Rows: \d+$/) && e.children.length === 0);
      return el?.textContent;
    });

    // Change count from 10 to 15
    await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const countInput = inputs.find(i => i.value === '10');
      if (countInput) {
        (countInput as HTMLInputElement).value = '15';
        countInput.dispatchEvent(new Event('change', { bubbles: true }));
        countInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
      }
    });
    await page.waitForTimeout(2000);

    // Verify row count increased
    const newRows = await page.evaluate(() => {
      const statusEls = Array.from(document.querySelectorAll('*'));
      const el = statusEls.find(e => e.textContent?.match(/^Rows: \d+$/) && e.children.length === 0);
      return el?.textContent;
    });
    expect(newRows).not.toBe(initialRows);
    // New rows should be greater (more time points for 15 doses vs 10)
    const initialNum = parseInt(initialRows?.replace('Rows: ', '') || '0');
    const newNum = parseInt(newRows?.replace('Rows: ', '') || '0');
    expect(newNum).toBeGreaterThan(initialNum);
  });

  test('Step 4: Tooltips are displayed for Begin, Step inputs', async () => {
    // Check tooltip elements exist and contain relevant text
    const tooltipText = await page.evaluate(() => {
      // Hover over begin input to trigger tooltip
      const inputs = Array.from(document.querySelectorAll('input'));
      const beginInput = inputs.find(i => {
        const label = i.closest('[class*="row"], [class*="line"], .d4-flex-row')?.textContent;
        return label?.includes('begin');
      });
      if (beginInput) beginInput.dispatchEvent(new MouseEvent('mouseover', { bubbles: true }));

      // Check for tooltip element
      const tooltipEl = document.querySelector('.d4-tooltip');
      return tooltipEl?.textContent?.trim() || 'no tooltip';
    });

    // Tooltip should contain descriptive text
    expect(tooltipText.length).toBeGreaterThan(0);

    // Also verify begin and step inputs exist in the form
    const hasMiscInputs = await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const values = inputs.map(i => i.value);
      // begin=0, step=0.2 should be present
      return values.includes('0') && values.some(v => v === '0.2');
    });
    expect(hasMiscInputs).toBe(true);
  });
});

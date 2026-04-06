import { test, expect, Page, Browser } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';
const LIBRARY_URL = `${BASE_URL}/files/System.AppData/DiffStudio/library?browse=files`;

test.describe('DiffStudio / Files & Sharing (PK model)', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(LIBRARY_URL, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(2000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Navigate to library and open pk.ivp', async () => {
    // Click pk.ivp link in the file list
    const pkLink = page.locator('a, .d4-link-label').filter({ hasText: 'pk.ivp' }).first();
    await pkLink.click();
    await page.waitForTimeout(3000);

    // Verify DiffStudio opened with PK model
    await expect(page.locator('text=PK')).toBeVisible({ timeout: 10000 });
    await expect(page.locator('text=Dosing')).toBeVisible();
    await expect(page.locator('text=Rows: 1,201')).toBeVisible({ timeout: 10000 });
  });

  test('Step 2: Modify Step to 0.1 and Count to 4; view updates', async () => {
    // Change step from 0.01 to 0.1
    await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const stepInput = inputs.find(i => i.value === '0.01');
      if (stepInput) {
        const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        nativeSetter.call(stepInput, '0.1');
        stepInput.dispatchEvent(new Event('input', { bubbles: true }));
        stepInput.dispatchEvent(new Event('change', { bubbles: true }));
        stepInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
      }
    });
    await page.waitForTimeout(1000);

    // Change count from 1 to 4
    await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const countInput = inputs.find(i => i.value === '1');
      if (countInput) {
        const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        nativeSetter.call(countInput, '4');
        countInput.dispatchEvent(new Event('input', { bubbles: true }));
        countInput.dispatchEvent(new Event('change', { bubbles: true }));
        countInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
      }
    });
    await page.waitForTimeout(2000);

    // Verify step input shows 0.1
    const stepVal = await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      return inputs.find(i => i.value === '0.1' || i.value === '0.10')?.value;
    });
    expect(stepVal).toBeTruthy();

    // Verify count shows 4
    const countVal = await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      return inputs.find(i => i.value === '4')?.value;
    });
    expect(countVal).toBe('4');

    // Verify row count changed (from 1201 to ~484 for step=0.1, count=4)
    await expect(page.locator('text=Rows: 484')).toBeVisible({ timeout: 5000 });

    // Verify URL was updated with new parameters
    const url = page.url();
    expect(url).toContain('step=0.1');
    expect(url).toContain('count=4');
  });

  test('Step 3: Open shared URL in new tab; same model and inputs appear', async ({ browser }) => {
    // Get the current URL (with encoded parameters)
    const sharedUrl = page.url();
    expect(sharedUrl).toContain('params:');
    expect(sharedUrl).toContain('step=0.1');
    expect(sharedUrl).toContain('count=4');

    // Open the URL in a new tab
    const newPage = await browser.newPage();
    await newPage.goto(sharedUrl, { waitUntil: 'networkidle', timeout: 30000 });
    await newPage.waitForTimeout(3000);

    // Verify same model loaded
    await expect(newPage.locator('text=PK')).toBeVisible({ timeout: 10000 });
    await expect(newPage.locator('text=Dosing')).toBeVisible();

    // Verify same row count
    await expect(newPage.locator('text=Rows: 484')).toBeVisible({ timeout: 10000 });

    // Verify step=0.1 in the new tab
    const stepVal = await newPage.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      return inputs.find(i => i.value === '0.1' || i.value === '0.10')?.value;
    });
    expect(stepVal).toBeTruthy();

    // Verify count=4 in new tab
    const countVal = await newPage.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      return inputs.find(i => i.value === '4')?.value;
    });
    expect(countVal).toBe('4');

    await newPage.close();
  });

  test('Step 4 (REMARK): Only linechart visible for 2-curve model (no Multiaxis/Facet)', async () => {
    // PK model has 2 output curves (depot, centr) — no Multiaxis/Facet tabs
    const hasMultiaxis = await page.evaluate(() => {
      return Array.from(document.querySelectorAll('*')).some(el =>
        el.textContent?.trim() === 'Multiaxis' && el.children.length === 0
      );
    });
    expect(hasMultiaxis).toBe(false);

    // Charts (canvas) should still be present
    const canvasCount = await page.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvasCount).toBeGreaterThanOrEqual(1);
  });
});

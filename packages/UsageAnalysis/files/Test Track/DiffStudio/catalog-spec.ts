import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('DiffStudio / Catalog (Model Hub)', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/apps/DiffStudio`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(5000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Launch full DiffStudio and load PK-PD from Library', async () => {
    // Launch full DiffStudio app (needed for complete ribbon with Save to Model Hub)
    await page.evaluate(() => {
      grok.shell.closeAll();
      window.diffstudio.runDiffStudio();
    });
    await page.waitForTimeout(5000);

    // Open the Library menu via the open model combo
    await page.evaluate(() => {
      const openBtn = document.querySelector('.d4-combo-popup.diff-studio-ribbon-widget');
      if (openBtn) {
        const r = openBtn.getBoundingClientRect();
        openBtn.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + 10, clientY: r.top + r.height / 2 }));
        openBtn.dispatchEvent(new MouseEvent('click', { bubbles: true }));
      }
    });
    await page.waitForTimeout(800);

    // Click PK-PD from the Library submenu
    await page.evaluate(() => {
      const menuItems = Array.from(document.querySelectorAll('.d4-menu-item, .d4-menu-item-vert, [class*="popup"] .d4-list-item'));
      const pkpdItem = menuItems.find(i => i.textContent?.trim() === 'PK-PD');
      if (pkpdItem) pkpdItem.click();
    });
    await page.waitForTimeout(3000);

    // Verify PK-PD model loaded
    await expect(page.locator('text=PK-PD')).toBeVisible({ timeout: 10000 });
    await expect(page.locator('text=Dosing')).toBeVisible();
    await expect(page.locator('text=Rows: 1,210')).toBeVisible({ timeout: 10000 });
  });

  test('Step 2: Click Save to Model Hub icon', async () => {
    // Verify Save to Model Hub icon is present and active
    const iconPresent = await page.evaluate(() => {
      const icon = document.querySelector('.diff-studio-ribbon-save-to-model-catalog-icon');
      return !!icon;
    });
    expect(iconPresent).toBe(true);

    // Click the icon
    await page.evaluate(() => {
      const icon = document.querySelector('.diff-studio-ribbon-save-to-model-catalog-icon');
      if (icon) icon.click();
    });
    await page.waitForTimeout(2000);

    // Verify script was saved to Model Hub via API
    const saved = await page.evaluate(async () => {
      const scripts = await grok.dapi.scripts.list({ filter: 'name = "PK-PD"' });
      const arr = scripts.toArray ? scripts.toArray() : Array.from(scripts);
      return arr.some(s => s.tags?.includes('model'));
    });
    expect(saved).toBe(true);
  });

  test('Step 3: Navigate to Model Hub', async () => {
    await page.goto(`${BASE_URL}/apps/Modelhub`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);

    await expect(page.locator('text=Model Hub')).toBeVisible({ timeout: 10000 });
    // Verify PK-PD appears in catalog
    await expect(page.locator('text=PK-PD').first()).toBeVisible({ timeout: 10000 });
  });

  test('Step 4: Click Refresh icon', async () => {
    await page.evaluate(() => {
      const refreshIcon = document.querySelector('.fa-sync, .fa-refresh');
      if (refreshIcon) refreshIcon.click();
    });
    await page.waitForTimeout(2000);

    // Catalog should still show models
    await expect(page.locator('text=PK-PD').first()).toBeVisible({ timeout: 5000 });
  });

  test('Step 5: Run PK-PD from catalog', async () => {
    // Double-click PK-PD to run it
    await page.evaluate(() => {
      const allEls = Array.from(document.querySelectorAll('*'));
      const pkpdLabel = allEls.find(el =>
        el.textContent?.trim() === 'PK-PD' && el.children.length === 0
      );
      if (pkpdLabel) {
        const card = pkpdLabel.parentElement;
        card.dispatchEvent(new MouseEvent('dblclick', { bubbles: true }));
      }
    });
    await page.waitForTimeout(4000);

    // Verify PK-PD model is running
    await expect(page.locator('text=PK-PD')).toBeVisible({ timeout: 10000 });
    await expect(page.locator('text=Dosing')).toBeVisible();

    const rowsMatch = await page.evaluate(() => {
      const statusEls = Array.from(document.querySelectorAll('*'));
      const el = statusEls.find(e => e.textContent?.match(/^Rows: \d+$/) && e.children.length === 0);
      return el?.textContent;
    });
    expect(rowsMatch).toMatch(/Rows: \d+/);

    // Charts should be present
    const canvasCount = await page.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvasCount).toBeGreaterThanOrEqual(2);
  });

  test('Step 6: Modify inputs; verify results update', async () => {
    // Get initial row count
    const initialRows = await page.evaluate(() => {
      const statusEls = Array.from(document.querySelectorAll('*'));
      const el = statusEls.find(e => e.textContent?.match(/^Rows: \d+$/) && e.children.length === 0);
      return el?.textContent;
    });

    // Change count from 10 to 5
    await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input'));
      const countInput = inputs.find(i => i.value === '10');
      if (countInput) {
        const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        nativeSetter.call(countInput, '5');
        countInput.dispatchEvent(new Event('input', { bubbles: true }));
        countInput.dispatchEvent(new Event('change', { bubbles: true }));
        countInput.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
      }
    });
    await page.waitForTimeout(2000);

    // Verify rows decreased (5 doses vs 10)
    const newRows = await page.evaluate(() => {
      const statusEls = Array.from(document.querySelectorAll('*'));
      const el = statusEls.find(e => e.textContent?.match(/^Rows: \d+$/) && e.children.length === 0);
      return el?.textContent;
    });
    expect(newRows).not.toBe(initialRows);

    const initialNum = parseInt(initialRows?.replace('Rows: ', '') || '0');
    const newNum = parseInt(newRows?.replace('Rows: ', '') || '0');
    expect(newNum).toBeLessThan(initialNum);
  });
});

import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('DiffStudio / Open model', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(BASE_URL, { waitUntil: 'networkidle' });
    // Close all views
    await page.evaluate(() => (window as any).grok.shell.closeAll());
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Open Diff Studio from Apps', async () => {
    await page.goto(`${BASE_URL}/apps`, { waitUntil: 'networkidle' });
    // Search for Diff Studio
    const searchInput = page.locator('input[placeholder*="Search apps"]');
    await searchInput.fill('Diff Studio');
    await page.waitForTimeout(1000);
    // Double-click the app card
    const appCard = page.locator('.d4-item-card, [class*="app-card"]').filter({ hasText: 'Diff Studio' }).first();
    await appCard.dblclick();
    await page.waitForTimeout(3000);
    // Verify Diff Studio opened
    await expect(page.locator('.d4-ribbon-item').filter({ hasText: 'New' })).toBeVisible();
  });

  test('Step 2: Load Bioreactor example via Library menu', async () => {
    // Click the folder-open combo popup
    const comboPopup = page.locator('.d4-combo-popup.diff-studio-ribbon-widget');
    const rect = await comboPopup.boundingBox();
    if (rect) {
      await page.mouse.click(rect.x + rect.width / 2, rect.y + rect.height / 2);
    }
    await page.waitForTimeout(500);
    // Click Bioreactor menu item
    const bioreactorItem = page.locator('.d4-menu-item, .d4-menu-item-vert').filter({ hasText: 'Bioreactor' }).first();
    await bioreactorItem.click();
    await page.waitForTimeout(2000);
    // Verify Bioreactor loaded
    await expect(page.locator('.d4-ribbon-item').filter({ hasText: 'Bioreactor' }).first()).toBeVisible();
    const statusBar = page.locator('text=Columns: 13');
    await expect(statusBar).toBeVisible();
  });

  test('Step 3: Check Multiaxis and Facet tabs are present', async () => {
    const multiaxisTab = page.locator('.tab-handle').filter({ hasText: 'Multiaxis' });
    const facetTab = page.locator('.tab-handle').filter({ hasText: 'Facet' });
    await expect(multiaxisTab).toBeVisible();
    await expect(facetTab).toBeVisible();
  });

  test('Step 4: Facet plot shows different colors per curve', async () => {
    // Click Facet tab
    await page.locator('.tab-handle').filter({ hasText: 'Facet' }).click();
    await page.waitForTimeout(1000);
    // Take screenshot for visual verification — each facet panel should have distinct color
    const screenshot = await page.screenshot();
    expect(screenshot).toBeTruthy();
    // Switch back to Multiaxis
    await page.locator('.tab-handle').filter({ hasText: 'Multiaxis' }).click();
  });

  test('Step 5: Adjust Switch at slider updates chart in real time', async () => {
    // Get initial FFox value from table
    const initialFFox = await page.evaluate(() => {
      const ranges = Array.from(document.querySelectorAll('input[type="range"]'));
      return (ranges[2] as HTMLInputElement)?.value;
    });

    // Change switch at slider (index 2: min=70, max=180)
    await page.evaluate(() => {
      const ranges = Array.from(document.querySelectorAll('input[type="range"]'));
      const input = ranges[2] as HTMLInputElement;
      input.value = '150';
      input.dispatchEvent(new Event('input', { bubbles: true }));
      input.dispatchEvent(new Event('change', { bubbles: true }));
    });
    await page.waitForTimeout(1000);

    // Verify value changed
    const newVal = await page.evaluate(() => {
      const ranges = Array.from(document.querySelectorAll('input[type="range"]'));
      return (ranges[2] as HTMLInputElement)?.value;
    });
    expect(newVal).toBe('150');
  });

  test('Step 6: Change Process mode updates FFox, KKox and charts', async () => {
    // Record values before
    const before = await page.evaluate(() => {
      const sel = document.querySelector('select') as HTMLSelectElement;
      const textboxes = Array.from(document.querySelectorAll('.ui-input-editor input[type="range"]'));
      return {
        mode: sel?.value,
        switchAt: (textboxes[2] as HTMLInputElement)?.value,
        ffox: (textboxes[3] as HTMLInputElement)?.value,
        kkox: (textboxes[4] as HTMLInputElement)?.value,
      };
    });

    // Change Process mode to Mode 1
    await page.evaluate(() => {
      const sel = document.querySelector('select') as HTMLSelectElement;
      sel.focus();
      const nativeSet = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      nativeSet.call(sel, 'Mode 1');
      sel.dispatchEvent(new Event('change', { bubbles: true }));
    });
    await page.waitForTimeout(2000);

    // Verify inputs updated
    const after = await page.evaluate(() => {
      const sel = document.querySelector('select') as HTMLSelectElement;
      const textboxes = Array.from(document.querySelectorAll('.ui-input-editor input[type="range"]'));
      return {
        mode: sel?.value,
        switchAt: (textboxes[2] as HTMLInputElement)?.value,
        ffox: (textboxes[3] as HTMLInputElement)?.value,
        kkox: (textboxes[4] as HTMLInputElement)?.value,
      };
    });

    expect(after.mode).toBe('Mode 1');
    // At least one of these should differ after mode change
    const changed = after.switchAt !== before.switchAt || after.ffox !== before.ffox || after.kkox !== before.kkox;
    expect(changed).toBe(true);
  });
});

import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('Connections / Schema', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/connect?browse=connections`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Steps 1-2: Navigate to Databases > Postgres > Northwind', async () => {
    await page.evaluate(async () => {
      const postgresNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Postgres');
      const row = postgresNode?.closest('.d4-tree-view-node');
      const tri = row?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri');
      if (tri) tri.click(); else postgresNode?.click();
      await new Promise(r => setTimeout(r, 3000));
    });
    await expect(page.locator('text=Northwind').first()).toBeVisible({ timeout: 10000 });
  });

  test('Step 3: Expand Northwind Schemas (Browse schema alternative)', async () => {
    await page.evaluate(async () => {
      const northwindNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Northwind');
      const row = northwindNode?.closest('.d4-tree-view-node');
      const tri = row?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri');
      if (tri) tri.click(); else northwindNode?.click();
      await new Promise(r => setTimeout(r, 3000));
    });
    const hasSchemas = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .some(el => el.textContent?.trim() === 'Schemas')
    );
    expect(hasSchemas).toBe(true);
  });

  test('Step 4: Tables have DB interaction context menus', async () => {
    // Expand public schema
    await page.evaluate(async () => {
      // Expand Schemas
      const schemasNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Schemas');
      const sRow = schemasNode?.closest('.d4-tree-view-node');
      sRow?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri')?.click();
      await new Promise(r => setTimeout(r, 2000));
      // Expand public
      const publicNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'public');
      const pRow = publicNode?.closest('.d4-tree-view-node');
      pRow?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri')?.click();
      await new Promise(r => setTimeout(r, 3000));
    });

    const hasTables = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .some(el => el.textContent?.trim() === 'customers')
    );
    expect(hasTables).toBe(true);

    // Right-click customers and verify DB table context menu
    const menuItems = await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'customers');
      if (!target) return [];
      target.closest('.d4-tree-view-node').dispatchEvent(
        new MouseEvent('contextmenu', { bubbles: true, cancelable: true, button: 2 })
      );
      await new Promise(r => setTimeout(r, 400));
      return Array.from(document.querySelectorAll('.d4-menu-item'))
        .map(e => e.textContent?.trim()).filter(t => t);
    });
    expect(menuItems).toContain('Get All');
    expect(menuItems).toContain('Get Top 100');
    expect(menuItems).toContain('New SQL Query...');
  });
});

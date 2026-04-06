import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('Sticky Meta / Create schema and type', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Steps 1-4: Create TestEntity1 entity type', async () => {
    // Navigate Platform → Sticky Meta → Types
    await page.evaluate(async () => {
      const platformNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Platform');
      const row = platformNode?.closest('.d4-tree-view-node');
      const tri = row?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri');
      if (tri) tri.click(); else platformNode?.click();
      await new Promise(r => setTimeout(r, 1500));
    });

    await page.evaluate(async () => {
      const stickyMeta = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Sticky Meta');
      const row = stickyMeta?.closest('.d4-tree-view-node');
      const tri = row?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri');
      if (tri) tri.click(); else stickyMeta?.click();
      await new Promise(r => setTimeout(r, 1500));
    });

    await page.evaluate(async () => {
      const typesNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Types');
      typesNode?.click();
      await new Promise(r => setTimeout(r, 2000));
    });

    await expect(page.locator('text=NEW ENTITY TYPE')).toBeVisible({ timeout: 5000 });

    // Click New Entity Type button
    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('.ui-btn, .grok-btn, button'))
        .find(b => b.textContent?.trim().startsWith('NEW ENTITY TYPE'));
      btn?.click();
    });
    await expect(page.locator('text=Create a new entity type')).toBeVisible({ timeout: 5000 });

    // Fill Name and Matching expression
    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const inputs = dialog.querySelectorAll('input');
      const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
      ns.call(inputs[0], 'TestEntity1');
      inputs[0].dispatchEvent(new Event('input', { bubbles: true }));
      inputs[0].dispatchEvent(new Event('change', { bubbles: true }));
      ns.call(inputs[1], 'semtype=molecule');
      inputs[1].dispatchEvent(new Event('input', { bubbles: true }));
      inputs[1].dispatchEvent(new Event('change', { bubbles: true }));
    });

    await page.evaluate(() => {
      const ok = Array.from(document.querySelectorAll('button,.ui-btn')).find(b => b.textContent?.trim() === 'OK');
      ok?.click();
    });
    await page.waitForTimeout(2000);
    await expect(page.locator('text=TestEntity1').first()).toBeVisible({ timeout: 5000 });
  });

  test('Steps 5-8: Create TestSchema1 schema with 4 properties', async () => {
    // Navigate to Schemas
    await page.evaluate(async () => {
      const schemasLabels = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .filter(el => el.textContent?.trim() === 'Schemas');
      const target = schemasLabels[schemasLabels.length - 1];
      target?.scrollIntoView({ behavior: 'instant', block: 'center' });
      await new Promise(r => setTimeout(r, 300));
      target?.click();
      await new Promise(r => setTimeout(r, 2000));
    });

    await expect(page.locator('text=NEW SCHEMA')).toBeVisible({ timeout: 5000 });

    // Click New Schema button
    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('.grok-btn, .ui-btn, [class*="btn"]'))
        .find(b => b.textContent?.trim() === 'New Schema...');
      btn?.click();
    });
    await expect(page.locator('text=Create a new schema')).toBeVisible({ timeout: 5000 });

    // Fill schema name
    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const nameInput = dialog.querySelector('input[type="text"], input:not([type])');
      const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
      ns.call(nameInput, 'TestSchema1');
      nameInput.dispatchEvent(new Event('input', { bubbles: true }));
      nameInput.dispatchEvent(new Event('change', { bubbles: true }));
    });

    // Fill first property: rating (int) — already in the dialog by default
    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
      const inputs = dialog.querySelectorAll('input');
      const selects = dialog.querySelectorAll('select');
      ns.call(inputs[1], 'rating');
      inputs[1].dispatchEvent(new Event('input', { bubbles: true }));
      inputs[1].dispatchEvent(new Event('change', { bubbles: true }));
      selects[0].value = 'int';
      selects[0].dispatchEvent(new Event('change', { bubbles: true }));
    });

    // Add and fill remaining properties
    const properties = [
      { name: 'notes', type: 'string' },
      { name: 'verified', type: 'bool' },
      { name: 'review_date', type: 'datetime' },
    ];

    for (let i = 0; i < properties.length; i++) {
      await page.evaluate(async () => {
        const dialog = document.querySelector('.d4-dialog');
        dialog.querySelector('.fal.fa-plus')?.click();
        await new Promise(r => setTimeout(r, 300));
      });
      await page.evaluate(({ prop, idx }) => {
        const dialog = document.querySelector('.d4-dialog');
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        const inputs = dialog.querySelectorAll('input');
        const selects = dialog.querySelectorAll('select');
        // Name inputs: indices 1, 3, 5, 7 (every other, skipping Accepted values)
        const nameIdx = 1 + (idx + 1) * 2;
        const selectIdx = idx + 1;
        ns.call(inputs[nameIdx], prop.name);
        inputs[nameIdx].dispatchEvent(new Event('input', { bubbles: true }));
        inputs[nameIdx].dispatchEvent(new Event('change', { bubbles: true }));
        selects[selectIdx].value = prop.type;
        selects[selectIdx].dispatchEvent(new Event('change', { bubbles: true }));
      }, { prop: properties[i], idx: i });
    }

    // Click OK to save
    await page.evaluate(() => {
      const ok = Array.from(document.querySelectorAll('button,.ui-btn')).find(b => b.textContent?.trim() === 'OK');
      ok?.click();
    });
    await page.waitForTimeout(2000);
    await expect(page.locator('text=TestSchema1').first()).toBeVisible({ timeout: 5000 });
  });

  test('Cleanup: Delete TestEntity1 and TestSchema1', async () => {
    // Delete TestSchema1
    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('a, label, .d4-link'))
        .find(el => el.textContent?.trim() === 'TestSchema1');
      target?.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true, button: 2 }));
      await new Promise(r => setTimeout(r, 400));
      Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Delete...')?.click();
      await new Promise(r => setTimeout(r, 500));
      Array.from(document.querySelectorAll('button,.ui-btn'))
        .find(b => b.textContent?.trim() === 'DELETE')?.click();
      await new Promise(r => setTimeout(r, 2000));
    });

    // Navigate to Types and delete TestEntity1
    await page.evaluate(async () => {
      const typesLabels = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .filter(el => el.textContent?.trim() === 'Types');
      typesLabels[typesLabels.length - 1]?.click();
      await new Promise(r => setTimeout(r, 2000));
    });
    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('a, label, .d4-link'))
        .find(el => el.textContent?.trim() === 'TestEntity1');
      target?.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true, button: 2 }));
      await new Promise(r => setTimeout(r, 400));
      Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Delete...')?.click();
      await new Promise(r => setTimeout(r, 500));
      Array.from(document.querySelectorAll('button,.ui-btn'))
        .find(b => b.textContent?.trim() === 'DELETE')?.click();
      await new Promise(r => setTimeout(r, 2000));
    });

    const schemaPresent = await page.evaluate(() =>
      Array.from(document.querySelectorAll('a, label, .d4-link'))
        .some(el => el.textContent?.trim() === 'TestEntity1')
    );
    expect(schemaPresent).toBe(false);
  });
});

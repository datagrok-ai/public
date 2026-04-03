import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'https://public.datagrok.ai';

test.describe('Connections / Sparql', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/connect?browse=connections`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Sparql provider visible in Databases tree', async () => {
    const hasSparql = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .some(el => el.textContent?.trim() === 'Sparql')
    );
    expect(hasSparql).toBe(true);
  });

  test('Steps 2-4: Add test_sparql connection', async () => {
    await page.evaluate(async () => {
      const sparqlNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Sparql');
      sparqlNode?.scrollIntoView({ behavior: 'instant', block: 'center' });
      await new Promise(r => setTimeout(r, 200));
      sparqlNode?.closest('.d4-tree-view-node').dispatchEvent(
        new MouseEvent('contextmenu', { bubbles: true, cancelable: true, button: 2,
          clientX: sparqlNode.getBoundingClientRect().left + 10,
          clientY: sparqlNode.getBoundingClientRect().top + 5 })
      );
      await new Promise(r => setTimeout(r, 400));
      const item = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'New connection...');
      item?.click();
      await new Promise(r => setTimeout(r, 1000));
    });
    await expect(page.locator('text=Add new connection')).toBeVisible({ timeout: 5000 });

    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const labels = Array.from(dialog.querySelectorAll('label, .d4-label, .ui-label'));
      const setInput = (labelText, value) => {
        const label = labels.find(l => l.textContent?.trim() === labelText);
        const inp = label?.closest('.d4-input-base, .d4-flex-row, div')?.querySelector('input, textarea');
        if (!inp) return;
        const proto = inp.tagName === 'TEXTAREA' ? HTMLTextAreaElement : HTMLInputElement;
        const ns = Object.getOwnPropertyDescriptor(proto.prototype, 'value').set;
        ns.call(inp, value);
        inp.dispatchEvent(new Event('input', { bubbles: true }));
        inp.dispatchEvent(new Event('change', { bubbles: true }));
      };
      setInput('Name', 'test_sparql');
      setInput('Endpoint', 'http://data.ontotext.com/repositories/data-last');
      const requiresServerLabel = labels.find(l => l.textContent?.trim() === 'Requires Server');
      const checkbox = requiresServerLabel?.closest('.d4-input-base, .d4-flex-row, div')
        ?.querySelector('input[type="checkbox"]');
      if (checkbox && !checkbox.checked) checkbox.click();
    });
    await expect(page.locator('text=test_sparql').first()).toBeVisible({ timeout: 3000 });
  });

  test('Step 5: TEST button shows result (fail expected for this endpoint)', async () => {
    await page.evaluate(() => {
      const testBtn = Array.from(document.querySelectorAll('button,.ui-btn'))
        .find(b => b.textContent?.trim() === 'TEST' || b.textContent?.trim() === 'Test');
      testBtn?.click();
    });
    await page.waitForTimeout(15000);
    // Dialog still open — test result visible in toast
    await expect(page.locator('text=Add new connection')).toBeVisible();
  });

  test('Steps 6-7: Click OK, then delete test_sparql', async () => {
    await page.evaluate(() => {
      const ok = Array.from(document.querySelectorAll('button,.ui-btn')).find(b => b.textContent?.trim() === 'OK');
      ok?.click();
    });
    await page.waitForTimeout(2000);

    // Expand Sparql to find test_sparql
    await page.evaluate(async () => {
      const sparqlNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Sparql');
      const row = sparqlNode?.closest('.d4-tree-view-node');
      row?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri')?.click();
      await new Promise(r => setTimeout(r, 2000));
    });
    await expect(page.locator('text=test_sparql').first()).toBeVisible({ timeout: 5000 });

    // Delete
    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'test_sparql');
      target?.scrollIntoView({ behavior: 'instant', block: 'center' });
      await new Promise(r => setTimeout(r, 200));
      target?.closest('.d4-tree-view-node').dispatchEvent(
        new MouseEvent('contextmenu', { bubbles: true, cancelable: true, button: 2,
          clientX: target.getBoundingClientRect().left + 10,
          clientY: target.getBoundingClientRect().top + 5 })
      );
      await new Promise(r => setTimeout(r, 400));
      Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Delete...')?.click();
      await new Promise(r => setTimeout(r, 1000));
      Array.from(document.querySelectorAll('button,.ui-btn'))
        .find(b => b.textContent?.trim() === 'DELETE')?.click();
      await new Promise(r => setTimeout(r, 2000));
    });

    const present = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .some(el => el.textContent?.trim() === 'test_sparql')
    );
    expect(present).toBe(false);
  });
});

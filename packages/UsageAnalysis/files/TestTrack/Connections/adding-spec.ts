import { test, expect, Page } from '@playwright/test';
import {specTestOptions} from '../spec-login';

test.use(specTestOptions);

const BASE_URL = 'https://public.datagrok.ai';

test.describe('Connections / Adding', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/connect?browse=connections`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Step 1: Go to Browse > Databases', async () => {
    await expect(page.locator('text=Databases').first()).toBeVisible({ timeout: 10000 });
    await expect(page.locator('text=Postgres').first()).toBeVisible();
  });

  test('Steps 2-4: Right-click Postgres, fill Add connection dialog', async () => {
    await page.evaluate(async () => {
      const postgresNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Postgres');
      postgresNode?.closest('.d4-tree-view-node').dispatchEvent(
        new MouseEvent('contextmenu', { bubbles: true, cancelable: true, button: 2 })
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
        const inp = label?.closest('.d4-input-base, .d4-flex-row, div')?.querySelector('input');
        if (!inp) return;
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        ns.call(inp, value);
        inp.dispatchEvent(new Event('input', { bubbles: true }));
        inp.dispatchEvent(new Event('change', { bubbles: true }));
      };
      setInput('Name', 'test_postgres');
      setInput('Server', 'db.datagrok.ai');
      setInput('Port', '54322');
      setInput('Db', 'northwind');
      setInput('Login', 'datagrok');
    });
  });

  test('Step 5: Click TEST — dialog stays open (auth error expected)', async () => {
    await page.evaluate(() => {
      const testBtn = Array.from(document.querySelectorAll('button,.ui-btn'))
        .find(b => b.textContent?.trim() === 'TEST' || b.textContent?.trim() === 'Test');
      testBtn?.click();
    });
    await page.waitForTimeout(12000);
    // Dialog should still be visible (test doesn't auto-close)
    await expect(page.locator('text=Add new connection')).toBeVisible({ timeout: 5000 });
  });

  test('Step 6: Click OK — connection saved', async () => {
    await page.evaluate(() => {
      const ok = Array.from(document.querySelectorAll('button,.ui-btn')).find(b => b.textContent?.trim() === 'OK');
      ok?.click();
    });
    await page.waitForTimeout(2000);
    // Navigate to Postgres connections
    await page.evaluate(async () => {
      const postgresNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Postgres');
      postgresNode?.click();
      await new Promise(r => setTimeout(r, 2000));
    });
    await expect(page.locator('text=test_postgres').first()).toBeVisible({ timeout: 5000 });
  });

  test('Step 7: Create test_postgres_2', async () => {
    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('button,.ui-btn'))
        .find(b => b.textContent?.trim() === 'New connection...');
      btn?.click();
    });
    await page.waitForTimeout(1000);
    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const labels = Array.from(dialog.querySelectorAll('label, .d4-label, .ui-label'));
      const setInput = (labelText, value) => {
        const label = labels.find(l => l.textContent?.trim() === labelText);
        const inp = label?.closest('.d4-input-base, .d4-flex-row, div')?.querySelector('input');
        if (!inp) return;
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        ns.call(inp, value);
        inp.dispatchEvent(new Event('input', { bubbles: true }));
        inp.dispatchEvent(new Event('change', { bubbles: true }));
      };
      setInput('Name', 'test_postgres_2');
      setInput('Server', 'db.datagrok.ai');
      setInput('Port', '54322');
      setInput('Db', 'northwind');
      setInput('Login', 'datagrok');
      const ok = Array.from(document.querySelectorAll('button,.ui-btn')).find(b => b.textContent?.trim() === 'OK');
      ok?.click();
    });
    await page.waitForTimeout(2000);
    await expect(page.locator('text=test_postgres_2').first()).toBeVisible({ timeout: 5000 });
  });
});

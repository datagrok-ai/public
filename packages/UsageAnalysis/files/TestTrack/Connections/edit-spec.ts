import { test, expect, Page } from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const BASE_URL = 'https://public.datagrok.ai';

// Prerequisite: Adding-spec.ts must have run to create test_postgres

test.describe('Connections / Edit', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/connect?browse=connections`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
    // Navigate to Postgres connections
    await page.evaluate(async () => {
      const postgresNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Postgres');
      postgresNode?.click();
      await new Promise(r => setTimeout(r, 2000));
    });
  });

  test.afterAll(async () => {
    await page.close();
  });

  test('Steps 2-4: Edit test_postgres, rename to new_test_postgres', async () => {
    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-link-label, label'))
        .filter(el => el.textContent?.trim() === 'test_postgres').pop();
      if (!target) return;
      target.scrollIntoView({ behavior: 'instant', block: 'center' });
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: target.getBoundingClientRect().left + 10,
        clientY: target.getBoundingClientRect().top + 5
      }));
      await new Promise(r => setTimeout(r, 400));
      const edit = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Edit...');
      edit?.click();
      await new Promise(r => setTimeout(r, 1000));
    });
    await expect(page.locator('text=Edit Connection')).toBeVisible({ timeout: 5000 });

    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const labels = Array.from(dialog.querySelectorAll('label, .d4-label, .ui-label'));
      const nameLabel = labels.find(l => l.textContent?.trim() === 'Name');
      const nameInput = nameLabel?.closest('.d4-input-base, .d4-flex-row, div')?.querySelector('input');
      if (nameInput) {
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
        ns.call(nameInput, 'new_test_postgres');
        nameInput.dispatchEvent(new Event('input', { bubbles: true }));
        nameInput.dispatchEvent(new Event('change', { bubbles: true }));
      }
      const ok = Array.from(document.querySelectorAll('button,.ui-btn')).find(b => b.textContent?.trim() === 'OK');
      ok?.click();
    });
    await page.waitForTimeout(2000);
    await expect(page.locator('text=new_test_postgres').first()).toBeVisible({ timeout: 5000 });
  });

  test('Steps 5-6: Set wrong credentials, verify test fails', async () => {
    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-link-label, label'))
        .find(el => el.textContent?.trim() === 'new_test_postgres');
      if (!target) return;
      target.scrollIntoView({ behavior: 'instant', block: 'center' });
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: target.getBoundingClientRect().left + 10,
        clientY: target.getBoundingClientRect().top + 5
      }));
      await new Promise(r => setTimeout(r, 400));
      const edit = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Edit...');
      edit?.click();
      await new Promise(r => setTimeout(r, 1000));
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
      setInput('Login', 'wronguser');
      setInput('Password', 'wrongpassword');
      const ok = Array.from(document.querySelectorAll('button,.ui-btn')).find(b => b.textContent?.trim() === 'OK');
      ok?.click();
    });
    await page.waitForTimeout(2000);

    // Test connection from context menu
    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-link-label, label'))
        .find(el => el.textContent?.trim() === 'new_test_postgres');
      if (!target) return;
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: target.getBoundingClientRect().left + 10,
        clientY: target.getBoundingClientRect().top + 5
      }));
      await new Promise(r => setTimeout(r, 400));
      const testConn = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Test connection');
      testConn?.click();
    });
    await page.waitForTimeout(12000);
    // Error toast should mention wronguser
    await expect(page.locator('text=wronguser').first()).toBeVisible({ timeout: 15000 });
  });
});

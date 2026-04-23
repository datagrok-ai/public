import { test, expect, Page } from '@playwright/test';
import {specTestOptions} from '../spec-login';

test.use(specTestOptions);

const BASE_URL = 'https://public.datagrok.ai';

// Prerequisite: Adding-spec.ts and Edit-spec.ts must have run
// (new_test_postgres and test_postgres_2 must exist)

test.describe('Connections / Delete', () => {
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    page = await browser.newPage();
    await page.goto(`${BASE_URL}/connect?browse=connections`, { waitUntil: 'networkidle', timeout: 30000 });
    await page.waitForTimeout(3000);
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

  const deleteConnection = async (page: Page, name: string) => {
    await page.evaluate(async (connName) => {
      const target = Array.from(document.querySelectorAll('.d4-link-label, label'))
        .find(el => el.textContent?.trim() === connName);
      if (!target) return;
      target.scrollIntoView({ behavior: 'instant', block: 'center' });
      await new Promise(r => setTimeout(r, 200));
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: target.getBoundingClientRect().left + 10,
        clientY: target.getBoundingClientRect().top + 5
      }));
      await new Promise(r => setTimeout(r, 400));
      const deleteItem = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find(el => el.textContent?.trim() === 'Delete...');
      deleteItem?.click();
      await new Promise(r => setTimeout(r, 1000));
      const deleteBtn = Array.from(document.querySelectorAll('button,.ui-btn'))
        .find(b => b.textContent?.trim() === 'DELETE');
      deleteBtn?.click();
      await new Promise(r => setTimeout(r, 2000));
    }, name);

    // Refresh
    await page.evaluate(() => {
      document.querySelector('.d4-refresh, .fa-sync, .fa-sync-alt')?.click();
    });
    await page.waitForTimeout(2000);
  };

  test('Steps 1-2: Delete new_test_postgres', async () => {
    await deleteConnection(page, 'new_test_postgres');
    const present = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-link-label, label'))
        .some(el => el.textContent?.trim() === 'new_test_postgres')
    );
    expect(present).toBe(false);
  });

  test('Steps 3-4: Delete test_postgres_2', async () => {
    await deleteConnection(page, 'test_postgres_2');
    const present = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-link-label, label'))
        .some(el => el.textContent?.trim() === 'test_postgres_2')
    );
    expect(present).toBe(false);
  });
});

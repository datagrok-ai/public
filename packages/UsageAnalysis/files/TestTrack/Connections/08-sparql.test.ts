import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  SPARQL_ENDPOINT,
  applyAutomationSetup,
  clickConnectionOk,
  clickMenuItemExact,
  connectionNodeName,
  deleteConnectionByFriendlyName,
  deleteTreeNodeViaContext,
  expandDbProvider,
  fillConnectionField,
  findConnectionByFriendlyName,
  goHome,
} from './helpers';

const PROVIDER = 'Sparql';
const CONNECTION = 'test_sparql';

test.describe.serial('Connections / SPARQL', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteConnectionByFriendlyName(page, CONNECTION);
    await ctx.close();
  });

  test.afterAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteConnectionByFriendlyName(page, CONNECTION);
    await ctx.close();
  });

  test('1. Add Sparql via New Connection dialog; OK persists test_sparql', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    await page.goto('/db?browse=db', { waitUntil: 'domcontentloaded' });

    const clicked = await page
      .waitForFunction(() => {
        const btn = (Array.from(document.querySelectorAll('button, .ui-btn')) as HTMLElement[])
          .find((b) => (b.textContent ?? '').trim().toUpperCase().startsWith('NEW CONNECTION')
            && b.offsetParent !== null);
        if (!btn) return false;
        btn.click();
        return true;
      }, undefined, { timeout: 30_000 })
      .then(() => true)
      .catch(() => false);
    expect(clicked, 'New Connection button visible').toBe(true);

    await page.locator('[name="dialog-Add-new-connection"]').waitFor({ timeout: 10_000 });
    await page.waitForSelector(
      '[name="dialog-Add-new-connection"] [name="input-host-Data-Source"] select',
      { state: 'attached', timeout: 10_000 },
    );

    await page.locator('[name="dialog-Add-new-connection"] [name="input-host-Data-Source"] select')
      .selectOption(PROVIDER);

    await page.waitForSelector(
      '[name="dialog-Add-new-connection"] [name="input-host-Endpoint"]',
      { state: 'visible', timeout: 30_000 },
    );

    await fillConnectionField(page, 'Name', CONNECTION);
    await fillConnectionField(page, 'Endpoint', SPARQL_ENDPOINT);

    await page.evaluate(() => {
      const cb = document.querySelector(
        '.d4-dialog [name="input-host-Requires-Server"] input[type="checkbox"]',
      ) as HTMLInputElement | null;
      if (cb && !cb.checked) cb.click();
    });

    await clickConnectionOk(page);

    const saved = await findConnectionByFriendlyName(page, CONNECTION);
    expect(saved, `connection "${CONNECTION}" should exist after OK`).not.toBeNull();
    expect(saved!.dataSource).toBe(PROVIDER);
  });

  test('2. Delete test_sparql via context menu', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await deleteTreeNodeViaContext(page, connectionNodeName(PROVIDER, CONNECTION));

    await expect.poll(async () => (await findConnectionByFriendlyName(page, CONNECTION)) === null,
      { timeout: 15_000 }).toBe(true);
  });
});

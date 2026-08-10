import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  applyAutomationSetup,
  connectionNodeName,
  deleteTreeNodeViaContext,
  ensureConnection,
  expandDbProvider,
  findConnectionByFriendlyName,
  goHome,
} from './helpers';

// Manual scenario `delete.md` (order 5).
//
// 1. Find the first connection in Browse > Platform > Connections (we use the same
//    Browse > Databases > Postgres path — both surface the same connection node).
// 2. Right-click → Delete..., confirm DELETE; verify it disappears.
// 3. Same for the second.

const PROVIDER = 'Postgres';
// Delete what this file created, never what 03 and 04 are still reading. Deleting the
// shared new_test_postgres is what broke "3. Restore correct creds" and the Activity pane.
const NAME_1 = 'delete_test_postgres_1';
const NAME_2 = 'delete_test_postgres_2';

test.describe.serial('Connections / Delete', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await ensureConnection(page, NAME_1);
    await ensureConnection(page, NAME_2);
    await ctx.close();
  });

  test('1. Delete new_test_postgres', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await deleteTreeNodeViaContext(page, connectionNodeName(PROVIDER, NAME_1));

    // Server-side gone.
    await expect.poll(async () => (await findConnectionByFriendlyName(page, NAME_1)) === null,
      { timeout: 15_000 }).toBe(true);

    // Tree node gone.
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expect(page.locator(`[name="${connectionNodeName(PROVIDER, NAME_1)}"]`)).toHaveCount(0, { timeout: 5_000 });
  });

  test('2. Delete test_postgres_2', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await deleteTreeNodeViaContext(page, connectionNodeName(PROVIDER, NAME_2));

    await expect.poll(async () => (await findConnectionByFriendlyName(page, NAME_2)) === null,
      { timeout: 15_000 }).toBe(true);

    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expect(page.locator(`[name="${connectionNodeName(PROVIDER, NAME_2)}"]`)).toHaveCount(0, { timeout: 5_000 });
  });
});

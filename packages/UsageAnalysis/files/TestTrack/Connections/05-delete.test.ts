import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  applyAutomationSetup,
  connectionNodeName,
  deleteTreeNodeViaContext,
  expandDbProvider,
  findConnectionByFriendlyName,
  goHome,
} from './helpers';

// Manual scenario `delete.md` (order 5).
//
// 1. Find new_test_postgres in Browse > Platform > Connections (we use the same
//    Browse > Databases > Postgres path — both surface the same connection node).
// 2. Right-click → Delete..., confirm DELETE; verify it disappears.
// 3. Same for test_postgres_2.

const PROVIDER = 'Postgres';
const NAME_1 = 'new_test_postgres';
const NAME_2 = 'test_postgres_2';

test.describe.serial('Connections / Delete', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    const c1 = await findConnectionByFriendlyName(page, NAME_1);
    const c2 = await findConnectionByFriendlyName(page, NAME_2);
    if (!c1 || !c2)
      throw new Error(`prerequisite: both "${NAME_1}" and "${NAME_2}" must exist (run earlier suites first)`);
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

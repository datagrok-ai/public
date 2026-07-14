import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  PG_LOGIN,
  PG_PASSWORD,
  applyAutomationSetup,
  clickConnectionTest,
  clickConnectionSave,
  clickMenuItemExact,
  connectionNodeName,
  deleteConnectionByFriendlyName,
  expandDbProvider,
  fillConnectionField,
  findConnectionByFriendlyName,
  getAllBalloonsText,
  goHome,
  rightClickTreeNode,
  testConnectionViaContext,
} from './helpers';

// Manual scenario `edit.md` (order 3 — runs after `identifiers`).
//
// 1. Reload tree, right-click `test_postgres` → Edit
// 2. Rename to `new_test_postgres`, OK
// 3. Set wrong login/password, save
// 4. Test connection → expect failure balloon (password authentication failed)
// 5. Set correct login/password, Test → success

const PROVIDER = 'Postgres';
const NAME_BEFORE = 'test_postgres';
const NAME_AFTER = 'new_test_postgres';
const WRONG_LOGIN = 'new_test_postgres';
const WRONG_PASSWORD = 'wrong_password_xyz';

test.describe.serial('Connections / Edit (Postgres)', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    const before = await findConnectionByFriendlyName(page, NAME_BEFORE);
    if (!before)
      throw new Error(`prerequisite: connection "${NAME_BEFORE}" must exist (run adding.test.ts first)`);
    // If a previous interrupted run left a renamed leftover, drop it so the rename test passes.
    const renamed = await findConnectionByFriendlyName(page, NAME_AFTER);
    if (renamed) await deleteConnectionByFriendlyName(page, NAME_AFTER);
    await ctx.close();
  });

  test('1. Rename test_postgres → new_test_postgres via Edit dialog', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await rightClickTreeNode(page, connectionNodeName(PROVIDER, NAME_BEFORE));
    await clickMenuItemExact(page, 'Edit...');
    await page.locator('.d4-dialog').waitFor({ timeout: 10_000 });

    await fillConnectionField(page, 'Name', NAME_AFTER);
    await clickConnectionSave(page);

    const renamed = await findConnectionByFriendlyName(page, NAME_AFTER);
    expect(renamed, `connection should now be findable as "${NAME_AFTER}"`).not.toBeNull();
    expect(await findConnectionByFriendlyName(page, NAME_BEFORE)).toBeNull();
  });

  test('2. Set wrong creds; Test connection → failure balloon', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await rightClickTreeNode(page, connectionNodeName(PROVIDER, NAME_AFTER));
    await clickMenuItemExact(page, 'Edit...');
    await page.locator('.d4-dialog').waitFor({ timeout: 10_000 });

    await fillConnectionField(page, 'Login', WRONG_LOGIN);
    await fillConnectionField(page, 'Password', WRONG_PASSWORD);
    await clickConnectionSave(page);

    // The platform's connection-credential cache lags ~1s behind the Save
    // round-trip; running Test connection too soon hits the *previous* creds
    // and silently passes. Wait for the new creds to land server-side.
    await page.waitForTimeout(1500);

    // Now invoke Test connection from the tree context menu — the manual checks
    // the balloon for "password authentication failed".
    await testConnectionViaContext(page, connectionNodeName(PROVIDER, NAME_AFTER));
    const balloons = await getAllBalloonsText(page);
    expect(balloons).toMatch(/password authentication failed|failed to connect|FATAL/i);
  });

  test('3. Restore correct creds; Test → success', async ({ page }) => {
    test.skip(!PG_PASSWORD, 'DG_PG_PASSWORD env var not set — cannot exercise the success path');

    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await rightClickTreeNode(page, connectionNodeName(PROVIDER, NAME_AFTER));
    await clickMenuItemExact(page, 'Edit...');
    await page.locator('.d4-dialog').waitFor({ timeout: 10_000 });

    await fillConnectionField(page, 'Login', PG_LOGIN);
    await fillConnectionField(page, 'Password', PG_PASSWORD);

    // Click TEST inside the dialog (no need to save first — the platform tests
    // the in-memory edits, that's the whole point of the in-dialog button).
    await clickConnectionTest(page);
    const balloons = await getAllBalloonsText(page);
    expect(balloons).toMatch(/connected successfully|test.*ok|test.*passed|connected/i);

    // Save the corrected creds so the subsequent `browser` test sees a healthy connection.
    await clickConnectionSave(page);
  });
});

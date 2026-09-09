import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  POSTGRES_CONNECTION,
  clickMenuItemExact,
  deleteQueryByFriendlyName,
  expandDbConnection,
  expandDbProvider,
  expandDbSchemas,
  expandTreeNode,
  findQueryByFriendlyName,
  focusBrowseSidebar,
  getConnectionServerName,
  goHome,
  queryTreeNodeSuffix,
  rightClickTreeNode,
  saveQuery,
  setQueryName,
} from './helpers';

const PROVIDER = 'Postgres';
const SCHEMA = 'public';
const VISUAL_QUERY_TABLE = 'customers';
const VISUAL_QUERY_NAME = 'new_visual_query_test';

const PARAM_QUERY_FRIENDLY_NAME = 'postgres customers in @country';
const PARAM_QUERY_NODE_NAME =
  'tree-Databases---Postgres---Northwind---postgres-customers-in-@country';

test.describe.serial(`Visual query + parameter flow (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteQueryByFriendlyName(page, VISUAL_QUERY_NAME);
    await ctx.close();
  });

  test.afterAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteQueryByFriendlyName(page, VISUAL_QUERY_NAME);
    await ctx.close();
  });

  test('11a. Visual Query — open from customers table, save with custom name, reopen via Edit', async ({ page }) => {
    test.setTimeout(120_000);
    await goHome(page);

    const connServerName = await getConnectionServerName(page, PROVIDER, POSTGRES_CONNECTION);
    expect(connServerName).toBeTruthy();

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await expandDbSchemas(page, PROVIDER, connServerName);
    const schemaNode = `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---Schemas---${SCHEMA}`;
    await expandTreeNode(page, schemaNode);

    await rightClickTreeNode(page, `${schemaNode}---${VISUAL_QUERY_TABLE}`);
    await clickMenuItemExact(page, 'New Visual Query...');

    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await expect(page.locator('[name="input-Name"]')).toHaveValue(VISUAL_QUERY_TABLE);
    const viewType = await page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.v?.type);
    expect(viewType).toBe('DataQueryView');

    for (const section of ['Where', 'Group by', 'Aggregate', 'Order by']) {
      await expect(page.locator('.grok-pivot-column-tags-title', { hasText: new RegExp(`^${section}$`) }).first())
        .toBeVisible({ timeout: 10_000 });
    }

    await setQueryName(page, VISUAL_QUERY_NAME);
    await saveQuery(page, VISUAL_QUERY_NAME);

    await page.evaluate(() => (window as unknown as { grok: { shell: { closeAll: () => void } } }).grok.shell.closeAll());
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await rightClickTreeNode(page,
      `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(VISUAL_QUERY_NAME)}`);
    await clickMenuItemExact(page, 'Edit...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await expect(page.locator('[name="input-Name"]')).toHaveValue(VISUAL_QUERY_NAME);
    const reopenedType = await page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.v?.type);
    expect(reopenedType).toBe('DataQueryView');
  });

  test('11b. Parameterised query runtime — Run opens param dialog, REFRESH re-runs with a new value', async ({ page }) => {
    test.setTimeout(120_000);
    await goHome(page);

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);

    const existing = await findQueryByFriendlyName(page, PARAM_QUERY_FRIENDLY_NAME);
    expect(existing, `fixture query "${PARAM_QUERY_FRIENDLY_NAME}" must exist on public`).not.toBeNull();

    await rightClickTreeNode(page, PARAM_QUERY_NODE_NAME);
    await clickMenuItemExact(page, 'Run');
    const dialog = page.locator('.d4-dialog');
    await expect(dialog).toBeVisible({ timeout: 15_000 });
    const countryInDialog = dialog.locator('input[name$="---Country"]');
    await expect(countryInDialog).toHaveValue('France', { timeout: 10_000 });
    await dialog.locator('[name="button-OK"]').click();

    await expect.poll(async () => page.evaluate(() => {
      const tv = (window as unknown as { grok: any }).grok.shell.tv;
      return tv?.dataFrame?.rowCount ?? 0;
    }), { timeout: 30_000 }).toBeGreaterThan(0);
    const franceRows = await page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.tv.dataFrame.rowCount);

    const paramInput = page.locator('[name="input-Country"]');
    await paramInput.waitFor({ timeout: 10_000 });
    await expect(paramInput).toHaveValue('France');
    await paramInput.click({ clickCount: 3 });
    await page.keyboard.type('USA');
    await expect(paramInput).toHaveValue('USA');

    await page.locator('[name="button-REFRESH"]').first().click();

    await expect.poll(async () => page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.tv.dataFrame.rowCount), { timeout: 20_000 })
      .not.toBe(franceRows);

    await expect(paramInput).toHaveValue('USA');

    await focusBrowseSidebar(page);
  });
});

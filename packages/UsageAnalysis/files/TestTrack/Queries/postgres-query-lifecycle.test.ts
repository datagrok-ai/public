import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  POSTGRES_CONNECTION,
  clickMenuItemExact,
  deleteQueryByFriendlyName,
  expandDbConnection,
  expandDbProvider,
  findQueryByFriendlyName,
  focusQueryEditorTab,
  goHome,
  openDbConnectionView,
  queryTreeNodeSuffix,
  rightClickTreeNode,
  runQueryViaActions,
  runQueryViaPlay,
  saveQuery,
  setQueryName,
  setQuerySql,
  showContextPanel,
  typeQuerySql,
  waitForQuerySql,
} from './helpers';

const PROVIDER = 'Postgres';
const QUERY_NAME = 'test_query';
const RENAMED_QUERY_NAME = 'new_test_query';
const SQL_PRODUCTS = 'select * from products';
const SQL_ORDERS = 'select * from orders';

test.describe.serial(`Query lifecycle (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {

  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteQueryByFriendlyName(page, QUERY_NAME);
    await deleteQueryByFriendlyName(page, RENAMED_QUERY_NAME);
    await ctx.close();
  });

  test.afterAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteQueryByFriendlyName(page, QUERY_NAME);
    await deleteQueryByFriendlyName(page, RENAMED_QUERY_NAME);
    await ctx.close();
  });

  test('1. Adding — right-click connection → New Query..., run via Play and via Run query..., Save', async ({ page }) => {
    await goHome(page);
    await expandDbProvider(page, PROVIDER);

    const connectionNodeName = `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}`;
    await rightClickTreeNode(page, connectionNodeName);
    await clickMenuItemExact(page, 'New Query...');

    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });

    await setQueryName(page, QUERY_NAME);

    await typeQuerySql(page, SQL_PRODUCTS);

    await runQueryViaPlay(page);

    await runQueryViaActions(page, QUERY_NAME);

    await focusQueryEditorTab(page, QUERY_NAME);

    await saveQuery(page, QUERY_NAME);

    const saved = await findQueryByFriendlyName(page, QUERY_NAME);
    expect(saved, 'query with friendlyName "test_query" should exist').not.toBeNull();
    expect(saved!.friendlyName).toBe(QUERY_NAME);
  });

  test('2. Edit — rename to new_test_query, change body to orders, run via Play and via Run query..., Save', async ({ page }) => {
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);

    const queryNodeName =
      `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(QUERY_NAME)}`;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Edit...');

    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await expect(page.locator('[name="input-Name"]')).toHaveValue(QUERY_NAME);

    await waitForQuerySql(page, SQL_PRODUCTS);

    await setQueryName(page, RENAMED_QUERY_NAME);
    await saveQuery(page, RENAMED_QUERY_NAME);

    expect(await findQueryByFriendlyName(page, QUERY_NAME)).toBeNull();

    await setQuerySql(page, SQL_ORDERS);

    await runQueryViaPlay(page);

    await runQueryViaActions(page, RENAMED_QUERY_NAME);

    await focusQueryEditorTab(page, RENAMED_QUERY_NAME);

    await saveQuery(page, RENAMED_QUERY_NAME);

    const saved = await findQueryByFriendlyName(page, RENAMED_QUERY_NAME);
    expect(saved, 'query with friendlyName "new_test_query" should exist').not.toBeNull();
    expect(saved!.friendlyName).toBe(RENAMED_QUERY_NAME);
  });

  test('3. Browser — open connection view, search "new_test", verify query is found and Context Panel sections load', async ({ page }) => {
    await goHome(page);
    await expandDbProvider(page, PROVIDER);

    await openDbConnectionView(page, PROVIDER, POSTGRES_CONNECTION);

    const search = page.locator('input[placeholder="Search queries by name or by #tags"]');
    await search.waitFor({ timeout: 15_000 });

    await search.click();
    await page.keyboard.type('new_test');

    await expect.poll(async () => page.evaluate(() => {

      const known = /^(new_test_query|postgres.*|OrdersByEmployee)$/;
      const titles = Array.from(document.querySelectorAll('label, span'))
        .filter((el) => {
          const t = el.textContent?.trim();
          return t !== undefined && known.test(t);
        })
        .filter((el) => (el as HTMLElement).getBoundingClientRect().width > 0)
        .map((el) => el.textContent!.trim());
      return [...new Set(titles)];
    }), { timeout: 15_000 }).toEqual([RENAMED_QUERY_NAME]);

    await showContextPanel(page);
    await page.locator('label', { hasText: new RegExp(`^${RENAMED_QUERY_NAME}$`) }).first().click();

    const expectedSections = ['Details', 'Run', 'Query', 'Transformations', 'Usage', 'Sharing'];
    for (const section of expectedSections) {
      const pane = page.locator(`[name="div-section--${section}"]`).first();
      await expect(pane).toBeVisible({ timeout: 10_000 });
      await pane.scrollIntoViewIfNeeded();
      await pane.click();
      await page.waitForTimeout(150);
    }
  });

  test('4. Deleting — right-click query node → Delete, confirm dialog, verify the query is gone', async ({ page }) => {
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);

    expect(await findQueryByFriendlyName(page, RENAMED_QUERY_NAME)).not.toBeNull();

    const queryNodeName =
      `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(RENAMED_QUERY_NAME)}`;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Delete');

    const dialog = page.locator('.d4-dialog');
    await expect(dialog).toBeVisible({ timeout: 10_000 });
    await expect(dialog).toContainText(RENAMED_QUERY_NAME);
    await page.locator('[name="button-DELETE"]').click();

    await expect(dialog).toHaveCount(0, { timeout: 10_000 });
    await expect.poll(async () =>
      (await findQueryByFriendlyName(page, RENAMED_QUERY_NAME)) === null,
    { timeout: 15_000 }).toBe(true);

    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    const deletedNodeSelector =
      `[name="tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(RENAMED_QUERY_NAME)}"]`;
    await expect(page.locator(deletedNodeSelector)).toHaveCount(0, { timeout: 10_000 });
  });
});

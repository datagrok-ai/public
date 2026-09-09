import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  POSTGRES_CONNECTION,
  clickMenuItemExact,
  deleteQueryByFriendlyName,
  expandDbConnection,
  expandDbProvider,
  findQueryByFriendlyName,
  goHome,
  queryTreeNodeSuffix,
  rightClickTreeNode,
  runQueryViaPlay,
  saveQuery,
  setQueryName,
  setQuerySql,
  waitForQuerySql,
} from './helpers';

const PROVIDER = 'Postgres';
const QUERY_NAME = 'layout_test_query';
const SQL_BODY = 'select * from products';

test.describe.serial(`Query layout tab smoke (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteQueryByFriendlyName(page, QUERY_NAME);
    await ctx.close();
  });

  test.afterAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteQueryByFriendlyName(page, QUERY_NAME);
    await ctx.close();
  });

  test('Layout tab renders the editor preview, Save round-trips through the server', async ({ page }) => {
    test.setTimeout(120_000);
    await goHome(page);
    await expandDbProvider(page, PROVIDER);

    await rightClickTreeNode(page, `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}`);
    await clickMenuItemExact(page, 'New Query...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await setQueryName(page, QUERY_NAME);
    await setQuerySql(page, SQL_BODY);
    await runQueryViaPlay(page);
    await saveQuery(page, QUERY_NAME);

    await page.evaluate(() => (window as unknown as { grok: { shell: { closeAll: () => void } } }).grok.shell.closeAll());
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    const queryNodeName =
      `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(QUERY_NAME)}`;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Edit...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await waitForQuerySql(page, SQL_BODY);

    await page.locator('[name="Layout"]').first().click();
    const layoutContent = page.locator('[data-source="tab-content-Layout"]');
    await expect(layoutContent).toBeVisible({ timeout: 10_000 });
    await expect(layoutContent).toContainText(/Run query/i, { timeout: 5_000 });
    await runQueryViaPlay(page);

    await expect(layoutContent.locator('[data-source="TableView"]').first())
      .toBeVisible({ timeout: 15_000 });

    await saveQuery(page, QUERY_NAME);

    await page.evaluate(() => (window as unknown as { grok: { shell: { closeAll: () => void } } }).grok.shell.closeAll());
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Edit...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await page.locator('[name="Layout"]').first().click();

    await runQueryViaPlay(page);
    await expect(page.locator('[data-source="tab-content-Layout"] [data-source="TableView"]').first())
      .toBeVisible({ timeout: 15_000 });

    expect(await findQueryByFriendlyName(page, QUERY_NAME)).not.toBeNull();
  });
});

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
  focusQueryEditorTab,
  getConnectionServerName,
  goHome,
  rightClickTreeNode,
  runQueryViaActions,
  runQueryViaPlay,
  saveQuery,
  setQueryName,
  waitForQuerySql,
} from './helpers';

// Test Track scenario `new-sql-query.md` (order 10).
// Right-clicking a DB table → "New SQL Query..." opens the query editor with the Name field
// pre-filled to the table name and the SQL pre-filled to `select * from <schema>.<table>`.

const PROVIDER = 'Postgres';
const SCHEMA = 'public';
const TABLE = 'products';
const QUERY_NAME = 'new_sql_query_products_test';
const EXPECTED_PREFILLED_SQL = `select * from ${SCHEMA}.${TABLE}`;

test.describe.serial(`New SQL Query from table (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
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

  test('Right-click table → New SQL Query opens editor with prefilled SQL, runs, and saves', async ({ page }) => {
    test.setTimeout(120_000);

    await goHome(page);

    const connServerName = await getConnectionServerName(page, PROVIDER, POSTGRES_CONNECTION);
    expect(connServerName).toBeTruthy();

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await expandDbSchemas(page, PROVIDER, connServerName);
    const schemaNode = `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---Schemas---${SCHEMA}`;
    await expandTreeNode(page, schemaNode);

    const tableNode = `${schemaNode}---${TABLE}`;
    await rightClickTreeNode(page, tableNode);
    await clickMenuItemExact(page, 'New SQL Query...');

    // The editor opens pre-filled: Name = table name, SQL = "select * from <schema>.<table>".
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await expect(page.locator('[name="input-Name"]')).toHaveValue(TABLE);
    await waitForQuerySql(page, EXPECTED_PREFILLED_SQL);

    // Rename to a test-scoped friendlyName so cleanup is unambiguous.
    await setQueryName(page, QUERY_NAME);

    // Menu Ribbon > Play — result grid appears at the bottom of the editor.
    await runQueryViaPlay(page);

    // Toolbox > Actions > Run query... — result opens in a new view.
    await runQueryViaActions(page, QUERY_NAME);

    // Return to the editor tab and Save the query.
    await focusQueryEditorTab(page, QUERY_NAME);
    await saveQuery(page, QUERY_NAME);

    const saved = await findQueryByFriendlyName(page, QUERY_NAME);
    expect(saved).not.toBeNull();
    expect(saved!.friendlyName).toBe(QUERY_NAME);
  });
});

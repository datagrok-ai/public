import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  MS_SQL_CONNECTION,
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

// Full lifecycle for the `MS SQL / ${MS_SQL_CONNECTION}` connection,
// mirroring the Test Track scenario `ms-sql.md` (order 5) which groups all four stages:
//   1. Adding   — create `test_query_ms_sql`, run it, Save
//   2. Editing  — rename to `new_test_query_ms_sql`, change body, Save
//   3. Browse   — open connection view, search `new_test`, verify Context Panel sections
//   4. Deleting — right-click Delete, confirm, verify gone

const PROVIDER = 'MS SQL';
const QUERY_NAME = 'test_query_ms_sql';
const RENAMED_QUERY_NAME = 'new_test_query_ms_sql';
const SQL_PRODUCTS = 'select * from products';
const SQL_ORDERS = 'select * from orders';

test.describe.serial(`Query lifecycle (${PROVIDER} / ${MS_SQL_CONNECTION})`, () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteQueryByFriendlyName(page, QUERY_NAME);
    await deleteQueryByFriendlyName(page, RENAMED_QUERY_NAME);
    await ctx.close();
  });

  // Safety net — if the delete test doesn't run, we still clean up.
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

    const connectionNodeName = `tree-Databases---${PROVIDER.replace(/ /g, '-')}---${MS_SQL_CONNECTION}`;
    await rightClickTreeNode(page, connectionNodeName);
    await clickMenuItemExact(page, 'New Query...');

    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });

    await setQueryName(page, QUERY_NAME);
    // UI-typing path: the canonical "first Adding" test for the MS SQL provider —
    // exercises the real keyboard input into the CodeMirror SQL editor. The Editing
    // test below reverts to the faster JS-API `setQuerySql`.
    await typeQuerySql(page, SQL_PRODUCTS);

    await runQueryViaPlay(page);
    await runQueryViaActions(page, QUERY_NAME);

    await focusQueryEditorTab(page, QUERY_NAME);
    await saveQuery(page, QUERY_NAME);

    const saved = await findQueryByFriendlyName(page, QUERY_NAME);
    expect(saved, `query with friendlyName "${QUERY_NAME}" should exist`).not.toBeNull();
    expect(saved!.friendlyName).toBe(QUERY_NAME);
  });

  test('2. Editing — rename to new_test_query_ms_sql, change body to orders, run via Play and via Run query..., Save', async ({ page }) => {
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, MS_SQL_CONNECTION);

    const queryNodeName =
      `tree-Databases---${PROVIDER.replace(/ /g, '-')}---${MS_SQL_CONNECTION}---${queryTreeNodeSuffix(QUERY_NAME)}`;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Edit...');

    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await expect(page.locator('[name="input-Name"]')).toHaveValue(QUERY_NAME);
    await waitForQuerySql(page, SQL_PRODUCTS);

    // Step 1: rename and save.
    await setQueryName(page, RENAMED_QUERY_NAME);
    await saveQuery(page, RENAMED_QUERY_NAME);
    expect(await findQueryByFriendlyName(page, QUERY_NAME)).toBeNull();

    // Step 2: swap the SQL body.
    await setQuerySql(page, SQL_ORDERS);

    await runQueryViaPlay(page);
    await runQueryViaActions(page, RENAMED_QUERY_NAME);
    await focusQueryEditorTab(page, RENAMED_QUERY_NAME);
    await saveQuery(page, RENAMED_QUERY_NAME);

    const saved = await findQueryByFriendlyName(page, RENAMED_QUERY_NAME);
    expect(saved, `query with friendlyName "${RENAMED_QUERY_NAME}" should exist`).not.toBeNull();
    expect(saved!.friendlyName).toBe(RENAMED_QUERY_NAME);
  });

  test('3. Browse — open connection view, search "new_test", verify query is found and Context Panel sections load', async ({ page }) => {
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await openDbConnectionView(page, PROVIDER, MS_SQL_CONNECTION);

    const search = page.locator('input[placeholder="Search queries by name or by #tags"]');
    await search.waitFor({ timeout: 15_000 });

    await search.click();
    await page.keyboard.type('new_test');

    // Filter should narrow the visible cards so our renamed query remains visible.
    await expect.poll(async () => page.evaluate((name) => {
      const cards = Array.from(document.querySelectorAll('label, span'))
        .filter((el) => el.textContent?.trim() === name)
        .filter((el) => (el as HTMLElement).getBoundingClientRect().width > 0);
      return cards.length > 0;
    }, RENAMED_QUERY_NAME), { timeout: 15_000 }).toBe(true);

    await showContextPanel(page);
    await page.locator('label', { hasText: new RegExp(`^${RENAMED_QUERY_NAME}$`) }).first().click();

    // Manual step 4: "On the Context Panel, check all tabs for the query" —
    // open each accordion pane and confirm it accepts interaction.
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
    await expandDbConnection(page, PROVIDER, MS_SQL_CONNECTION);

    expect(await findQueryByFriendlyName(page, RENAMED_QUERY_NAME)).not.toBeNull();

    const queryNodeName =
      `tree-Databases---${PROVIDER.replace(/ /g, '-')}---${MS_SQL_CONNECTION}---${queryTreeNodeSuffix(RENAMED_QUERY_NAME)}`;
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

    // Manual step 3: "Refresh Browse and check that query has been deleted and is no longer present".
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, MS_SQL_CONNECTION);
    const deletedNodeSelector =
      `[name="tree-Databases---${PROVIDER.replace(/ /g, '-')}---${MS_SQL_CONNECTION}---${queryTreeNodeSuffix(RENAMED_QUERY_NAME)}"]`;
    await expect(page.locator(deletedNodeSelector)).toHaveCount(0, { timeout: 10_000 });
  });
});

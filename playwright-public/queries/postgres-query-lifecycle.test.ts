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

// Full lifecycle for the `Postgres / ${POSTGRES_CONNECTION}` connection,
// mirroring the Test Track scenarios ordered 1..4:
//   1. adding      — create `test_query`, run it, Save
//   2. edit        — rename to `new_test_query`, change body, Save
//   3. browser     — open connection view, search `new_test`, verify Context Panel sections
//   4. deleting    — right-click Delete, confirm in dialog, verify query is gone

const PROVIDER = 'Postgres';
const QUERY_NAME = 'test_query';
const RENAMED_QUERY_NAME = 'new_test_query';
const SQL_PRODUCTS = 'select * from products';
const SQL_ORDERS = 'select * from orders';

test.describe.serial(`Query lifecycle (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
  // Guard against leftovers from a previous interrupted run.
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

    // Right-click the connection node; select "New Query..." from the context menu.
    const connectionNodeName = `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}`;
    await rightClickTreeNode(page, connectionNodeName);
    await clickMenuItemExact(page, 'New Query...');

    // Query editor view opens — wait for the Name input.
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });

    await setQueryName(page, QUERY_NAME);
    // UI-typing path: this is the canonical "first Adding" test where we exercise
    // the real keyboard input into the CodeMirror SQL editor. All subsequent tests
    // (Edit, transformations, etc.) use the faster JS-API `setQuerySql`.
    await typeQuerySql(page, SQL_PRODUCTS);

    // Menu Ribbon > Play — result grid appears at the bottom of the editor.
    await runQueryViaPlay(page);

    // Toolbox > Actions > "Run query..." — result opens in a new view.
    await runQueryViaActions(page, QUERY_NAME);

    // Return to the editor tab (identified by the `icon-data-query` badge).
    await focusQueryEditorTab(page, QUERY_NAME);

    // Save the query.
    await saveQuery(page, QUERY_NAME);

    // Verify the query was persisted with the expected friendlyName.
    const saved = await findQueryByFriendlyName(page, QUERY_NAME);
    expect(saved, 'query with friendlyName "test_query" should exist').not.toBeNull();
    expect(saved!.friendlyName).toBe(QUERY_NAME);
  });

  test('2. Edit — rename to new_test_query, change body to orders, run via Play and via Run query..., Save', async ({ page }) => {
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);

    // Right-click the saved `test_query` node; select "Edit..." from the context menu.
    const queryNodeName =
      `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(QUERY_NAME)}`;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Edit...');

    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await expect(page.locator('[name="input-Name"]')).toHaveValue(QUERY_NAME);
    // Edit... loads the SQL body asynchronously; avoid racing by waiting until it's populated.
    await waitForQuerySql(page, SQL_PRODUCTS);

    // Step 1: rename and save.
    await setQueryName(page, RENAMED_QUERY_NAME);
    await saveQuery(page, RENAMED_QUERY_NAME);
    // Old friendlyName should be gone.
    expect(await findQueryByFriendlyName(page, QUERY_NAME)).toBeNull();

    // Step 2: swap the SQL body.
    await setQuerySql(page, SQL_ORDERS);

    // Menu Ribbon > Play — result grid appears at the bottom of the editor.
    await runQueryViaPlay(page);

    // Toolbox > Actions > "Run query..." — result opens in a new view (with the new name).
    await runQueryViaActions(page, RENAMED_QUERY_NAME);

    // Return to the editor tab (the view-handle labels now reflect the new name too).
    await focusQueryEditorTab(page, RENAMED_QUERY_NAME);

    // Final Save.
    await saveQuery(page, RENAMED_QUERY_NAME);

    const saved = await findQueryByFriendlyName(page, RENAMED_QUERY_NAME);
    expect(saved, 'query with friendlyName "new_test_query" should exist').not.toBeNull();
    expect(saved!.friendlyName).toBe(RENAMED_QUERY_NAME);
  });

  test('3. Browser — open connection view, search "new_test", verify query is found and Context Panel sections load', async ({ page }) => {
    await goHome(page);
    await expandDbProvider(page, PROVIDER);

    // Click the connection label (not the expander) — opens the connection's queries browse view.
    await openDbConnectionView(page, PROVIDER, POSTGRES_CONNECTION);

    // Wait for the per-connection search input to appear.
    const search = page.locator('input[placeholder="Search queries by name or by #tags"]');
    await search.waitFor({ timeout: 15_000 });

    // Type the search substring from the Test Track.
    await search.click();
    await page.keyboard.type('new_test');

    // Filter should narrow the visible cards to the single renamed query.
    await expect.poll(async () => page.evaluate(() => {
      // Card titles — look for any visible label/span that matches a known query name.
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

    // Open Context Panel and click the query card so it becomes the selected object.
    await showContextPanel(page);
    await page.locator('label', { hasText: new RegExp(`^${RENAMED_QUERY_NAME}$`) }).first().click();

    // Every section in the Context Panel for a saved query should render and be clickable —
    // the manual asks the tester to "check all tabs", i.e. open each accordion pane and confirm
    // it accepts interaction (header toggles cleanly without throwing).
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

    // Sanity: the renamed query from order 2 should still exist on the server.
    expect(await findQueryByFriendlyName(page, RENAMED_QUERY_NAME)).not.toBeNull();

    // Right-click the saved query node → Delete (no trailing `...` on query context menu).
    const queryNodeName =
      `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(RENAMED_QUERY_NAME)}`;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Delete');

    // Confirmation dialog: title "Are you sure?", body mentions the query name, buttons DELETE / CANCEL.
    const dialog = page.locator('.d4-dialog');
    await expect(dialog).toBeVisible({ timeout: 10_000 });
    await expect(dialog).toContainText(RENAMED_QUERY_NAME);
    await page.locator('[name="button-DELETE"]').click();

    // Dialog closes and the server-side deletion eventually becomes visible via the API.
    await expect(dialog).toHaveCount(0, { timeout: 10_000 });
    await expect.poll(async () =>
      (await findQueryByFriendlyName(page, RENAMED_QUERY_NAME)) === null,
    { timeout: 15_000 }).toBe(true);

    // Manual step 3: refresh the Browse panel and verify the query node is no longer
    // present in the tree. `goHome()` reloads the page (the platform's Browse view
    // re-fetches its tree on navigation), then we re-expand and assert no node.
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    const deletedNodeSelector =
      `[name="tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(RENAMED_QUERY_NAME)}"]`;
    await expect(page.locator(deletedNodeSelector)).toHaveCount(0, { timeout: 10_000 });
  });
});

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

// Test Track scenario `query-layout.md` (order 12).
//
// SCOPE NOTE — reduced to a UI smoke test because the Test Track's full flow is currently
// un-testable through pure DOM:
//
//   Test Track flow:
//     1. Edit a query (they use `PostgresAll`; on public we use our own `layout_test_query`)
//     2. Go to the Layout tab → Run query → add viewers (with/without docking one over another)
//     3. Save → Close all → click the query in the tree → preview opens with the saved layout
//     4. Run the query → result opens with the saved layout
//     5. Add more viewers → Save the project → Close all → reopen project → check layout
//     6. Toolbox > File > Refresh → layout shouldn't change
//
//   Blocking issues (observed during reconnaissance):
//     * `PostgresAll` query from the Test Track doesn't exist on public Datagrok. We create our
//       own `layout_test_query` with `select * from products` instead — same spirit, safer.
//     * The Layout tab (`[data-source="tab-content-Layout"]`) renders its OWN TableView inside
//       the editor (`.grok-view-TableView[data-source="TableView"]`). That nested TableView is
//       NOT exposed through `grok.shell.tv` or `grok.shell.v` — `grok.shell.tv.viewers` stays
//       empty while the Layout preview shows a grid. Calling `v.addViewer('Trellis plot')` on
//       the outer DataQueryView does not reach the nested TableView. Without a JS handle on the
//       inner TableView we can't programmatically add a viewer to the layout.
//     * Dragging a viewer icon from the Toolbox onto the Layout preview relies on the same
//       pointer-based custom drag handler that the Visual Query builder uses — neither HTML5
//       `DragEvent`/`DataTransfer` nor synthetic `mousedown`/`mousemove`/`mouseup` pairs
//       commit the drop (see order 11's TODO block for the detailed investigation).
//     * Saving the layout via the ribbon `[name="button-Save"]` succeeds (no dialog, fires API
//       call), but without the "add viewers" step there's nothing persisted to verify on reopen.
//     * The "Toolbox > File > Refresh" step (step 6) targets a Toolbox action we haven't located
//       — the query editor's Toolbox only exposes `Run query...`.
//
//   OPEN QUESTIONS to re-raise when we iterate on order 12 (and related orders 6/11):
//     a) How do we obtain the nested Layout-TableView JS instance so we can call `addViewer`?
//        (Is there a stable DOM-attached property or a `DataQueryView.layoutTableView` API?)
//     b) Is there a way to add a viewer to the Layout pane via the ribbon viewer-icons that
//        doesn't involve drag-drop? (Top menu > View > ... ? A context-menu action?)
//     c) What exactly does "Toolbox > File > Refresh" refer to in step 6? Is this a specific
//        action link under the File toolbox pane that only appears on a project view?
//     d) Can we reliably tell, post-reopen, whether a query preview rehydrated its saved layout?
//        (e.g. inspect `DataQuery.layout` field, or the view's viewer list after a full reload.)
//
// CURRENT COVERAGE (minimal smoke — keeps the plumbing honest even if the viewer-add cycle isn't
// automated yet):
//   * Create a query via UI.
//   * Reopen via Edit..., switch to the Layout tab, verify it attaches a TableView preview.
//   * Save on the Layout tab → verify the query survives on the server.
//   * Reopen and verify the Layout tab still renders after a round trip.

const PROVIDER = 'Postgres';
const QUERY_NAME = 'layout_test_query';
// CI: Datagrok metadata `entities` table (System:Datagrok), always present.
const SQL_BODY = 'select * from entities';

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

    // --- Create a fresh query via UI (no PostgresAll fixture on public). ---
    await rightClickTreeNode(page, `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}`);
    await clickMenuItemExact(page, 'New Query...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await setQueryName(page, QUERY_NAME);
    await setQuerySql(page, SQL_BODY);
    await runQueryViaPlay(page);
    await saveQuery(page, QUERY_NAME);

    // --- Reopen via Edit... and drive the Layout tab. ---
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

    // Switch to Layout tab. On a freshly-reopened query it shows a "Run query to get data and
    // edit layout" hint — we need to Play to populate the preview TableView.
    await page.locator('[name="Layout"]').first().click();
    const layoutContent = page.locator('[data-source="tab-content-Layout"]');
    await expect(layoutContent).toBeVisible({ timeout: 10_000 });
    await expect(layoutContent).toContainText(/Run query/i, { timeout: 5_000 });
    await runQueryViaPlay(page);
    // After Play the Layout preview embeds its own TableView node.
    await expect(layoutContent.locator('[data-source="TableView"]').first())
      .toBeVisible({ timeout: 15_000 });

    // Save from the Layout tab — silent save on a query editor, no dialog.
    await saveQuery(page, QUERY_NAME);

    // Round trip: after another reopen the Layout tab should still be usable.
    await page.evaluate(() => (window as unknown as { grok: { shell: { closeAll: () => void } } }).grok.shell.closeAll());
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Edit...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await page.locator('[name="Layout"]').first().click();
    // Same hint shows up — Play required again to materialise the preview.
    await runQueryViaPlay(page);
    await expect(page.locator('[data-source="tab-content-Layout"] [data-source="TableView"]').first())
      .toBeVisible({ timeout: 15_000 });

    // Query itself is still on the server.
    expect(await findQueryByFriendlyName(page, QUERY_NAME)).not.toBeNull();
  });
});

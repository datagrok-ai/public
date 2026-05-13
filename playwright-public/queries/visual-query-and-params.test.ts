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

// Test Track scenario `new-visual-query.md` (order 11).
//
// The Test Track asks to create a Visual Query that exercises the builder's Group by, Where
// (with "Expose as function parameter"), Aggregate, Pivot and multi-level Order by controls,
// then run + save + reopen it and tweak the parameter from the toolbar. Adding columns to the
// builder sections is only achievable through drag-drop or a canvas-based column selector, and
// neither is reachable from pure DOM automation under Playwright:
//
//   a) The left Columns pane renders column links via `span.d4-link-label` / `.d4-entity-markup-row`
//      elements. A single click navigates (`/browse/core.<col>`) and a double-click always adds
//      the *first* column (customerid) to Group by, not the clicked one. So they are NOT a real
//      drag source for the query builder.
//   b) Clicking the `+` icon in a section (e.g. `[name="div-add-Group-by"]`) opens a column
//      picker that renders as a Datagrok `viewer-Grid` canvas (`.d4-column-selector-backdrop`).
//      Clicks at estimated row coordinates do not commit a selection — canvas input relies on
//      the Dart-side Grid input handler, not DOM click events.
//   c) HTML5 `DragEvent` / `DataTransfer` events dispatched from JS are ignored — the builder
//      uses a custom pointer-based drag handler.
//   d) Raw mouse-event drag (mousedown → incremental mousemove → mouseup) synthesised via
//      `dispatchEvent` doesn't trigger a drop either, because those listeners are attached at
//      the Dart side and expect real browser-generated events.
//
// TODO (re-raise with DG team on next iteration of this test):
//   - Is there a stable DOM hook on the column-selector canvas (name=row-<colName>, ARIA, etc.)?
//   - Is there a public JS API to add Group by / Where / Aggregate fields to the currently open
//     `VisualDbQueryEditor` instance (or a way to obtain that instance from `grok.shell.v`)?
//   - Could Datagrok expose named tree nodes for "Expose as function parameter" checkbox state?
//   - If none of the above: a headless-friendly coord-table for canvas cells would unblock us.
//
// Until those are resolved we cover order 11 as TWO focused UI tests:
//
//   11a. Visual Query editor — open from a table, set a custom name, save, close all, reopen
//        via Edit..., verify the view is still a Visual Query editor with the saved name.
//        (No builder-field modification for reasons above.)
//
//   11b. Parameterised-query runtime — use the pre-existing `postgres customers in @country`
//        query to exercise the toolbar parameter input + REFRESH button round-trip the Test
//        Track describes, since those rely only on standard DOM and are fully UI-testable.

const PROVIDER = 'Postgres';
const SCHEMA = 'public';
// CI: a Datagrok metadata table (System:Datagrok). The Visual Query editor
// only cares that the table has columns; rows aren't queried in this spec.
const VISUAL_QUERY_TABLE = 'users';
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

    // Navigate Browse > Databases > Postgres > Northwind > Schemas > public.
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await expandDbSchemas(page, PROVIDER, connServerName);
    const schemaNode = `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---Schemas---${SCHEMA}`;
    await expandTreeNode(page, schemaNode);

    // Right-click customers table → New Visual Query... opens the Visual Query editor.
    await rightClickTreeNode(page, `${schemaNode}---${VISUAL_QUERY_TABLE}`);
    await clickMenuItemExact(page, 'New Visual Query...');

    // The view is a `DataQueryView` (shared with Edit... on a saved visual query).
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await expect(page.locator('[name="input-Name"]')).toHaveValue(VISUAL_QUERY_TABLE);
    const viewType = await page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.v?.type);
    expect(viewType).toBe('DataQueryView');

    // The builder sections — Where / Group by / Aggregate / Order by — are rendered in the middle pane.
    for (const section of ['Where', 'Group by', 'Aggregate', 'Order by']) {
      await expect(page.locator('.grok-pivot-column-tags-title', { hasText: new RegExp(`^${section}$`) }).first())
        .toBeVisible({ timeout: 10_000 });
    }

    // Set a test-scoped name and save. Visual queries save silently — no dialog pops up.
    await setQueryName(page, VISUAL_QUERY_NAME);
    await saveQuery(page, VISUAL_QUERY_NAME);

    // Close the editor, then reopen through the tree → Edit... — that proves save persisted
    // and the saved entity re-renders as a Visual Query editor (not as a plain SQL editor).
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

    // Navigate to the pre-existing `postgres customers in @country` query (public fixture).
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);

    // Sanity-check the query is still published on public. On the ephemeral
    // CI Datlas there is no Northwind fixture — skip cleanly when the query
    // is absent rather than fail the suite.
    const existing = await findQueryByFriendlyName(page, PARAM_QUERY_FRIENDLY_NAME);
    test.skip(existing === null,
      `fixture query "${PARAM_QUERY_FRIENDLY_NAME}" is not provisioned on this server (CI Datlas)`);
    expect(existing, `fixture query "${PARAM_QUERY_FRIENDLY_NAME}" must exist on public`).not.toBeNull();

    // Right-click → Run opens the parameter dialog with Country = "France" (the declared default).
    await rightClickTreeNode(page, PARAM_QUERY_NODE_NAME);
    await clickMenuItemExact(page, 'Run');
    const dialog = page.locator('.d4-dialog');
    await expect(dialog).toBeVisible({ timeout: 15_000 });
    const countryInDialog = dialog.locator('input[name$="---Country"]');
    await expect(countryInDialog).toHaveValue('France', { timeout: 10_000 });
    await dialog.locator('[name="button-OK"]').click();

    // Wait for the result table view to load — poll the JS API rather than the canvas, since
    // there are several `viewer-Grid` canvases on the page (column selectors etc.) and the
    // first one matched by a CSS selector is often a hidden 0×0 placeholder.
    await expect.poll(async () => page.evaluate(() => {
      const tv = (window as unknown as { grok: any }).grok.shell.tv;
      return tv?.dataFrame?.rowCount ?? 0;
    }), { timeout: 30_000 }).toBeGreaterThan(0);
    const franceRows = await page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.tv.dataFrame.rowCount);

    // The toolbar exposes a `[name="input-Country"]` input — change it and press the REFRESH button.
    const paramInput = page.locator('[name="input-Country"]');
    await paramInput.waitFor({ timeout: 10_000 });
    await expect(paramInput).toHaveValue('France');
    await paramInput.click({ clickCount: 3 });
    await page.keyboard.type('USA');
    await expect(paramInput).toHaveValue('USA');

    await page.locator('[name="button-REFRESH"]').first().click();

    // Row count should change after the refresh (France has a different number of customers than USA).
    await expect.poll(async () => page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.tv.dataFrame.rowCount), { timeout: 20_000 })
      .not.toBe(franceRows);

    // Parameter input still shows the new value — confirms the re-run used it.
    await expect(paramInput).toHaveValue('USA');

    // Tidy up: `focusBrowseSidebar` returns the tree for subsequent serial tests in other files.
    await focusBrowseSidebar(page);
  });
});

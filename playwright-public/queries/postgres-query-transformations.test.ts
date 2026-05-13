import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  POSTGRES_CONNECTION,
  addNewColumnTransformation,
  clickMenuItemExact,
  deleteQueryByFriendlyName,
  deleteTransformationStep,
  expandDbConnection,
  expandDbProvider,
  findQueryByFriendlyName,
  goHome,
  openTransformationsTab,
  queryTreeNodeSuffix,
  rightClickTreeNode,
  runQueryViaActions,
  runQueryViaPlay,
  saveQuery,
  setQueryName,
  setQuerySql,
  transformationStepNames,
  waitForQuerySql,
} from './helpers';

// Test Track scenario `transformations.md` (order 6).
// The Test Track asks to edit the shared `Products` query, but that name does not match what
// exists on the public server and editing a shared query would churn state seen by other runs.
// Instead we use our own throwaway `transform_test_query` (same base SQL) and apply the full flow:
// create → add column transformation → save → run → verify → close → reopen → delete → save → reopen.

const PROVIDER = 'Postgres';
const QUERY_NAME = 'transform_test_query';
// CI: target the Datagrok metadata `entities` table (System:Datagrok). It
// always exists and exposes an `id` UUID column used by the column
// transformation. Measure base column count dynamically (init_db's `entities`
// has 9 columns today, but additive db_up migrations could change that —
// don't hardcode the baseline).
const SQL_PRODUCTS = 'select * from entities';
const NEW_COLUMN_EXPRESSION = '${id}';

test.describe.serial(`Query transformations (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
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

  test('Add a ${productid} column transformation, save, run, verify persistence, then delete and verify gone', async ({ page }) => {
    // This test exercises the full transformation lifecycle — create → transform → run → reopen → delete → reopen.
    // The default 60s budget isn't enough; give it more room.
    test.setTimeout(180_000);
    // --- 1. Set up the test query (acts as stand-in for the shared `Products` query from the Test Track). ---
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await rightClickTreeNode(page, `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}`);
    await clickMenuItemExact(page, 'New Query...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await setQueryName(page, QUERY_NAME);
    await setQuerySql(page, SQL_PRODUCTS);
    // Play populates the in-memory DataFrame that the Transformations tab operates on.
    // IMPORTANT: save() flushes that in-memory state, causing "Previous step was not completed"
    // when we later add a transformation step — so do Play → Transformations → Save, not save-in-the-middle.
    await runQueryViaPlay(page);

    // Capture the base column count of `select * from entities` so the
    // post-transform assertion is independent of incidental schema changes.
    const baseColumnCount = await page.evaluate(() =>
      (window as unknown as { grok: { shell: { tv: { dataFrame: { columns: { length: number } } } } } })
        .grok.shell.tv.dataFrame.columns.length);
    expect(baseColumnCount).toBeGreaterThan(0);

    // --- 2. Add an "Add New Column" transformation step with the ${id} expression. ---
    await openTransformationsTab(page);
    expect(await transformationStepNames(page)).toEqual([QUERY_NAME]);
    await addNewColumnTransformation(page, NEW_COLUMN_EXPRESSION);
    expect(await transformationStepNames(page)).toEqual([QUERY_NAME, 'Add New Column']);

    // --- 3. Save and verify the column count expands when running via Toolbox > Actions > Run query... ---
    await saveQuery(page, QUERY_NAME);
    await runQueryViaActions(page, QUERY_NAME);
    const columnCount = await page.evaluate(() =>
      (window as unknown as { grok: { shell: { tv: { dataFrame: { columns: { length: number } } } } } })
        .grok.shell.tv.dataFrame.columns.length);
    expect(columnCount).toBe(baseColumnCount + 1);

    // --- 4. Close everything, reopen the query, and verify the transformation persists. ---
    await page.evaluate(() => (window as unknown as { grok: { shell: { closeAll: () => void } } }).grok.shell.closeAll());
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    const queryNodeName = `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---${queryTreeNodeSuffix(QUERY_NAME)}`;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Edit...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await waitForQuerySql(page, SQL_PRODUCTS);
    await openTransformationsTab(page);
    expect(await transformationStepNames(page)).toEqual([QUERY_NAME, 'Add New Column']);

    // --- 5. Remove the transformation step, save, and verify it's gone after a fresh reopen. ---
    await deleteTransformationStep(page, 1);
    expect(await transformationStepNames(page)).toEqual([QUERY_NAME]);
    await saveQuery(page, QUERY_NAME);

    await page.evaluate(() => (window as unknown as { grok: { shell: { closeAll: () => void } } }).grok.shell.closeAll());
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Edit...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await waitForQuerySql(page, SQL_PRODUCTS);
    await openTransformationsTab(page);
    expect(await transformationStepNames(page)).toEqual([QUERY_NAME]);

    // Sanity: the query itself still exists on the server.
    expect(await findQueryByFriendlyName(page, QUERY_NAME)).not.toBeNull();
  });
});

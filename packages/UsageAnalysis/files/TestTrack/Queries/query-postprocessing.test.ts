import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  POSTGRES_CONNECTION,
  clickMenuItemExact,
  deleteQueryByFriendlyName,
  expandDbProvider,
  findQueryByFriendlyName,
  goHome,
  rightClickTreeNode,
  runQueryViaPlay,
  saveQuery,
  setQueryName,
  setQueryPostProcessScript,
  setQuerySql,
} from './helpers';

// Test Track scenario `query-postprocessing.md` (order 13).
//
// SCOPE — this autotest is intentionally a server-roundtrip smoke. The runtime assertion
// from the manual ("Run the query → green info balloon `77` appears") cannot be reached
// from DOM automation under any of the paths we tried, all documented in
// `query-postprocessing-ui.md`:
//
//   1. Adding scatterplot/correlation viewers to the Layout tab requires drag-drop or
//      a JS handle on the nested TableView — same Dart pointer-handler / canvas blocker
//      as `query-layout.test.ts`.
//
//   2. Typing the post-process body into the Post-Process tab's `ScriptView` editor:
//      the CodeMirror inside the embedded ScriptView only updates
//      `DataQueryView.postProcessScript` through `scriptView.onScriptChanged`
//      (`data_query_view.dart:850-861`); JS-injected values (`cm.setValue`,
//      `keyboard.insertText`) don't re-fire that listener under headless Playwright.
//      We bypass the editor and patch via dapi so the field still round-trips through
//      the REST API. UI typing into CM is regression-tested by `typeQuerySql` and
//      `addNewColumnTransformation`; the Tabs widget by `openTransformationsTab`.
//
//   3. Running the saved query end-to-end so the post-process emits the green info
//      balloon. Right-click → Run AND `q.prepare().call()` both raise red
//      "Handler for null is not registered" — the client function registry never
//      rebinds the handler for a query whose `postProcessScript` was patched in-session
//      (full `page.reload()` + 2 s settle does not help, neither does Edit + Play in
//      the editor — the latter only re-runs the SQL preview without applying the
//      post-process body). The right-click → Run UI is regression-tested in
//      `chembl-parameterized-and-project.test.ts` against a stable fixture; balloon-on-Run
//      is exercised by the manual checks in `query-postprocessing-ui.md`.
//
// COVERED HERE (autotest smoke):
//   * Create `Test_Postprocessing` query via UI (right-click → New Query, Name + SQL).
//   * Run via Play (so the editor's preview grid populates).
//   * Save via UI Save button.
//   * Patch `postProcessScript` via dapi (JS API fallback — see SCOPE 2).
//   * Read the saved query back from the server and assert `postProcessScript` contains
//     the expected snippet.

const PROVIDER = 'Postgres';
const QUERY_NAME = 'Test_Postprocessing';
const SQL_BODY = 'select * from products';
const POST_PROCESS_SCRIPT = `//name: postprocess
//tags: postprocessing
//input: dataframe result
//output: dataframe out

grok.shell.info(result.rowCount);
out = result;
`;

test.describe.serial(`Query post-processing (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
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

  test('Create query via UI, save, patch postProcessScript via dapi, verify server round-trip', async ({ page }) => {
    test.setTimeout(120_000);

    // 1-2. Create the query via right-click on the connection.
    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await rightClickTreeNode(page, `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}`);
    await clickMenuItemExact(page, 'New Query...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await setQueryName(page, QUERY_NAME);
    await setQuerySql(page, SQL_BODY);

    // 3. Run via Play — populates the editor preview.
    await runQueryViaPlay(page);

    // 6. Save the query via the UI Save button.
    await saveQuery(page, QUERY_NAME);

    // 4. (JS API fallback) Patch the post-process body. See SCOPE 2 above.
    await setQueryPostProcessScript(page, QUERY_NAME, POST_PROCESS_SCRIPT);

    // Server-side assertion: dapi save persisted the post-process body end-to-end.
    const persisted = await page.evaluate(async (fn) => {
      const grok = (window as unknown as { grok: any }).grok;
      const q = (await grok.dapi.queries.filter(`friendlyName = "${fn}"`).list())[0];
      return { has: !!q?.postProcessScript, content: q?.postProcessScript ?? '' };
    }, QUERY_NAME);
    expect(persisted.has, 'postProcessScript should be persisted on the server after dapi save').toBe(true);
    expect(persisted.content).toContain('grok.shell.info(result.rowCount)');

    // Sanity: the query exists on the server and is queryable by friendlyName.
    expect(await findQueryByFriendlyName(page, QUERY_NAME)).not.toBeNull();
  });
});

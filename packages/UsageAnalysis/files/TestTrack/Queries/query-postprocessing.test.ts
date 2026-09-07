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

    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await rightClickTreeNode(page, `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}`);
    await clickMenuItemExact(page, 'New Query...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await setQueryName(page, QUERY_NAME);
    await setQuerySql(page, SQL_BODY);

    await runQueryViaPlay(page);

    await saveQuery(page, QUERY_NAME);

    await setQueryPostProcessScript(page, QUERY_NAME, POST_PROCESS_SCRIPT);

    const persisted = await page.evaluate(async (fn) => {
      const grok = (window as unknown as { grok: any }).grok;
      const q = (await grok.dapi.queries.filter(`friendlyName = "${fn}"`).list())[0];
      return { has: !!q?.postProcessScript, content: q?.postProcessScript ?? '' };
    }, QUERY_NAME);
    expect(persisted.has, 'postProcessScript should be persisted on the server after dapi save').toBe(true);
    expect(persisted.content).toContain('grok.shell.info(result.rowCount)');

    expect(await findQueryByFriendlyName(page, QUERY_NAME)).not.toBeNull();
  });
});

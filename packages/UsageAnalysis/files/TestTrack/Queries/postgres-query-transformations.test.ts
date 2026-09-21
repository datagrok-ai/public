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

const PROVIDER = 'Postgres';
const QUERY_NAME = 'transform_test_query';
const SQL_PRODUCTS = 'select * from products';
const NEW_COLUMN_EXPRESSION = '${productid}';
const PRODUCTS_COLUMN_COUNT = 10; 

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

    test.setTimeout(180_000);

    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await rightClickTreeNode(page, `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}`);
    await clickMenuItemExact(page, 'New Query...');
    await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });
    await setQueryName(page, QUERY_NAME);
    await setQuerySql(page, SQL_PRODUCTS);

    await runQueryViaPlay(page);

    await openTransformationsTab(page);
    expect(await transformationStepNames(page)).toEqual([QUERY_NAME]);
    await addNewColumnTransformation(page, NEW_COLUMN_EXPRESSION);
    expect(await transformationStepNames(page)).toEqual([QUERY_NAME, 'Add New Column']);

    await saveQuery(page, QUERY_NAME);
    await runQueryViaActions(page, QUERY_NAME);
    const columnCount = await page.evaluate(() =>
      (window as unknown as { grok: { shell: { tv: { dataFrame: { columns: { length: number } } } } } })
        .grok.shell.tv.dataFrame.columns.length);
    expect(columnCount).toBe(PRODUCTS_COLUMN_COUNT + 1);

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

    expect(await findQueryByFriendlyName(page, QUERY_NAME)).not.toBeNull();
  });
});

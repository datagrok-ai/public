import { test, expect } from '@playwright/test';
import {
  POSTGRES_CONNECTION,
  clickMenuItemExact,
  expandDbConnection,
  expandDbProvider,
  expandDbSchemas,
  expandTreeNode,
  focusBrowseSidebar,
  getConnectionServerName,
  goHome,
  rightClickTreeNode,
  showContextPanel,
} from './helpers';

const PROVIDER = 'Postgres';
const SCHEMA = 'public';
const TABLE = 'orders';

test.describe.serial(`Table Get All / Get Top 100 (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
  test('Right-click orders table → Get All returns the full dataset, Get Top 100 returns exactly 100 rows', async ({ page }) => {
    test.setTimeout(120_000);

    await goHome(page);
    await showContextPanel(page);

    const connServerName = await getConnectionServerName(page, PROVIDER, POSTGRES_CONNECTION);
    expect(connServerName).toBeTruthy();

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await expandDbSchemas(page, PROVIDER, connServerName);
    const schemaNode = `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---Schemas---${SCHEMA}`;
    await expandTreeNode(page, schemaNode);

    const tableNode = `${schemaNode}---${TABLE}`;

    await rightClickTreeNode(page, tableNode);
    await clickMenuItemExact(page, 'Get All');
    await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 30_000 });
    const getAllRowCount = await page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.tv.dataFrame.rowCount);

    expect(getAllRowCount).toBeGreaterThan(100);

    await focusBrowseSidebar(page);

    await rightClickTreeNode(page, tableNode);
    await clickMenuItemExact(page, 'Get Top 100');

    await expect.poll(async () => page.evaluate(() => {
      const grok = (window as unknown as { grok: any }).grok;
      const tableViews = Array.from(grok.shell.views).filter((v: any) => v.type === 'TableView');
      if (tableViews.length < 2) return null;
      return (tableViews as any[]).map((v) => v.dataFrame?.rowCount ?? -1);
    }), { timeout: 30_000 })
      .toEqual(expect.arrayContaining([getAllRowCount, 100]));
  });
});

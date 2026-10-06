import { test, expect } from '@playwright/test';
import {
  PG_PASSWORD,
  applyAutomationSetup,
  closeMenuPopup,
  connectionNodeName,
  expandDbConnection,
  expandDbGroupWrapper,
  expandDbProvider,
  goHome,
  readMenuItems,
  rightClickTreeNode,
} from './helpers';

const PROVIDER = 'Postgres';
const CONNECTION = 'Datagrok';

const CONN_SERVER_NAME = 'Datagrok';
const CATALOG = 'datagrok';
const SCHEMA = 'public';

test.describe.serial('Connections / Schema (Postgres / Datagrok substitution)', () => {
  test('1. Connection context menu offers Browse / Browse queries', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await rightClickTreeNode(page, connectionNodeName(PROVIDER, CONNECTION));
    const items = await readMenuItems(page);
    expect(items.some((t) => /^(Browse schema|Browse)$/.test(t))).toBe(true);
    expect(items).toEqual(expect.arrayContaining(['Browse queries']));
    await closeMenuPopup(page);
  });

  test('2. Expand Catalogs → catalog → schema; first table exposes DB-table context menu', async ({ page }) => {
    test.skip(!PG_PASSWORD, 'DG_PG_PASSWORD not set — schema browsing requires live DB credentials');

    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);

    await expandDbGroupWrapper(page, PROVIDER, CONN_SERVER_NAME, 'Catalogs');

    const catalogsRoot = `tree-Databases---${PROVIDER}---${CONNECTION}---Catalogs`;
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      catalogsRoot,
      { timeout: 30_000 },
    );

    const clickExpanderByName = async (name: string) => page.evaluate((n) => {
      const el = document.querySelector(`[name="${n}"]`) as HTMLElement | null;
      if (el && !el.classList.contains('d4-tree-view-tri-expanded')) el.click();
    }, name);

    await clickExpanderByName(`tree-expander-Databases---${PROVIDER}---${CONNECTION}---Catalogs---${CATALOG}`);
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      `${catalogsRoot}---${CATALOG}`,
      { timeout: 30_000 },
    );
    await clickExpanderByName(`tree-expander-Databases---${PROVIDER}---${CONNECTION}---Catalogs---${CATALOG}---${SCHEMA}`);
    await page.waitForFunction(
      (prefix) => document.querySelectorAll(`[name^="${prefix}---"]`).length > 0,
      `${catalogsRoot}---${CATALOG}---${SCHEMA}`,
      { timeout: 30_000 },
    );

    const schemaPrefix = `${catalogsRoot}---${CATALOG}---${SCHEMA}---`;
    const firstTable = page.locator(`[name^="${schemaPrefix}"]`).first();
    await firstTable.waitFor({ state: 'visible', timeout: 20_000 });
    const firstTableName = (await firstTable.getAttribute('name'))!;

    await rightClickTreeNode(page, firstTableName);
    const tableMenu = await readMenuItems(page);
    expect(tableMenu).toEqual(expect.arrayContaining([
      'Get All', 'Get Top 100', 'New SQL Query...', 'New Visual Query...',
    ]));
    await closeMenuPopup(page);
  });
});

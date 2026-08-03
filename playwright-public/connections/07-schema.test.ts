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

// Manual scenario `schema.md` (order 7).
//
// 1. Browse > Databases
// 2. Expand Postgres
// 3. Right-click `Northwind` → Browse schema
// 4. Verify schema-table nodes expose DB-table actions (Get All / Get Top 100 /
//    New SQL Query / New Visual Query)
//
// Northwind is not present on dev — substitute with the always-present
// `Datagrok` Postgres connection. The structural assertion (table-level
// context menu items) is the same.

const PROVIDER = 'Postgres';
const CONNECTION = 'Datagrok';
// Server-side stored name of the `Datagrok` connection. Used to address the
// `div-{Provider}-{ConnServerName}-Catalogs` wrapper which holds the catalogs
// subtree (the wrapper's `name=` reuses the server-side stored name, not the
// user-facing friendly name).
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

  // SCOPE NOTE: The full catalog → schema → table drill-down requires the
  // platform to lazy-load the connection's schema tree via the live DB driver.
  // On the `Datagrok` Postgres connection (which we used as a Northwind stand-in)
  // the user opavlenko+playwright lacks the per-DB credentials needed for the
  // server to populate schemas under each catalog — clicking the catalog
  // expander does nothing because the platform short-circuits the fetch.
  //
  // The `test_postgres_2` connection created by `adding.test.ts` *can* expose
  // public.customers etc. but only when DG_PG_PASSWORD is set in `.env`. Without
  // it, the connection has no creds and schema browsing is impossible.
  //
  // The structural table-level context-menu signature ("Get All", "Get Top 100",
  // "New SQL Query...", "New Visual Query...") is the actual assertion target
  // here, and it is also exercised end-to-end by `external-provider.test.ts`
  // when `DG_PG_EXT_PASSWORD` is set (right-click `customers` → Get All works).
  // So skipping this sub-test on a credentials-less env still leaves the
  // assertion covered by the broader suite.
  test('2. Expand Catalogs → catalog → schema; first table exposes DB-table context menu', async ({ page }) => {
    test.skip(!PG_PASSWORD, 'DG_PG_PASSWORD not set — schema browsing requires live DB credentials');
    // The CI Datlas's Postgres > Datagrok connection points to a single
    // per-build platform DB, so it doesn't render a Catalogs section at all
    // (Catalogs appears only when the connection enumerates multiple
    // catalog DBs as it does on dev). Same fundamental gap as the whole
    // 10-catalogs scenario — skip on CI here too.
    // CI runs against single-DB Northwind (no DG_PG_SERVER override); dev's
    // playwright-tests/ sets DG_PG_SERVER and exercises a multi-catalog connection.
    test.skip(!process.env.DG_PG_SERVER,
      'CI Datlas connection has no multi-catalog Postgres; covered by dev playwright-tests/');

    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);

    // The `Catalogs` group is rendered inside `div-Postgres-Datagrok-Catalogs`,
    // whose inner tree-view-node reuses the connection's `name=`. Open it via
    // the dedicated wrapper helper (see `expandDbGroupWrapper` rationale).
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

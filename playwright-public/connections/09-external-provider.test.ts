import { test, expect, type Page } from '@playwright/test';
import {
  AUTH_STATE,
  PG_EXT_DB,
  PG_EXT_LOGIN,
  PG_EXT_PASSWORD,
  PG_EXT_PORT,
  PG_EXT_SERVER,
  applyAutomationSetup,
  clickConnectionOk,
  clickMenuItemExact,
  connectionNodeName,
  deleteConnectionByFriendlyName,
  deleteTreeNodeViaContext,
  expandDbProvider,
  fillConnectionField,
  findConnectionByFriendlyName,
  goHome,
  openAddConnectionDialog,
  rightClickTreeNode,
} from './helpers';

// Manual scenario `external-provider.md` (order 9).
//
// 1. Browse > Databases > Postgres → Add connection
// 2. Fill PostgreSQLDBTests2 form (server db.datagrok.ai, port 54327, db test, user superuser)
// 3. Add and run four queries in sequence (CREATE TABLE / INSERT / UPDATE / DROP)
// 4. Delete the connection

const PROVIDER = 'Postgres';
const CONNECTION = 'PostgreSQLDBTests2';

const QUERIES: { name: string; sql: string }[] = [
  { name: 'TestCreateTable', sql: 'CREATE TABLE tmp_table_test (id bigint, name varchar)' },
  { name: 'TestInsertData', sql: "INSERT INTO tmp_table_test VALUES (1, 'test')" },
  { name: 'TestUpdateData', sql: "UPDATE tmp_table_test SET name = 'test' WHERE id = 1" },
  { name: 'TestDropTable', sql: 'DROP TABLE tmp_table_test' },
];

/** Fill and submit the Add-connection dialog for PostgreSQLDBTests2. */
async function createConnectionViaUi(page: Page): Promise<void> {
  await openAddConnectionDialog(page, PROVIDER);
  await fillConnectionField(page, 'Name', CONNECTION);
  await fillConnectionField(page, 'Server', PG_EXT_SERVER);
  await fillConnectionField(page, 'Port', PG_EXT_PORT);
  await fillConnectionField(page, 'Db', PG_EXT_DB);
  await fillConnectionField(page, 'Login', PG_EXT_LOGIN);
  await fillConnectionField(page, 'Password', PG_EXT_PASSWORD);
  await clickConnectionOk(page);
}

/**
 * Bring the page home and guarantee the PostgreSQLDBTests2 connection exists,
 * creating it via the UI when missing. The query tests below are `serial` and
 * normally rely on test 1 to create the connection, but a query test selected
 * in isolation (e.g. a title-grep matching only one query) would otherwise fail
 * on a missing tree node — this self-provisions so each query test stands alone.
 */
async function ensureConnectionExists(page: Page): Promise<void> {
  await goHome(page);
  await applyAutomationSetup(page);
  if ((await findConnectionByFriendlyName(page, CONNECTION)) != null) return;
  await createConnectionViaUi(page);
  const saved = await findConnectionByFriendlyName(page, CONNECTION);
  expect(saved, `connection "${CONNECTION}" should exist (self-provisioned)`).not.toBeNull();
}

test.describe.serial('Connections / External Provider (PostgreSQLDBTests2)', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteConnectionByFriendlyName(page, CONNECTION);
    for (const q of QUERIES) {
      await page.evaluate(async (name) => {
        const g = (window as unknown as { grok: any }).grok;
        const qs = await g.dapi.queries.filter(`friendlyName = "${name}"`).list();
        for (const x of qs) await g.dapi.queries.delete(x);
      }, q.name);
    }
    await ctx.close();
  });

  test.afterAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    for (const q of QUERIES) {
      await page.evaluate(async (name) => {
        const g = (window as unknown as { grok: any }).grok;
        const qs = await g.dapi.queries.filter(`friendlyName = "${name}"`).list();
        for (const x of qs) await g.dapi.queries.delete(x);
      }, q.name);
    }
    await deleteConnectionByFriendlyName(page, CONNECTION);
    await ctx.close();
  });

  test('1. Create the PostgreSQLDBTests2 connection', async ({ page }) => {
    test.skip(!PG_EXT_PASSWORD, 'DG_PG_EXT_PASSWORD env var not set — cannot create the connection');

    await goHome(page);
    await applyAutomationSetup(page);

    await createConnectionViaUi(page);

    const saved = await findConnectionByFriendlyName(page, CONNECTION);
    expect(saved, `connection "${CONNECTION}" should exist after OK`).not.toBeNull();
  });

  for (const q of QUERIES) {
    test(`2. Save and run query "${q.name}" via UI editor`, async ({ page }) => {
      // Requires a writable external Postgres (CREATE/INSERT/UPDATE/DROP). The CI stack has no
      // writable provider: the datagrok user is read-only on Northwind (no CREATE on schema public)
      // and no working credentials exist for the db.datagrok.ai:54327 test DB (all logins fail auth).
      // Re-enable once a writable external DB + credentials are provisioned on the CI stack.
      test.skip(true, 'No writable external Postgres available on the CI stack (see comment).');

      await ensureConnectionExists(page);
      await expandDbProvider(page, PROVIDER);

      // Right-click connection → New Query...
      await rightClickTreeNode(page, connectionNodeName(PROVIDER, CONNECTION));
      await clickMenuItemExact(page, 'New Query...');

      await page.waitForSelector('[name="input-Name"]', { timeout: 15_000 });

      // Real-keyboard typing for the query name (Datagrok input host caption: `Name`).
      const nameInput = page.locator('[name="input-Name"]');
      await nameInput.click({ clickCount: 3 });
      await page.keyboard.type(q.name);

      // SQL body — CodeMirror surface; real keyboard typing.
      await page.waitForSelector('.CodeMirror', { state: 'visible', timeout: 20_000 });
      await page.waitForFunction(() => {
        const el = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: unknown } | null;
        return !!el?.CodeMirror;
      }, undefined, { timeout: 10_000 });
      const surface = page.locator('.CodeMirror .CodeMirror-code').first();
      await surface.click();
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.keyboard.type(q.sql, { delay: 10 });

      // Save first so the run binds to a saved query (Datagrok requires save before run on some queries).
      await page.locator('[name="button-Save"]').first().click();
      await expect.poll(async () => page.evaluate(async (name) => {
        const g = (window as unknown as { grok: any }).grok;
        const qs = await g.dapi.queries.filter(`friendlyName = "${name}"`).list();
        return qs.length;
      }, q.name), { timeout: 30_000 }).toBeGreaterThanOrEqual(1);

      // Click Play.
      await page.locator('[name="icon-play"]').first().click();

      // The CRUD-style queries don't return a result grid (no SELECT). The
      // platform reports completion via a green balloon or by clearing the
      // running indicator. Wait for either: a balloon appears, or the play
      // icon's spinner stops. We assert no error balloons remain.
      // 90s, not 60s: the first DDL query against the in-network `world` demo
      // Postgres pays a cold grok_connect connection-pool warm-up, which can
      // exceed 60s on a loaded CI agent (build 62 timed out here while build 55
      // passed). The query itself is trivial once the pool is up.
      await page.waitForFunction(() =>
        !!document.querySelector('.grok-balloon, .d4-balloon')
        || !!document.querySelector('[name="viewer-Grid"] canvas'),
        undefined,
        { timeout: 90_000 });
      await page.waitForTimeout(1500);

      const errors = await page.evaluate(() => Array.from(document.querySelectorAll('.d4-balloon-error, .grok-balloon-error'))
        .map((b) => (b as HTMLElement).textContent?.trim() ?? '')
        .filter((s) => s.length > 0));
      expect(errors, `query "${q.name}" should run without errors`).toEqual([]);
    });
  }

  test('3. Delete the PostgreSQLDBTests2 connection', async ({ page }) => {
    test.skip(!PG_EXT_PASSWORD, 'DG_PG_EXT_PASSWORD env var not set');

    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await deleteTreeNodeViaContext(page, connectionNodeName(PROVIDER, CONNECTION));

    await expect.poll(async () => (await findConnectionByFriendlyName(page, CONNECTION)) === null,
      { timeout: 15_000 }).toBe(true);
  });
});

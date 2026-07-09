import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  PG_PASSWORD,
  applyAutomationSetup,
  clickMenuItemExact,
  expandDbConnection,
  expandDbGroupWrapper,
  expandDbProvider,
  expandTreeNode,
  fillConnectionField,
  findConnectionByFriendlyName,
  goHome,
  rightClickTreeNode,
  selectConnectionField,
  showContextPanel,
} from './helpers';

// The Configure Identifiers dialog reads the live DB schema at open time. Without
// real Postgres credentials the platform can't enumerate tables/columns and the
// dialog never paints — see existing `identifiers-spec.ts` (it skips on the same
// condition). Set `DG_PG_PASSWORD` in `.env` to enable the full pass.
test.skip(!PG_PASSWORD, 'DG_PG_PASSWORD not set — scenario requires reachable Northwind DB to load schemas');

// Manual scenario `identifiers.md` (order 2 — runs after `adding`).
//
// 1. Right-click `test_postgres` → Configure Identifiers...
// 2. Schema = public, OK
// 3. Add identifier (CUSTOMER_ID / customers / customerid / [A-Z]{5}), SAVE
// 4. Reload
// 5. Open `customers` table; values in `customerid` should be blue (visual — see -ui.md)
// 6. Click column header → Context Panel > Details shows semantic type
// 7. Remove identifiers config from `test_postgres`
// 8. Reload, verify column no longer carries the semantic type

const PROVIDER = 'Postgres';
const CONNECTION = 'test_postgres';
// Server-side stored name (PascalCase'd from the friendly name) — used by the
// platform when naming the inline Schemas/Catalogs wrapper under the connection.
const CONN_SERVER_NAME = 'TestPostgres';
const CONNECTION_DASH = CONNECTION.replace(/_/g, '-');
const SCHEMA = 'public';
const TABLE = 'customers';
const COLUMN = 'customerid';
const SEM_TYPE = 'CUSTOMER_ID';
const REGEXP = '[A-Z]{5}';

const SCHEMA_NODE = `tree-Databases---${PROVIDER}---${CONNECTION_DASH}---Schemas---${SCHEMA}`;
const TABLE_NODE = `${SCHEMA_NODE}---${TABLE}`;
const ID_PUBLIC_NODE = `tree-Databases---${PROVIDER}---${CONNECTION_DASH}---Identifiers---${SCHEMA}`;

test.describe.serial('Connections / Identifiers', () => {
  test.beforeAll(async ({ browser }) => {
    // The identifiers scenario depends on `adding` having created `test_postgres`.
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    const conn = await findConnectionByFriendlyName(page, CONNECTION);
    if (!conn)
      throw new Error(`prerequisite: connection "${CONNECTION}" must exist (run adding.test.ts first)`);
    await ctx.close();
  });

  test('1. Configure identifier on test_postgres / public.customers.customerid', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    // Right-click the connection node → Configure Identifiers...
    const nodeName = `tree-Databases---${PROVIDER}---${CONNECTION.replace(/_/g, '-')}`;
    await rightClickTreeNode(page, nodeName);
    await clickMenuItemExact(page, 'Configure Identifiers...');

    // Stage 1: schema picker dialog. `Schema` renders as a `<select>` (combobox),
    // not a text input — `fillConnectionField` would miss the inner `input` and
    // never resolve. The OK button then advances to the per-connection
    // Identifiers view (see Stage 2).
    await page.locator('.d4-dialog').waitFor({ timeout: 15_000 });
    await selectConnectionField(page, 'Schema', SCHEMA);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 15_000 });

    // Stage 2: the platform opens a full-page "TestPostgres.public Identifiers"
    // view with an accordion (Identifiers / Joins / ... / Renderers) and global
    // IMPORT JSON / EXPORT JSON / SAVE buttons at the bottom. There is no second
    // dialog. Click the inline "Add a new identifier" affordance inside the
    // Identifiers section to spawn the per-row editor dialog.
    const addRowBtn = page.locator('[aria-label="Add a new identifier"]').first();
    await addRowBtn.waitFor({ state: 'visible', timeout: 30_000 });
    await addRowBtn.click();

    // Stage 3: per-row dialog (`name="dialog-Add-Identifier"`) with five fields.
    // Table and Column are `<select>`s; Semantic-Type and Match-Regexp are text
    // inputs. Column options refresh asynchronously after Table changes.
    const addDialog = page.locator('.d4-dialog[name="dialog-Add-Identifier"]');
    await addDialog.waitFor({ timeout: 15_000 });
    await fillConnectionField(page, 'Semantic-Type', SEM_TYPE);
    await selectConnectionField(page, 'Table', TABLE);
    await page.waitForFunction(
      (col) => {
        const sel = document.querySelector(
          '.d4-dialog[name="dialog-Add-Identifier"] [name="input-host-Column"] select',
        ) as HTMLSelectElement | null;
        return !!sel && Array.from(sel.options).some((o) => o.value === col);
      },
      COLUMN,
      { timeout: 10_000 },
    );
    await selectConnectionField(page, 'Column', COLUMN);
    await fillConnectionField(page, 'Match-Regexp', REGEXP);

    await addDialog.locator('[name="button-Add"]').click();
    await addDialog.waitFor({ state: 'detached', timeout: 10_000 });

    // Stage 4: click the global SAVE button on the Identifiers view (NOT a
    // dialog button). The platform confirms with a "Identifier configuration
    // saved..." balloon — but if a config already exists for this
    // connection/schema, an "Overwrite existing configuration?" dialog gates
    // the save first.
    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('button'))
        .find((b) => b.textContent?.trim().toLowerCase() === 'save' && !b.closest('.d4-dialog'));
      if (!btn) throw new Error('Global SAVE button not found on Identifiers view');
      (btn as HTMLElement).click();
    });
    // The overwrite dialog renders slightly after the SAVE click — poll for it
    // up to 5s with a real wait (Locator.isVisible does NOT wait, only checks
    // current state regardless of the timeout option).
    const overwriteOk = page.locator('.d4-dialog [name="button-OK"]').first();
    try {
      await overwriteOk.waitFor({ state: 'visible', timeout: 5_000 });
      await overwriteOk.click();
      await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 5_000 }).catch(() => null);
    }
    catch {
      // No overwrite dialog this run (first time saving for this connection/schema).
    }
    await page.waitForFunction(
      () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
        .some((b) => (b as HTMLElement).textContent?.includes('Identifier configuration saved')),
      undefined,
      { timeout: 15_000 },
    );
  });

  test('2. After reload, customerid carries the configured semantic type', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    // Open the customers table via the schema tree:
    // Postgres → test_postgres → Schemas → public → customers → Get All.
    // The Schemas wrapper's `name=` uses the *server-side* stored name
    // (`TestPostgres`), not the friendly name (`test_postgres`); the schema
    // tree nodes themselves use the dash-version of the friendly name.
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONN_SERVER_NAME, 'Schemas');
    await expandTreeNode(page, SCHEMA_NODE);

    await rightClickTreeNode(page, TABLE_NODE);
    await clickMenuItemExact(page, 'Get All');

    await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 60_000 });
    await showContextPanel(page);

    // Verify semType on the column server-side. This is the runtime answer to
    // "is the identifier active?" — UI then renders blue text from this value.
    // SCOPE NOTE: visual blue-highlight verification belongs in identifiers-ui.md;
    // grid cells are canvas-rendered, so we cannot read pixel colours from the DOM.
    const semType = await page.evaluate((col) => {
      const g = (window as unknown as { grok: any }).grok;
      const tv = g.shell.tv;
      if (!tv) return null;
      const c = tv.dataFrame.col(col);
      return c ? c.semType : null;
    }, COLUMN);
    expect(semType).toBe(SEM_TYPE);
  });

  test('3. Remove identifiers config and verify the column no longer carries the type', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    // The saved per-schema identifier configuration appears in the tree under
    // the connection as `Identifiers > public`. Right-clicking that node
    // exposes `Edit... / Remove...` — picking Remove drops the entire schema's
    // identifier config (the platform takes care of the persistence side).
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await rightClickTreeNode(page, ID_PUBLIC_NODE);
    await clickMenuItemExact(page, 'Remove...');

    // The platform shows a confirm dialog titled "Remove Identifier
    // Configuration". Scope to that title — `.d4-dialog .first()` would
    // otherwise pick up a stale, hidden dialog left from earlier steps.
    const confirmDlg = page.locator('.d4-dialog').filter({ hasText: 'Remove Identifier Configuration' });
    await confirmDlg.waitFor({ state: 'visible', timeout: 10_000 });
    await confirmDlg.locator('[name="button-OK"]').click();
    await confirmDlg.waitFor({ state: 'detached', timeout: 10_000 });
    // Wait for the platform's "Identifier configuration removed..." balloon —
    // confirms the Dapi round-trip completed before we reload.
    await page.waitForFunction(
      () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
        .some((b) => (b as HTMLElement).textContent?.includes('Identifier configuration removed')),
      undefined,
      { timeout: 15_000 },
    );

    // Reload and reopen the customers table to confirm semType cleared.
    await page.reload({ waitUntil: 'domcontentloaded' });
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONN_SERVER_NAME, 'Schemas');
    await expandTreeNode(page, SCHEMA_NODE);
    await rightClickTreeNode(page, TABLE_NODE);
    await clickMenuItemExact(page, 'Get All');
    await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 60_000 });

    const semType = await page.evaluate((col) => {
      const g = (window as unknown as { grok: any }).grok;
      const tv = g.shell.tv;
      if (!tv) return null;
      const c = tv.dataFrame.col(col);
      return c ? c.semType : null;
    }, COLUMN);
    // After config removal, the column should have no semantic type (empty/null/undefined).
    expect(semType ?? '').toBe('');
  });
});

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

test.skip(!PG_PASSWORD, 'DG_PG_PASSWORD not set — scenario requires reachable Northwind DB to load schemas');

const PROVIDER = 'Postgres';
const CONNECTION = 'test_postgres';

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

    const nodeName = `tree-Databases---${PROVIDER}---${CONNECTION.replace(/_/g, '-')}`;
    await rightClickTreeNode(page, nodeName);
    await clickMenuItemExact(page, 'Configure Identifiers...');

    await page.locator('.d4-dialog').waitFor({ timeout: 15_000 });
    await selectConnectionField(page, 'Schema', SCHEMA);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 15_000 });

    const addRowBtn = page.locator('[aria-label="Add a new identifier"]').first();
    await addRowBtn.waitFor({ state: 'visible', timeout: 30_000 });
    await addRowBtn.click();

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

    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('button'))
        .find((b) => b.textContent?.trim().toLowerCase() === 'save' && !b.closest('.d4-dialog'));
      if (!btn) throw new Error('Global SAVE button not found on Identifiers view');
      (btn as HTMLElement).click();
    });

    const overwriteOk = page.locator('.d4-dialog [name="button-OK"]').first();
    try {
      await overwriteOk.waitFor({ state: 'visible', timeout: 5_000 });
      await overwriteOk.click();
      await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 5_000 }).catch(() => null);
    }
    catch {

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

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await expandDbGroupWrapper(page, PROVIDER, CONN_SERVER_NAME, 'Schemas');
    await expandTreeNode(page, SCHEMA_NODE);

    await rightClickTreeNode(page, TABLE_NODE);
    await clickMenuItemExact(page, 'Get All');

    await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 60_000 });
    await showContextPanel(page);

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

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, CONNECTION);
    await rightClickTreeNode(page, ID_PUBLIC_NODE);
    await clickMenuItemExact(page, 'Remove...');

    const confirmDlg = page.locator('.d4-dialog').filter({ hasText: 'Remove Identifier Configuration' });
    await confirmDlg.waitFor({ state: 'visible', timeout: 10_000 });
    await confirmDlg.locator('[name="button-OK"]').click();
    await confirmDlg.waitFor({ state: 'detached', timeout: 10_000 });

    await page.waitForFunction(
      () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
        .some((b) => (b as HTMLElement).textContent?.includes('Identifier configuration removed')),
      undefined,
      { timeout: 15_000 },
    );

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

    expect(semType ?? '').toBe('');
  });
});

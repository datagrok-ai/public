import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  PG_DB,
  PG_LOGIN,
  PG_PASSWORD,
  PG_PORT,
  PG_SERVER,
  applyAutomationSetup,
  clickConnectionOk,
  clickConnectionTest,
  deleteConnectionByFriendlyName,
  fillConnectionField,
  findConnectionByFriendlyName,
  goHome,
  openAddConnectionDialog,
  readConnectionField,
} from './helpers';

const PROVIDER = 'Postgres';
const NAME_1 = 'test_postgres';
const NAME_2 = 'test_postgres_2';

test.describe.serial('Connections / Adding (Postgres)', () => {

  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteConnectionByFriendlyName(page, NAME_1);
    await deleteConnectionByFriendlyName(page, NAME_2);
    await ctx.close();
  });

  test('1. Add `test_postgres` — fill fields with real keyboard, TEST, OK; verify saved', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    await openAddConnectionDialog(page, PROVIDER);

    await fillConnectionField(page, 'Name', NAME_1);

    await fillConnectionField(page, 'Server', PG_SERVER);
    await fillConnectionField(page, 'Port', PG_PORT);
    await fillConnectionField(page, 'Db', PG_DB);
    await fillConnectionField(page, 'Login', PG_LOGIN);
    if (PG_PASSWORD) await fillConnectionField(page, 'Password', PG_PASSWORD);

    expect(await readConnectionField(page, 'Name')).toBe(NAME_1);
    expect(await readConnectionField(page, 'Server')).toBe(PG_SERVER);
    expect(await readConnectionField(page, 'Port')).toBe(PG_PORT);
    expect(await readConnectionField(page, 'Db')).toBe(PG_DB);
    expect(await readConnectionField(page, 'Login')).toBe(PG_LOGIN);

    await clickConnectionTest(page);

    await clickConnectionOk(page);

    const saved = await findConnectionByFriendlyName(page, NAME_1);
    expect(saved, `connection "${NAME_1}" should exist server-side after OK`).not.toBeNull();
    expect(saved!.friendlyName).toBe(NAME_1);
    expect(saved!.dataSource).toBe(PROVIDER);
  });

  test('2. Add `test_postgres_2` — same flow, second connection', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    await openAddConnectionDialog(page, PROVIDER);

    await fillConnectionField(page, 'Name', NAME_2);
    await fillConnectionField(page, 'Server', PG_SERVER);
    await fillConnectionField(page, 'Port', PG_PORT);
    await fillConnectionField(page, 'Db', PG_DB);
    await fillConnectionField(page, 'Login', PG_LOGIN);
    if (PG_PASSWORD) await fillConnectionField(page, 'Password', PG_PASSWORD);

    await clickConnectionOk(page);

    const saved = await findConnectionByFriendlyName(page, NAME_2);
    expect(saved, `connection "${NAME_2}" should exist server-side after OK`).not.toBeNull();
    expect(saved!.friendlyName).toBe(NAME_2);
  });

});

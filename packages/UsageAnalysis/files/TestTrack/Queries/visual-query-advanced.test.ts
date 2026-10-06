import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  POSTGRES_CONNECTION,
  clickMenuItemExact,
  deleteProjectByFriendlyName,
  expandDbConnection,
  expandDbProvider,
  findProjectByFriendlyName,
  findQueryByFriendlyName,
  goHome,
  rightClickTreeNode,
} from './helpers';

const PROVIDER = 'Postgres';
const FIXTURE_QUERY_FN = 'postgres customers in @country';
const FIXTURE_QUERY_NODE =
  'tree-Databases---Postgres---Northwind---postgres-customers-in-@country';
const PROJECT_NAME = 'test_visual_advanced_playwright';

test.describe.serial(`Visual query advanced runtime (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteProjectByFriendlyName(page, PROJECT_NAME);
    await ctx.close();
  });

  test.afterAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteProjectByFriendlyName(page, PROJECT_NAME);
    await ctx.close();
  });

  test('Run fixture query, parameter refresh, add result viewer, save project', async ({ page }) => {
    test.setTimeout(180_000);

    await goHome(page);
    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);

    expect(await findQueryByFriendlyName(page, FIXTURE_QUERY_FN),
      `fixture query "${FIXTURE_QUERY_FN}" must exist on public`).not.toBeNull();

    await rightClickTreeNode(page, FIXTURE_QUERY_NODE);
    await clickMenuItemExact(page, 'Run');

    const dialog = page.locator('.d4-dialog');
    await expect(dialog).toBeVisible({ timeout: 15_000 });
    const country = dialog.locator('input[name$="---Country"]');
    await expect(country).toHaveValue('France', { timeout: 10_000 });
    await dialog.locator('[name="button-OK"]').click();

    await expect.poll(async () => page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.tv?.dataFrame?.rowCount ?? 0),
    { timeout: 30_000 }).toBeGreaterThan(0);
    const franceRows = await page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.tv.dataFrame.rowCount);

    const paramInput = page.locator('[name="input-Country"]');
    await paramInput.waitFor({ timeout: 10_000 });
    await expect(paramInput).toHaveValue('France');
    await paramInput.click({ clickCount: 3 });
    await page.keyboard.type('USA');
    await expect(paramInput).toHaveValue('USA');
    await page.locator('[name="button-REFRESH"]').first().click();
    await expect.poll(async () => page.evaluate(() =>
      (window as unknown as { grok: any }).grok.shell.tv.dataFrame.rowCount),
    { timeout: 30_000 }).not.toBe(franceRows);

    await page.evaluate(() => {
      const grok = (window as unknown as { grok: { shell: { tv: { addViewer: (type: string) => void } } } }).grok;
      grok.shell.tv.addViewer('Bar chart');
    });
    await expect(page.locator('[name="viewer-Bar-chart"]').first()).toBeVisible({ timeout: 10_000 });

    await page.locator('[name="button-Save"]').first().click();
    const saveDialog = page.locator('.d4-dialog').filter({ hasText: 'Save project' });
    await expect(saveDialog).toBeVisible({ timeout: 10_000 });
    const projectNameInput = saveDialog.locator('input[type="text"]').first();
    await projectNameInput.click({ clickCount: 3 });
    await page.keyboard.type(PROJECT_NAME);
    await saveDialog.locator('[name="button-OK"]').click();

    const shareDialog = page.locator('.d4-dialog').filter({ hasText: `Share ${PROJECT_NAME}` });
    await expect(shareDialog).toBeVisible({ timeout: 10_000 });
    await shareDialog.locator('[name="button-CANCEL"]').click();
    await expect(page.locator('.d4-dialog')).toHaveCount(0, { timeout: 10_000 });

    await expect.poll(async () => findProjectByFriendlyName(page, PROJECT_NAME),
    { timeout: 15_000 }).not.toBeNull();

    const inMemoryViewers = await page.evaluate(() => {
      const grok = (window as unknown as { grok: any }).grok;
      const tv = grok.shell.tv;
      return tv && tv.viewers ? Array.from(tv.viewers).map((v: any) => v.type) : [];
    });
    expect(inMemoryViewers).toEqual(expect.arrayContaining(['Grid', 'Bar chart']));
  });
});

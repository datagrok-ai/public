import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  clickMenuItemExact,
  deleteProjectByFriendlyName,
  expandDbProvider,
  expandTreeNode,
  findProjectByFriendlyName,
  goHome,
  rightClickTreeNode,
} from './helpers';

const CHEMBL_PROVIDER = 'Postgres';
const CHEMBL_CONNECTION = 'CHEMBL';
const FRAC_NODE_NAME =
  'tree-Databases---Postgres---CHEMBL---Search-----FRAC-classification-with-substructure-search';
const FRAC_QUERY_FN = 'Search | FRAC classification with substructure search';
const LEVEL1_DEFAULT = 'STEROL BIOSYNTHESIS IN MEMBRANES';
const FRAC_RESULT_VIEW = 'FRAC classification with substructure search';

const PROJECT_NAME = 'test_project_frac_playwright';

test.describe.serial('CHEMBL FRAC — parameterized query and save project', () => {
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

  test('7a. Parameterized query — FRAC classification dialog has level defaults, OK runs and returns rows', async ({ page }) => {
    test.setTimeout(120_000);
    await goHome(page);
    await expandDbProvider(page, CHEMBL_PROVIDER);
    await expandTreeNode(page, `tree-Databases---${CHEMBL_PROVIDER}---${CHEMBL_CONNECTION}`);

    await expandTreeNode(page, `tree-Databases---${CHEMBL_PROVIDER}---${CHEMBL_CONNECTION}---Search-`);

    await rightClickTreeNode(page, FRAC_NODE_NAME);
    await clickMenuItemExact(page, 'Run');

    const dialog = page.locator('.d4-dialog');
    await expect(dialog).toBeVisible({ timeout: 15_000 });
    await expect(dialog.locator('.d4-dialog-title')).toHaveText(FRAC_QUERY_FN);

    const level1 = dialog.locator('select[name$="---Level1"]');
    const level2 = dialog.locator('select[name$="---Level2"]');
    const level3 = dialog.locator('select[name$="---Level3"]');
    const level4 = dialog.locator('select[name$="---Level4"]');
    await expect(level1).toBeVisible({ timeout: 60_000 });
    await expect(level1).toHaveValue(LEVEL1_DEFAULT, { timeout: 30_000 });
    await expect(level2).toHaveValue('');
    await expect(level3).toHaveValue('');
    await expect(level4).toHaveValue('');

    await level1.selectOption({ index: 1 });
    await expect(level2).toHaveValue('');
    await expect(level3).toHaveValue('');
    await expect(level4).toHaveValue('');

    await level1.selectOption(LEVEL1_DEFAULT);

    await page.locator('.d4-dialog [name="button-OK"]').click();

    await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 30_000 });
    const result = await page.evaluate(() => {
      const tv = (window as unknown as { grok: { shell: { tv: { name: string; dataFrame: { rowCount: number; columns: { length: number } } } } } })
        .grok.shell.tv;
      return { name: tv.name, rowCount: tv.dataFrame.rowCount, columnCount: tv.dataFrame.columns.length };
    });
    expect(result.name).toBe(FRAC_RESULT_VIEW);
    expect(result.rowCount).toBeGreaterThan(0);
    expect(result.columnCount).toBeGreaterThan(0);
  });

  test('7b. Save project — run FRAC, add Trellis plot, save, close, reopen, verify the viewer is restored', async ({ page }) => {
    test.setTimeout(180_000);
    await goHome(page);
    await expandDbProvider(page, CHEMBL_PROVIDER);
    await expandTreeNode(page, `tree-Databases---${CHEMBL_PROVIDER}---${CHEMBL_CONNECTION}`);
    await expandTreeNode(page, `tree-Databases---${CHEMBL_PROVIDER}---${CHEMBL_CONNECTION}---Search-`);

    await rightClickTreeNode(page, FRAC_NODE_NAME);
    await clickMenuItemExact(page, 'Run');
    const paramDialog = page.locator('.d4-dialog');
    await expect(paramDialog).toBeVisible({ timeout: 15_000 });

    await expect(paramDialog.locator('select[name$="---Level1"]')).toBeVisible({ timeout: 60_000 });
    await expect(paramDialog.locator('select[name$="---Level1"]')).toHaveValue(LEVEL1_DEFAULT, { timeout: 30_000 });
    await paramDialog.locator('[name="button-OK"]').click();
    await page.waitForSelector('[name="viewer-Grid"] canvas', { timeout: 30_000 });

    await page.evaluate(() => {
      const grok = (window as unknown as { grok: { shell: { tv: { addViewer: (type: string) => void } } } }).grok;
      grok.shell.tv.addViewer('Trellis plot');
    });
    await expect(page.locator('[name="viewer-Trellis-plot"]').first()).toBeVisible({ timeout: 10_000 });

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

    await expect.poll(async () => findProjectByFriendlyName(page, PROJECT_NAME), { timeout: 15_000 })
      .not.toBeNull();

    const inMemoryViewers = await page.evaluate(() => {
      const grok = (window as unknown as { grok: any }).grok;
      const tv = grok.shell.tv;
      return tv && tv.viewers ? Array.from(tv.viewers).map((v: any) => v.type) : [];
    });
    expect(inMemoryViewers).toEqual(expect.arrayContaining(['Grid', 'Trellis plot']));

  });
});

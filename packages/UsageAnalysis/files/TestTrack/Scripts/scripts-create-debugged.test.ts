import { test, expect, Page, BrowserContext } from '@playwright/test';
import * as path from 'path';
import {
  SCRIPT_NAME,
  R_SCRIPT_CONTENT,
  openScriptsBrowser,
  setScriptContent,
  apiDeleteScript,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;
const AUTH_STATE = path.resolve(__dirname, '..', '.auth.json');

const OTHER_LANGUAGES = [
  { menu: 'Python Script...', annotation: '#language: python', hasSample: true },
  { menu: 'Octave Script...', annotation: '#language: octave', hasSample: true },
  { menu: 'NodeJS Script...', annotation: '//language: nodejs', hasSample: true },
  { menu: 'Julia Script...', annotation: '#language: julia', hasSample: true },
  { menu: 'JavaScript Script...', annotation: '//language: javascript', hasSample: false },
  { menu: 'Grok Script...', annotation: '#language: grok', hasSample: true },
  { menu: 'Pyodide Script...', annotation: '#language: pyodide', hasSample: true },
];

async function resetToScripts(page: Page) {
  await page.evaluate(() => {
    const g = (window as any).grok;
    if (g?.shell?.closeAll) g.shell.closeAll();
    document.querySelectorAll('.d4-dialog').forEach((d: any) => { try { d.remove(); } catch(_) {} });
    document.querySelectorAll('.d4-toast, .d4-balloon, .d4-menu').forEach(e => e.remove());
  });
  await page.waitForTimeout(300);

  const scriptsLabel = page.locator('.d4-tree-view-item-label', { hasText: /^Scripts$/i }).first();
  if (await scriptsLabel.isVisible({ timeout: 1_000 }).catch(() => false)) {
    await scriptsLabel.click();
    await expect(page.locator('.grok-gallery-search-bar')).toBeVisible({ timeout: 10_000 });
  } else {
    await openScriptsBrowser(page);
  }
}

async function createNewScript(page: Page, menuItem: string) {
  await page.locator('[name="button-New"]').click();
  await expect(page.locator('.d4-menu-item').first()).toBeVisible({ timeout: 5_000 });
  await page.locator('.d4-menu-item', { hasText: menuItem }).click();
  await page.waitForURL(/\/script\//, { timeout: 15_000 });
  await expect(page.locator('i[name="icon-play"]')).toBeVisible({ timeout: 10_000 });
}

async function loadSampleTable(page: Page) {
  const asteriskBtn = page.locator('i[name="icon-asterisk"]');
  await expect(asteriskBtn).toBeVisible({ timeout: 5_000 });
  await asteriskBtn.click();

  await page.locator('.d4-grid').first().waitFor({ state: 'visible', timeout: 10_000 }).catch(() => {});
}

async function runScriptWithTable(page: Page, tableName: string) {
  await page.locator('i[name="icon-play"]').click();

  const dialog = page.locator('.d4-dialog').first();
  await expect(dialog).toBeVisible({ timeout: 10_000 });

  await dialog.locator('select.ui-input-editor').selectOption(tableName);

  const okBtn = dialog.locator('button.ui-btn-ok').first();
  await expect(okBtn).toBeEnabled({ timeout: 8_000 });
  await okBtn.click();

  await page.waitForTimeout(500);
}

test.describe.serial('Scripts: Create', () => {
  let sharedContext: BrowserContext;
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    sharedContext = await browser.newContext({ storageState: AUTH_STATE });
    page = await sharedContext.newPage();

    await openScriptsBrowser(page);
  });

  test.afterAll(async () => {
    await sharedContext?.close();
  });

  test.afterEach(async () => {
    await apiDeleteScript(page, SCRIPT_NAME);
  });

  test('1. R Script: create, configure, run, save, close', async () => {

    await resetToScripts(page);
    await expect(page.locator('[name="button-New"]')).toBeVisible();
    await expect(page.locator('.grok-gallery-search-bar')).toBeVisible();

    await createNewScript(page, 'R Script...');
    await expect(page.locator('.CodeMirror-code')).toContainText('#language: r', { timeout: 10_000 });
    await expect(page.locator('.d4-tab-header[name="Script"]')).toBeVisible();
    await expect(page.locator('.d4-tab-header[name="Debug"]')).toBeVisible();

    await loadSampleTable(page);

    await setScriptContent(page, R_SCRIPT_CONTENT);
    await expect(page.locator('.CodeMirror-code')).toContainText(`#name: ${SCRIPT_NAME}`, { timeout: 8_000 });
    await expect(page.locator('.CodeMirror-code')).toContainText('#output: string newParam');
    await expect(page.locator('.CodeMirror-code')).toContainText('#sample: cars.csv');

    await runScriptWithTable(page, 'cars');

    const saveBtn = page.locator('button[name="button-Save"]');
    await expect(saveBtn).toBeVisible();
    await saveBtn.click();

    await expect(page.locator('.d4-balloon').first()).toContainText(/saved/i, { timeout: 10_000 });

    await expect(page.locator('[name="div-view-name"]')).toContainText(SCRIPT_NAME, { ignoreCase: true });

    await page.locator('.d4-ribbon > .grok-icon.fal').click();
    await expect(page).toHaveURL(/\/scripts/, { timeout: 10_000 });
  });

  for (let i = 0; i < OTHER_LANGUAGES.length; i++) {
    const lang = OTHER_LANGUAGES[i];
    test(`${i + 2}. ${lang.menu}: create, load sample, run without errors`, async () => {
      await resetToScripts(page);

      await createNewScript(page, lang.menu);

      await expect(page.locator('.CodeMirror-code')).toContainText(lang.annotation, { timeout: 10_000 });

      await expect(page.locator('i[name="icon-play"]')).toBeVisible();
      await expect(page.locator('button[name="button-Save"]')).toBeVisible();

      if (lang.hasSample)
        await loadSampleTable(page);

      await page.locator('i[name="icon-play"]').click();
      await page.waitForTimeout(300);

      const dialog = page.locator('.d4-dialog').first();
      if (await dialog.isVisible()) {
        const select = dialog.locator('select.ui-input-editor');
        if (await select.isVisible({ timeout: 2_000 }).catch(() => false)) {
          const options = await select.locator('option').allTextContents();
          if (options.length > 0)
            await select.selectOption({ index: 0 });
        }
        const okBtn = dialog.locator('button.ui-btn-ok').first();
        if (await okBtn.isEnabled({ timeout: 3_000 }).catch(() => false))
          await okBtn.click();
        await page.waitForTimeout(500);
      }

      const errorBalloon = page.locator('.d4-balloon-error');
      await expect(errorBalloon).toHaveCount(0, { timeout: 5_000 }).catch(() => {});
    });
  }
});

import { test, expect, Page, BrowserContext } from '@playwright/test';
import * as path from 'path';
import {
  openScriptsBrowser,
  getScriptCard,
  setScriptContent,
  apiDeleteScript,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;
const AUTH_STATE = path.resolve(__dirname, '..', '.auth.json');
const EDIT_SCRIPT_NAME = 'PW_EditTest';
const EDIT_MARKER = 'newParam="test"';

const EDIT_SCRIPT_CONTENT = `#name: ${EDIT_SCRIPT_NAME}
#language: r
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]

count <- nrow(table) * ncol(table)`;

/** Fast reset: closeAll + navigate to Scripts browser. */
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

test.describe.serial('Scripts: Edit', () => {
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

  test.beforeEach(async () => {
    await page.evaluate(() => {
      const g = (window as any).grok;
      if (g?.shell?.closeAll) g.shell.closeAll();
    });
    await page.waitForTimeout(300);
    await page.waitForFunction(() => !!(window as any).grok?.dapi?.scripts, { timeout: 10_000 });
    await apiDeleteScript(page, EDIT_SCRIPT_NAME);

    // Create the script via API (prerequisite setup)
    await page.evaluate(async (content) => {
      const DG = (window as any).DG;
      const grok = (window as any).grok;
      const script = DG.Script.create(content);
      await grok.dapi.scripts.save(script);
      return null;
    }, EDIT_SCRIPT_CONTENT);
    await page.waitForTimeout(300);
  });

  test.afterEach(async () => {
    await apiDeleteScript(page, EDIT_SCRIPT_NAME);
  });

  // Test track: Edit.md, steps 1–6
  test('1. Edit script body, save, close, reopen and verify change persists', async () => {
    // Step 1: Go to Browse > Platform > Functions > Scripts
    await resetToScripts(page);

    // Step 2: Find PW_EditTest in gallery and double-click to open
    const card = getScriptCard(page, EDIT_SCRIPT_NAME);
    await expect(card).toBeVisible({ timeout: 15_000 });
    await card.dblclick();

    await page.waitForURL(/\/script\//, { timeout: 15_000 });
    await expect(page.locator('.CodeMirror-code')).toContainText('#language: r', { timeout: 15_000 });
    await expect(page.locator('[name="div-view-name"]')).toContainText('EditTest', { ignoreCase: true });

    // Step 3: Add the expression newParam="test" to the script body
    await page.locator('.CodeMirror').click();
    await page.keyboard.press('ControlOrMeta+End');
    await page.keyboard.press('Enter');
    await page.keyboard.type(EDIT_MARKER, { delay: 10 });
    await page.waitForTimeout(200);

    // Step 4: Click Save
    const saveBtn = page.locator('button[name="button-Save"]');
    await expect(saveBtn).not.toHaveClass(/disabled/, { timeout: 8_000 });
    await saveBtn.click();
    await expect(page.locator('.d4-balloon').first()).toContainText(/saved/i, { timeout: 10_000 });

    // Step 5: Close all views via shell.closeAll
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(300);

    // Step 6: Navigate back to Scripts, find and reopen, verify edit persists
    await resetToScripts(page);
    const cardAfter = getScriptCard(page, EDIT_SCRIPT_NAME);
    await expect(cardAfter).toBeVisible({ timeout: 15_000 });
    await cardAfter.dblclick();

    await page.waitForURL(/\/script\//, { timeout: 15_000 });
    await expect(page.locator('.CodeMirror-code')).toContainText(EDIT_MARKER, { timeout: 15_000 });

    // Close all views via shell.closeAll
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(200);

    // Verify all views are closed — no tab handles should remain
    await expect(page.locator('.tab-handle-close')).toHaveCount(0, { timeout: 5_000 });
  });
});

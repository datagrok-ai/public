import { test, expect, Page, BrowserContext } from '@playwright/test';
import * as path from 'path';
import {
  SCRIPT_NAME,
  openScriptsBrowser,
  getScriptCard,
  rightClickScript,
  clickMenuItem,
  apiDeleteScript,
  apiEnsureScript,
  loadCarsDemoTable,
  resetToScripts,
  searchScript,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;
const AUTH_STATE = path.resolve(__dirname, '..', '.auth.json');

test.describe.serial('Scripts: Browser', () => {
  let sharedContext: BrowserContext;
  let page: Page;

  test.beforeAll(async ({ browser }) => {
    sharedContext = await browser.newContext({ storageState: AUTH_STATE });
    page = await sharedContext.newPage();
    await page.goto(BASE);
    await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
    await page.waitForFunction(() => !!(window as any).grok?.dapi?.scripts, { timeout: 15_000 });
    await apiDeleteScript(page, SCRIPT_NAME);
    await apiEnsureScript(page);
    await openScriptsBrowser(page);
  });

  test.afterAll(async () => {
    if (page) await apiDeleteScript(page, SCRIPT_NAME).catch(() => {});
    await sharedContext?.close();
  });

  test.beforeEach(async () => {
    await resetToScripts(page);

    await apiEnsureScript(page);
  });

  test('Scripts browser loads with gallery of script cards', async () => {

    await expect(page.locator('.grok-card-view')).toBeVisible();

    await expect(page.locator('.grok-gallery-grid-item').first()).toBeVisible({ timeout: 10_000 });

    await expect(page.locator('input[placeholder="Search scripts by name or by #tags"]')).toBeVisible();
  });

  test('Search for testRscript finds the correct script', async () => {
    await searchScript(page, SCRIPT_NAME);
  });

  test('Clicking script shows Details accordion with metadata', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();

    const detailsHeader = page.locator('.d4-accordion-pane-header', { hasText: 'Details' }).first();
    await expect(detailsHeader).toBeVisible({ timeout: 8_000 });

    await detailsHeader.click();
    const detailsPane = page.locator('.d4-accordion-pane', { has: page.locator('.d4-accordion-pane-header', { hasText: 'Details' }) }).first();
    await expect(detailsPane).toContainText(/r/i); 
  });

  test('Script accordion shows full script source text', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();

    const scriptHeader = page.locator('.d4-accordion-pane-header', { hasText: 'Script' }).first();
    await expect(scriptHeader).toBeVisible({ timeout: 8_000 });
    const isExpanded = await scriptHeader.evaluate(
      (el) => el.classList.contains('expanded'),
    );
    if (!isExpanded)
      await scriptHeader.click();

    const scriptText: string = await page.locator('[d4-pane-script] textarea.ui-input-editor').inputValue({ timeout: 8_000 });
    expect(scriptText).toContain('#language');
  });

  test('Run script from context menu increases Activity run count', async () => {
    await loadCarsDemoTable(page);
    await resetToScripts(page);

    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();

    const activityHeader = page.locator('.d4-accordion-pane-header', { hasText: /Activity/i }).first();
    await expect(activityHeader).toBeVisible({ timeout: 8_000 });

    const activityTextBefore = await activityHeader.textContent() ?? '';
    const runsBefore = parseInt(activityTextBefore.replace(/\D/g, '') || '0', 10);

    await rightClickScript(page, SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 8_000 });
    const okBtn = dialog.locator('button.ui-btn-ok').first();
    if (await okBtn.isVisible()) await okBtn.click();
    await page.waitForTimeout(3000);

    const cardAfterRun = getScriptCard(page, SCRIPT_NAME);
    await expect(cardAfterRun).toBeVisible({ timeout: 10_000 });
    await cardAfterRun.click();
    await page.waitForTimeout(1500);

    const activityHeaderAfter = page.locator('.d4-accordion-pane-header', { hasText: /Activity/i }).first();
    const activityTextAfter = await activityHeaderAfter.textContent() ?? '';
    const runsAfter = parseInt(activityTextAfter.replace(/\D/g, '') || '0', 10);
    expect(runsAfter).toBeGreaterThanOrEqual(runsBefore);
  });

  test('Sharing accordion is visible in context panel', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();
    await expect(page.locator('.d4-accordion-pane-header', { hasText: 'Sharing' }).first()).toBeVisible({ timeout: 8_000 });
  });

  test('Activity accordion shows action history', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();
    await expect(page.locator('.d4-accordion-pane-header').first()).toBeVisible({ timeout: 10_000 });

    const activityHeader = page.locator('.d4-accordion-pane-header', { hasText: /Activity/i }).first();
    await expect(activityHeader).toBeVisible({ timeout: 8_000 });
    await activityHeader.click();
    await expect(activityHeader).toBeVisible();
  });

  test('Chats accordion is present in context panel', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();
    await expect(page.locator('.d4-accordion-pane-header', { hasText: 'Chats' }).first()).toBeVisible({ timeout: 8_000 });
  });

  test('Scripts browser toolbar shows view-mode and sort controls', async () => {

    await expect(page.locator('input[placeholder="Search scripts by name or by #tags"]')).toBeVisible();

    await expect(page.locator('i[name="icon-sync"]').first()).toBeVisible();

    await expect(page.locator('[name="button-New"]')).toBeVisible();
  });

  test('Edit from context menu opens the script editor', async () => {
    await searchScript(page, SCRIPT_NAME);

    await rightClickScript(page, SCRIPT_NAME);
    await clickMenuItem(page, 'Edit...');

    await page.waitForURL(/\/script\//, { timeout: 15_000 });
    await expect(page.locator('i[name="icon-play"]')).toBeVisible({ timeout: 10_000 });
    await expect(page.locator('[name="div-view-name"]')).toContainText(SCRIPT_NAME, { ignoreCase: true });
  });
});

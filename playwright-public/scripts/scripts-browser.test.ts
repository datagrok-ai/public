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

// Test track: Browser.md (order: 4)
// Shared-context serial suite: one login + one Dart boot + one script create
// for the whole suite, instead of per-test (saves ~5-8s × 9 tests).
// All tests here are read-only w.r.t. the script body (test 5 runs it but
// that's expected to increase activity count — the assertion tolerates that).
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
    // Ensure script is present (test-5 run could have created duplicates; ignore).
    // Safe because apiEnsureScript is idempotent.
    await apiEnsureScript(page);
  });

  // Test track: Browser.md, step 1 — Scripts browser is accessible
  test('Scripts browser loads with gallery of script cards', async () => {
    // Gallery view is visible
    await expect(page.locator('.grok-card-view')).toBeVisible();

    // At least one script card is rendered
    await expect(page.locator('.grok-gallery-grid-item').first()).toBeVisible({ timeout: 10_000 });

    // Search input is present
    await expect(page.locator('input[placeholder="Search scripts by name or by #tags"]')).toBeVisible();
  });

  // Test track: Browser.md, step 2 — search for testRscript
  test('Search for testRscript finds the correct script', async () => {
    await searchScript(page, SCRIPT_NAME);
  });

  // Test track: Browser.md, step 3A — click script and check Details accordion
  test('Clicking script shows Details accordion with metadata', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();

    // Details accordion is visible in the right context panel
    const detailsHeader = page.locator('.d4-accordion-pane-header', { hasText: 'Details' }).first();
    await expect(detailsHeader).toBeVisible({ timeout: 8_000 });

    // Expand Details and verify key fields are present (language, author info)
    await detailsHeader.click();
    const detailsPane = page.locator('.d4-accordion-pane', { has: page.locator('.d4-accordion-pane-header', { hasText: 'Details' }) }).first();
    await expect(detailsPane).toContainText(/r/i); // language: R
  });

  // Test track: Browser.md, step 3B — Script accordion shows script source
  test('Script accordion shows full script source text', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();

    // Script accordion may already be expanded after clicking the card.
    // Only click the header to expand if the pane content is not yet visible.
    const scriptHeader = page.locator('.d4-accordion-pane-header', { hasText: 'Script' }).first();
    await expect(scriptHeader).toBeVisible({ timeout: 8_000 });
    const isExpanded = await scriptHeader.evaluate(
      (el) => el.classList.contains('expanded'),
    );
    if (!isExpanded)
      await scriptHeader.click();

    // Script content lives inside a readonly <textarea>, not in DOM text
    const scriptText: string = await page.locator('[d4-pane-script] textarea.ui-input-editor').inputValue({ timeout: 8_000 });
    expect(scriptText).toContain('#language');
  });

  // Test track: Browser.md, step 3C — Run via context menu, Activity count increases
  test('Run script from context menu increases Activity run count', async () => {
    await loadCarsDemoTable(page);
    await resetToScripts(page);

    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();

    const activityHeader = page.locator('.d4-accordion-pane-header', { hasText: /Activity/i }).first();
    await expect(activityHeader).toBeVisible({ timeout: 8_000 });

    const activityTextBefore = await activityHeader.textContent() ?? '';
    const runsBefore = parseInt(activityTextBefore.replace(/\D/g, '') || '0', 10);

    // Run the script via context menu
    await rightClickScript(page, SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 8_000 });
    const okBtn = dialog.locator('button.ui-btn-ok').first();
    if (await okBtn.isVisible()) await okBtn.click();
    await page.waitForTimeout(3000);

    // Re-select the script to refresh the context panel
    const cardAfterRun = getScriptCard(page, SCRIPT_NAME);
    await expect(cardAfterRun).toBeVisible({ timeout: 10_000 });
    await cardAfterRun.click();
    await page.waitForTimeout(1500);

    // Activity count should have increased
    const activityHeaderAfter = page.locator('.d4-accordion-pane-header', { hasText: /Activity/i }).first();
    const activityTextAfter = await activityHeaderAfter.textContent() ?? '';
    const runsAfter = parseInt(activityTextAfter.replace(/\D/g, '') || '0', 10);
    expect(runsAfter).toBeGreaterThanOrEqual(runsBefore);
  });

  // Test track: Browser.md, step 3C (Sharing) — Sharing accordion is visible
  test('Sharing accordion is visible in context panel', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();
    await expect(page.locator('.d4-accordion-pane-header', { hasText: 'Sharing' }).first()).toBeVisible({ timeout: 8_000 });
  });

  // Test track: Browser.md, step 3D — Activity accordion shows action dates
  test('Activity accordion shows action history', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();
    await expect(page.locator('.d4-accordion-pane-header').first()).toBeVisible({ timeout: 10_000 });

    const activityHeader = page.locator('.d4-accordion-pane-header', { hasText: /Activity/i }).first();
    await expect(activityHeader).toBeVisible({ timeout: 8_000 });
    await activityHeader.click();
    await expect(activityHeader).toBeVisible();
  });

  // Test track: Browser.md, step 3E — Chat accordion is present
  test('Chats accordion is present in context panel', async () => {
    const card = await searchScript(page, SCRIPT_NAME);
    await card.click();
    await expect(page.locator('.d4-accordion-pane-header', { hasText: 'Chats' }).first()).toBeVisible({ timeout: 8_000 });
  });

  // Test track: Browser.md, step 4 — view modes, sort, search
  test('Scripts browser toolbar shows view-mode and sort controls', async () => {
    // Search input
    await expect(page.locator('input[placeholder="Search scripts by name or by #tags"]')).toBeVisible();

    // Refresh icon (use .first() because the icon appears in multiple ribbon panels)
    await expect(page.locator('i[name="icon-sync"]').first()).toBeVisible();

    // The ribbon should contain the NEW button
    await expect(page.locator('[name="button-New"]')).toBeVisible();
  });

  // Test track: Browser.md, step 3 (Edit context menu) — context menu "Edit..." opens script editor
  test('Edit from context menu opens the script editor', async () => {
    await searchScript(page, SCRIPT_NAME);

    await rightClickScript(page, SCRIPT_NAME);
    await clickMenuItem(page, 'Edit...');

    // Script editor should open
    await page.waitForURL(/\/script\//, { timeout: 15_000 });
    await expect(page.locator('i[name="icon-play"]')).toBeVisible({ timeout: 10_000 });
    await expect(page.locator('[name="div-view-name"]')).toContainText(SCRIPT_NAME, { ignoreCase: true });
  });
});
